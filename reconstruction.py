"""Surface reconstruction from oriented point clouds with two implicit functions.

1. Blended point-to-plane distances (method from class):

       Phi(x) = sum_i phi(|x - p_i|) n_i^T (x - p_i) / sum_i phi(|x - p_i|)
       phi(r) = 1 / (r^2 + epsilon^2)

2. Generalized winding number of a point cloud (Barill et al., SIGGRAPH 2018):

       w(x)   = 1 / (4 pi) sum_i A_i n_i^T (p_i - x) / |p_i - x|^3
       Phi(x) = 0.5 - w(x)

   with A_i = (mean distance from p_i to its k nearest neighbours)^2.

Both functions are negative inside and positive outside, which is the sign
convention geomproc.marching_cubes expects (it produces outward normals).
"""

import math
import os
import random

import numpy as np
from scipy.spatial import cKDTree

import geomproc  # symlink to the geomproc/ package of a GeomProc clone (see README)

REPO_DIR = os.path.dirname(os.path.abspath(__file__))
# Root of that GeomProc clone, whose meshes/ folder holds the input meshes
GEOMPROC_DIR = os.path.dirname(os.path.dirname(os.path.realpath(geomproc.__file__)))

# Reconstruction volume used everywhere (meshes are normalized to [-1, 1]).
BOUNDS_MIN = -1.5
BOUNDS_MAX = 1.5
NUM_CUBES = 64

# Experiment settings shared by all exp*.py scripts
MESH_NAMES = ['sphere', 'cube', 'bunny', 'camel', 'kid']  # closed GeomProc meshes (bunny has small holes)
NUM_RECON_POINTS = 2000          # same as GeomProc's RBF reconstruction demo
TEST_POINTS_MULTIPLIER = 10      # test set is one order of magnitude larger
WINDING_NEIGHBORS = 6            # k for A_i; k = 6 makes sum(A_i) ~ surface area (exp1 prints the check)
RECON_SEED, TEST_SEED, NOISE_SEED = 1, 2, 3
RESULTS_DIR = os.path.join(REPO_DIR, 'results')


# ---------------------------------------------------------------------------
# Input data
# ---------------------------------------------------------------------------

def load_normalized_mesh(mesh_name):
    """Load a mesh by GeomProc name ('bunny') or by path, normalize it to
    [-1, 1] and compute its normals."""
    if os.path.exists(mesh_name):
        mesh_path = mesh_name
    else:
        mesh_path = os.path.join(GEOMPROC_DIR, 'meshes', mesh_name + '.obj')
    mesh = geomproc.load(mesh_path)
    mesh.normalize()
    mesh.compute_vertex_and_face_normals()
    return mesh


def sample_point_cloud(mesh, num_points, seed):
    """Uniform area-weighted samples with face normals (geomproc.mesh.sample).

    geomproc.mesh.sample draws from Python's `random` module, not numpy,
    so that is the generator we seed.
    """
    random.seed(seed)
    return mesh.sample(num_points)


def add_position_noise(point_cloud, noise_sigma, seed):
    """Return a copy of the point cloud with zero-mean Gaussian noise of
    standard deviation noise_sigma added to every coordinate. Normals are
    left untouched so the experiment isolates positional noise.

    (geomproc.mesh.add_noise is not used: it adds uniform noise in
    [0, scale), which also shifts the whole shape by scale / 2.)
    """
    noise_generator = np.random.default_rng(seed)
    noisy_cloud = point_cloud.copy()
    noisy_cloud.point = noisy_cloud.point + noise_generator.normal(
        scale=noise_sigma, size=noisy_cloud.point.shape)
    return noisy_cloud


def knn_area_weights(points, num_neighbors):
    """A_i = (mean distance from p_i to its k nearest neighbours)^2, with the
    neighbours found by GeomProc's KD-tree (geomproc.KDTree).

    Returns (area_weights, mean_neighbor_distances). The query asks for
    k + 1 neighbours because the closest "neighbour" of p_i is p_i itself.
    """
    points = np.asarray(points, dtype=float)
    point_list = points.tolist()
    tree = geomproc.KDTree(point_list, list(range(len(point_list))))
    # With an index list, nn_query returns [point, index] pairs, nearest first
    neighbor_indices = np.array([[index for _, index in tree.nn_query(point, num_neighbors + 1)]
                                 for point in point_list])
    neighbor_distances = np.linalg.norm(points[neighbor_indices[:, 1:]] - points[:, None, :], axis=2)
    mean_neighbor_distances = neighbor_distances.mean(axis=1)
    return mean_neighbor_distances ** 2, mean_neighbor_distances


# ---------------------------------------------------------------------------
# Implicit functions
# ---------------------------------------------------------------------------

class PointCloudImplicit:
    """Shared machinery: stores the oriented samples and evaluates the
    function on many query points at once, in chunks to bound memory."""

    chunk_size = 512

    def __init__(self, point_cloud):
        self.sample_points = np.asarray(point_cloud.point, dtype=float)
        sample_normals = np.asarray(point_cloud.normal, dtype=float)
        self.sample_normals = sample_normals / np.linalg.norm(sample_normals, axis=1, keepdims=True)
        # n_i^T p_i, reused to get n_i^T (x - p_i) = n_i^T x - n_i^T p_i
        self.sample_plane_offsets = np.sum(self.sample_points * self.sample_normals, axis=1)
        self.sample_squared_norms = np.sum(self.sample_points ** 2, axis=1)

    def pairwise_terms(self, query_points):
        """For query points x_q and samples p_i, return the matrices
        |x_q - p_i|^2 and n_i^T (x_q - p_i), both of shape (queries, samples).
        Uses |x - p|^2 = |x|^2 - 2 x^T p + |p|^2 so the work is two matrix
        products instead of a (queries, samples, 3) difference array."""
        query_squared_norms = np.sum(query_points ** 2, axis=1)
        squared_distances = (query_squared_norms[:, None]
                             - 2.0 * (query_points @ self.sample_points.T)
                             + self.sample_squared_norms[None, :])
        np.maximum(squared_distances, 0.0, out=squared_distances)
        signed_plane_distances = query_points @ self.sample_normals.T - self.sample_plane_offsets[None, :]
        return squared_distances, signed_plane_distances

    def evaluate_many(self, query_points):
        query_points = np.atleast_2d(np.asarray(query_points, dtype=float))
        values = np.empty(query_points.shape[0])
        for chunk_start in range(0, query_points.shape[0], self.chunk_size):
            chunk_stop = chunk_start + self.chunk_size
            values[chunk_start:chunk_stop] = self.evaluate_chunk(query_points[chunk_start:chunk_stop])
        return values

    def evaluate(self, query_point):
        """Single-point evaluation with the signature geomproc expects."""
        return float(self.evaluate_many(np.asarray(query_point, dtype=float)[None, :])[0])

    def evaluate_chunk(self, query_points):
        raise NotImplementedError


# Exponent of the point-to-plane weight phi(r) = 1 / (r^2 + eps^2)^power.
# power = 1 is the weight given in the handout. Summed over ALL samples of a
# surface it is not local enough: sum_i phi_i behaves like the integral of
# 1/r^2 over a 2D surface, which diverges logarithmically, so distant samples
# dominate the blend. Their plane distances average to about -3V/A < 0 for a
# closed surface, so Phi ends up negative on BOTH sides of the surface and
# marching cubes finds nothing (see exp1_weight_locality.py).
# power = 2 decays like 1/r^4, the weight integral converges, nearby samples
# dominate, and Phi becomes a proper signed distance near the surface.
# The other common fix, keeping power = 1 but summing only over the k nearest
# samples, is available through num_neighbors (exp4 reports it for reference).
# Neither makes Phi reliable FAR from the samples: see exp4 / fig3 (camel, kid).
HANDOUT_WEIGHT_POWER = 1
LOCAL_WEIGHT_POWER = 2
NEAREST_NEIGHBOR_BLEND = 16

# eps for point-to-plane, as a multiple of the mean sample spacing h. exp2
# sweeps eps / h: the distance between the point-to-plane and winding-number
# surfaces is flat below ~0.25h and grows after that (sphere, cube, bunny,
# noisy bunny; camel and kid are dominated by spurious surfaces), so we use:
MATCHED_EPSILON_RATIO = 0.1
# delta for the (optional) regularized winding-number kernel, also times h
REGULARIZATION_RATIO = 0.5


class PointToPlaneImplicit(PointCloudImplicit):
    """Phi(x) = sum_i phi_i n_i^T (x - p_i) / sum_i phi_i,
    phi_i = 1 / (|x - p_i|^2 + eps^2)^weight_power.

    The sum runs over all samples, or over the num_neighbors nearest
    samples of x when num_neighbors is given."""

    def __init__(self, point_cloud, epsilon, weight_power=LOCAL_WEIGHT_POWER, num_neighbors=None):
        super().__init__(point_cloud)
        self.epsilon = float(epsilon)
        self.weight_power = weight_power
        self.num_neighbors = num_neighbors
        self.sample_tree = cKDTree(self.sample_points) if num_neighbors else None

    def evaluate_chunk(self, query_points):
        if self.num_neighbors:
            return self.evaluate_chunk_nearest(query_points)
        squared_distances, signed_plane_distances = self.pairwise_terms(query_points)
        # phi_i, computed in place to avoid extra (queries x samples) temporaries
        weights = squared_distances
        weights += self.epsilon ** 2
        np.reciprocal(weights, out=weights)
        if self.weight_power != 1:
            np.power(weights, self.weight_power, out=weights)
        weighted_plane_distances = np.einsum('qi,qi->q', weights, signed_plane_distances)
        return weighted_plane_distances / np.sum(weights, axis=1)

    def evaluate_chunk_nearest(self, query_points):
        neighbor_distances, neighbor_indices = self.sample_tree.query(query_points, k=self.num_neighbors)
        offsets = query_points[:, None, :] - self.sample_points[neighbor_indices]
        signed_plane_distances = np.einsum('qkj,qkj->qk', offsets, self.sample_normals[neighbor_indices])
        weights = 1.0 / (neighbor_distances ** 2 + self.epsilon ** 2) ** self.weight_power
        return np.sum(weights * signed_plane_distances, axis=1) / np.sum(weights, axis=1)


class GaussianPointToPlaneImplicit(PointToPlaneImplicit):
    """Point-to-plane blend with the exponential weight phi(r) = exp(-r^2 / eps^2)
    (reference variant; several people in the course switched to it).

    Far from the samples every exp(-r_i^2 / eps^2) underflows to 0 and Phi
    becomes 0/0. Using exp(-(r_i^2 - r_min^2) / eps^2) instead multiplies all
    weights by the same factor, so Phi is unchanged, but the nearest sample
    always has weight 1."""

    def __init__(self, point_cloud, epsilon):
        super().__init__(point_cloud, epsilon)

    def evaluate_chunk(self, query_points):
        squared_distances, signed_plane_distances = self.pairwise_terms(query_points)
        squared_distances -= squared_distances.min(axis=1, keepdims=True)
        squared_distances *= -1.0 / self.epsilon ** 2
        weights = np.exp(squared_distances, out=squared_distances)
        return np.einsum('qi,qi->q', weights, signed_plane_distances) / np.sum(weights, axis=1)


class WindingNumberImplicit(PointCloudImplicit):
    """Phi(x) = 0.5 - w(x), w(x) = 1/(4 pi) sum_i A_i n_i^T (p_i - x) / |p_i - x|^3.

    kernel_regularization (delta) is an optional extension, not part of the
    handout: |p_i - x|^3 becomes (|p_i - x|^2 + delta^2)^(3/2), which removes
    the singularity at each sample (the same role eps plays in phi). The
    default delta = 0 is exactly the formula from the handout.
    """

    def __init__(self, point_cloud, num_neighbors=WINDING_NEIGHBORS, kernel_regularization=0.0):
        super().__init__(point_cloud)
        self.num_neighbors = num_neighbors
        self.kernel_regularization = float(kernel_regularization)
        self.area_weights, self.mean_neighbor_distances = knn_area_weights(self.sample_points, num_neighbors)

    def winding_number_chunk(self, query_points):
        squared_distances, signed_plane_distances = self.pairwise_terms(query_points)
        if self.kernel_regularization > 0.0:
            squared_distances += self.kernel_regularization ** 2
        cubed_distances = squared_distances * np.sqrt(squared_distances)
        # Floor avoids 0/0 if a query lands exactly on a sample (that term is 0).
        np.maximum(cubed_distances, 1e-300, out=cubed_distances)
        # n_i^T (p_i - x) / |p_i - x|^3, using n_i^T (p_i - x) = -n_i^T (x - p_i)
        dipole_terms = np.divide(signed_plane_distances, cubed_distances, out=cubed_distances)
        dipole_terms *= -1.0
        # sum_i A_i * term_i as a matrix-vector product
        return (dipole_terms @ self.area_weights) / (4.0 * math.pi)

    def winding_number_many(self, query_points):
        return 0.5 - self.evaluate_many(query_points)

    def evaluate_chunk(self, query_points):
        return 0.5 - self.winding_number_chunk(query_points)


# ---------------------------------------------------------------------------
# Marching cubes
# ---------------------------------------------------------------------------

def geomproc_axis(start_value, cube_size, num_cubes):
    """Grid coordinates as geomproc.marching_cubes generates them: num_cubes + 1
    values obtained by repeatedly adding cube_size in floating point."""
    coordinates = np.empty(num_cubes + 1)
    value = start_value
    for index in range(num_cubes + 1):
        coordinates[index] = value
        value += cube_size
    return coordinates


def evaluate_on_grid(implicit_function, num_cubes=NUM_CUBES, bounds_min=BOUNDS_MIN, bounds_max=BOUNDS_MAX):
    """Evaluate the implicit function on the (num_cubes + 1)^3 grid corners
    used by geomproc.marching_cubes. Returns (values[ix, iy, iz], start, cube_size)."""
    start = np.full(3, bounds_min, dtype=float)
    end = np.full(3, bounds_max, dtype=float)
    cube_size = (end - start) / num_cubes
    axis_coordinates = [geomproc_axis(start[axis], cube_size[axis], num_cubes) for axis in range(3)]
    grid_x, grid_y, grid_z = np.meshgrid(*axis_coordinates, indexing='ij')
    grid_points = np.column_stack([grid_x.ravel(), grid_y.ravel(), grid_z.ravel()])
    grid_values = implicit_function.evaluate_many(grid_points).reshape(grid_x.shape)
    return grid_values, start, cube_size


def reconstruct_mesh(implicit_function, num_cubes=NUM_CUBES, bounds_min=BOUNDS_MIN, bounds_max=BOUNDS_MAX):
    """Run geomproc.marching_cubes on an implicit function.

    geomproc.marching_cubes (precompute=True) calls fun(point) once per grid
    corner, i.e. (num_cubes + 1)^3 Python calls. We evaluate the whole grid
    in one vectorized pass first and hand GeomProc a function that looks
    the value up, so the output is identical but much faster.

    GeomProc versions before the Sept 24, 2026 fix controlled these loops with
    floating-point comparisons and could skip or repeat a grid layer for some
    grid sizes (never for the default 3/64 cells, which are exact in binary).
    Counting the lookups turns that into an error instead of a wrong mesh.
    Returns (mesh, grid_values).
    """
    grid_values, start, cube_size = evaluate_on_grid(implicit_function, num_cubes, bounds_min, bounds_max)
    num_lookups = 0

    def precomputed_value(point):
        nonlocal num_lookups
        num_lookups += 1
        grid_index = np.rint((point - start) / cube_size).astype(int)
        return grid_values[grid_index[0], grid_index[1], grid_index[2]]

    end = np.full(3, bounds_max, dtype=float)
    mesh = geomproc.marching_cubes(start, end, num_cubes, precomputed_value)
    if num_lookups != grid_values.size:
        raise RuntimeError(f'geomproc.marching_cubes evaluated {num_lookups} grid corners instead of '
                           f'{grid_values.size}. Update GeomProc (marching cubes fix of Sept 24, 2026) or '
                           'use a cell size that is exact in binary, such as 3/64.')
    if mesh.vertex.shape[0] > 0:
        mesh.compute_vertex_and_face_normals()
    return mesh, grid_values


# ---------------------------------------------------------------------------
# Evaluation helpers
# ---------------------------------------------------------------------------

def closest_points_on_triangles(query_points, corner_a, corner_b, corner_c):
    """Vectorized closest point on triangle (Ericson, Real-Time Collision
    Detection, sec. 5.1.5). All inputs have shape (m, 3), one triangle per
    query point. Degenerate triangles produce NaN, which callers skip."""

    def row_dot(first, second):
        return np.einsum('ij,ij->i', first, second)

    edge_ab = corner_b - corner_a
    edge_ac = corner_c - corner_a
    a_to_query = query_points - corner_a
    b_to_query = query_points - corner_b
    c_to_query = query_points - corner_c
    d1, d2 = row_dot(edge_ab, a_to_query), row_dot(edge_ac, a_to_query)
    d3, d4 = row_dot(edge_ab, b_to_query), row_dot(edge_ac, b_to_query)
    d5, d6 = row_dot(edge_ab, c_to_query), row_dot(edge_ac, c_to_query)
    area_c = d1 * d4 - d3 * d2
    area_b = d5 * d2 - d1 * d6
    area_a = d3 * d6 - d5 * d4

    with np.errstate(divide='ignore', invalid='ignore'):
        # Interior of the face (default), then edge and vertex regions. The
        # regions are disjoint; later assignments mirror the early returns of
        # the scalar algorithm.
        inverse_total = 1.0 / (area_a + area_b + area_c)
        closest = corner_a + edge_ab * (area_b * inverse_total)[:, None] + edge_ac * (area_c * inverse_total)[:, None]

        on_edge_bc = (area_a <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
        fraction_bc = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        closest[on_edge_bc] = (corner_b + (corner_c - corner_b) * fraction_bc[:, None])[on_edge_bc]

        on_edge_ac = (area_b <= 0) & (d2 >= 0) & (d6 <= 0)
        fraction_ac = d2 / (d2 - d6)
        closest[on_edge_ac] = (corner_a + edge_ac * fraction_ac[:, None])[on_edge_ac]

        on_vertex_c = (d6 >= 0) & (d5 <= d6)
        closest[on_vertex_c] = corner_c[on_vertex_c]

        on_edge_ab = (area_c <= 0) & (d1 >= 0) & (d3 <= 0)
        fraction_ab = d1 / (d1 - d3)
        closest[on_edge_ab] = (corner_a + edge_ab * fraction_ab[:, None])[on_edge_ab]

        on_vertex_b = (d3 >= 0) & (d4 <= d3)
        closest[on_vertex_b] = corner_b[on_vertex_b]

        on_vertex_a = (d1 <= 0) & (d2 <= 0)
        closest[on_vertex_a] = corner_a[on_vertex_a]
    return closest


def point_to_mesh_distances(query_points, mesh, num_candidate_triangles=32):
    """Exact unsigned distance from each query point to a triangle mesh,
    checking the triangles whose centroids are nearest to the query."""
    triangle_corners = mesh.vertex[mesh.face]  # (faces, 3 corners, 3 coords)
    centroid_tree = cKDTree(triangle_corners.mean(axis=1))
    num_candidates = min(num_candidate_triangles, mesh.face.shape[0])
    _, candidate_triangles = centroid_tree.query(query_points, k=num_candidates)
    candidate_triangles = candidate_triangles.reshape(query_points.shape[0], num_candidates)

    best_distances = np.full(query_points.shape[0], np.inf)
    for candidate_column in range(num_candidates):
        corners = triangle_corners[candidate_triangles[:, candidate_column]]
        closest = closest_points_on_triangles(query_points, corners[:, 0], corners[:, 1], corners[:, 2])
        distances = np.linalg.norm(query_points - closest, axis=1)
        distances[~np.isfinite(distances)] = np.inf
        best_distances = np.minimum(best_distances, distances)
    return best_distances


def symmetric_mean_distance(mesh_a, mesh_b):
    """Average of the mean vertex-to-surface distances A->B and B->A
    (a Chamfer-style distance between two reconstructions). Marching-cubes
    vertices are roughly uniform, so vertices stand in for surface samples."""
    if mesh_a.face.ndim != 2 or mesh_b.face.ndim != 2 or mesh_a.face.shape[0] == 0 or mesh_b.face.shape[0] == 0:
        return np.inf
    a_to_b = point_to_mesh_distances(mesh_a.vertex, mesh_b).mean()
    b_to_a = point_to_mesh_distances(mesh_b.vertex, mesh_a).mean()
    return 0.5 * (a_to_b + b_to_a)


def mesh_is_empty(mesh):
    return mesh.face.ndim != 2 or mesh.face.shape[0] == 0


def output_path(*parts):
    """Path under results/, creating folders as needed."""
    path = os.path.join(RESULTS_DIR, *parts)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path
