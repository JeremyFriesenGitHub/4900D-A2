"""Experiment 4: quantitative comparison on five meshes.

For each mesh, reconstruct from NUM_RECON_POINTS samples, then draw an
independent test set, 10x larger, from the same surface. The error of a test
point is |Phi(x)|, the distance to the surface according to the implicit
function (zero on the true surface). We report the max and mean over the
test set, plus the median, because the winding number's max comes from test
points that land right next to a sample, where its kernel is singular.

|0.5 - w| is not a length, so two unit-consistent extras are also reported,
both computed on the marching-cubes output:
  test -> reconstruction   distance from each test point to the reconstructed mesh
  reconstruction -> truth  distance from each reconstructed vertex to the true mesh
                           (this is the one that exposes spurious surfaces)

Methods on clean samples: winding number (handout formula), point-to-plane
with the squared weight and matched eps, and for reference the handout weight
over all samples, the handout weight over the 16 nearest samples, the
exponential weight exp(-r^2 / eps^2) with eps = h, and the regularized
winding number (extension). Noisy samples (noise on the
reconstruction samples only): the two main methods.

Outputs: results/tables/quantitative.csv and results/tables/quantitative.md
"""

import csv

import numpy as np

from reconstruction import (HANDOUT_WEIGHT_POWER, LOCAL_WEIGHT_POWER, MATCHED_EPSILON_RATIO, MESH_NAMES,  # noqa: E402
                            NEAREST_NEIGHBOR_BLEND, NOISE_SEED, NUM_RECON_POINTS, RECON_SEED,
                            REGULARIZATION_RATIO, TEST_POINTS_MULTIPLIER, TEST_SEED, WINDING_NEIGHBORS,
                            GaussianPointToPlaneImplicit, PointToPlaneImplicit, WindingNumberImplicit,
                            add_position_noise, load_normalized_mesh,
                            mesh_is_empty, output_path, point_to_mesh_distances, reconstruct_mesh,
                            sample_point_cloud)

NOISE_SIGMA = 0.02
MAIN_METHODS = ['winding number', 'point-to-plane']
REFERENCE_METHODS = ['point-to-plane, handout weight', f'point-to-plane, handout weight, {NEAREST_NEIGHBOR_BLEND} nearest',
                     'point-to-plane, exponential weight', 'winding number, regularized']


def build_methods(cloud, include_reference):
    winding = WindingNumberImplicit(cloud, WINDING_NEIGHBORS)
    spacing = winding.mean_neighbor_distances.mean()
    epsilon = MATCHED_EPSILON_RATIO * spacing
    methods = [
        ('winding number', winding, f'k = {WINDING_NEIGHBORS}'),
        ('point-to-plane', PointToPlaneImplicit(cloud, epsilon, LOCAL_WEIGHT_POWER),
         f'eps = {MATCHED_EPSILON_RATIO:g}h, squared weight'),
    ]
    if include_reference:
        methods += [
            (REFERENCE_METHODS[0], PointToPlaneImplicit(cloud, epsilon, HANDOUT_WEIGHT_POWER),
             f'eps = {MATCHED_EPSILON_RATIO:g}h, all samples'),
            (REFERENCE_METHODS[1], PointToPlaneImplicit(cloud, epsilon, HANDOUT_WEIGHT_POWER, NEAREST_NEIGHBOR_BLEND),
             f'eps = {MATCHED_EPSILON_RATIO:g}h'),
            (REFERENCE_METHODS[2], GaussianPointToPlaneImplicit(cloud, spacing), 'eps = 1h, phi = exp(-r^2/eps^2)'),
            (REFERENCE_METHODS[3], WindingNumberImplicit(cloud, WINDING_NEIGHBORS, REGULARIZATION_RATIO * spacing),
             f'delta = {REGULARIZATION_RATIO:g}h'),
        ]
    return spacing, methods


def main():
    rows = []
    for mesh_name in MESH_NAMES:
        mesh = load_normalized_mesh(mesh_name)
        clean_cloud = sample_point_cloud(mesh, NUM_RECON_POINTS, RECON_SEED)
        # Independent test set: different seed, 10x more points
        test_points = sample_point_cloud(mesh, TEST_POINTS_MULTIPLIER * NUM_RECON_POINTS, TEST_SEED).point
        for noise_sigma in [0.0, NOISE_SIGMA]:
            cloud = add_position_noise(clean_cloud, noise_sigma, NOISE_SEED) if noise_sigma > 0 else clean_cloud
            spacing, methods = build_methods(cloud, include_reference=(noise_sigma == 0.0))
            for method_name, implicit, parameters in methods:
                test_errors = np.abs(implicit.evaluate_many(test_points))
                reconstruction = reconstruct_mesh(implicit)[0]
                if mesh_is_empty(reconstruction):
                    test_to_recon = np.array([np.inf])
                    recon_to_truth = np.array([np.inf])
                else:
                    test_to_recon = point_to_mesh_distances(test_points, reconstruction)
                    recon_to_truth = point_to_mesh_distances(reconstruction.vertex, mesh)
                rows.append(dict(
                    mesh=mesh_name, noise_sigma=noise_sigma, method=method_name, parameters=parameters,
                    num_recon_points=cloud.point.shape[0], num_test_points=test_points.shape[0],
                    spacing_h=round(spacing, 5),
                    phi_error_max=test_errors.max(), phi_error_mean=test_errors.mean(),
                    phi_error_median=np.median(test_errors),
                    test_to_recon_max=test_to_recon.max(), test_to_recon_mean=test_to_recon.mean(),
                    recon_to_truth_max=recon_to_truth.max(), recon_to_truth_mean=recon_to_truth.mean()))
                print(f'{mesh_name:7s} noise {noise_sigma:<5} {method_name:45s} |Phi| max {test_errors.max():10.4f} '
                      f'mean {test_errors.mean():8.4f} median {np.median(test_errors):7.4f} | test->recon max '
                      f'{test_to_recon.max():.4f} mean {test_to_recon.mean():.4f} | recon->truth max '
                      f'{recon_to_truth.max():.4f} mean {recon_to_truth.mean():.4f}', flush=True)

    with open(output_path('tables', 'quantitative.csv'), 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    write_markdown(rows)


def format_value(value):
    return 'no surface' if not np.isfinite(value) else f'{value:.4f}'


def write_markdown(rows):
    def table(noise_sigma, methods, columns, headers):
        lines = ['| mesh | method | ' + ' | '.join(headers) + ' |',
                 '|---|---|' + '---:|' * len(headers)]
        for row in rows:
            if row['noise_sigma'] == noise_sigma and row['method'] in methods:
                lines.append(f"| {row['mesh']} | {row['method']} | "
                             + ' | '.join(format_value(row[column]) for column in columns) + ' |')
        return '\n'.join(lines)

    phi_columns = ['phi_error_max', 'phi_error_mean', 'phi_error_median']
    phi_headers = ['max \\|Phi\\|', 'mean \\|Phi\\|', 'median \\|Phi\\|']
    distance_columns = ['test_to_recon_max', 'test_to_recon_mean', 'recon_to_truth_max', 'recon_to_truth_mean']
    distance_headers = ['test->recon max', 'test->recon mean', 'recon->truth max', 'recon->truth mean']
    num_test = rows[0]['num_test_points']
    sections = [
        f'# Quantitative results\n\n{NUM_RECON_POINTS} reconstruction samples and {num_test} independent test '
        f'points per mesh; meshes normalized to [-1, 1]; marching cubes on a 64^3 grid over [-1.5, 1.5]^3.\n'
        f'Point-to-plane: eps = {MATCHED_EPSILON_RATIO:g} h (h = mean kNN spacing of the samples), squared weight. '
        f'Winding number: k = {WINDING_NEIGHBORS}.',
        '## 1. Error |Phi(x)| at the test points, clean samples (required comparison)\n\n'
        'Units differ: point-to-plane Phi is a length, 0.5 - w is a fraction of the full solid angle.\n\n'
        + table(0.0, MAIN_METHODS, phi_columns, phi_headers),
        '## 2. Distances measured on the reconstructed meshes, clean samples\n\n'
        'test->recon: test point to reconstructed surface. recon->truth: reconstructed vertex to true surface '
        '(large values = spurious surfaces).\n\n'
        + table(0.0, MAIN_METHODS, distance_columns, distance_headers),
        '## 3. Reference variants, clean samples\n\n'
        + table(0.0, REFERENCE_METHODS, phi_columns + distance_columns, phi_headers + distance_headers),
        f'## 4. Noisy samples (sigma = {NOISE_SIGMA}); test points stay on the clean surface\n\n'
        + table(NOISE_SIGMA, MAIN_METHODS, phi_columns + distance_columns, phi_headers + distance_headers),
    ]
    with open(output_path('tables', 'quantitative.md'), 'w') as markdown_file:
        markdown_file.write('\n\n'.join(sections) + '\n')


if __name__ == '__main__':
    main()
