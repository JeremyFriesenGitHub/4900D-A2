"""Experiment 1: sanity checks that motivate two implementation choices.

(a) Area weights: sum_i A_i should match the surface area, otherwise w(x) is
    not ~1 inside and the 0.5 threshold is biased. Prints the ratio for k = 1..12.

(b) Point-to-plane weight: with the handout's phi(r) = 1 / (r^2 + eps^2),
    Phi stays negative on both sides of the surface, so marching cubes finds
    no surface. Squaring the weight fixes it. Plots Phi along a line that
    crosses the surface, and reports the sign of Phi over the whole grid.

Outputs: results/figures/fig1_weight_locality.png, results/tables/area_weight_check.csv
"""

import csv

import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from reconstruction import (HANDOUT_WEIGHT_POWER, LOCAL_WEIGHT_POWER, MESH_NAMES, NUM_RECON_POINTS,  # noqa: E402
                            RECON_SEED, TEST_SEED, WINDING_NEIGHBORS, PointToPlaneImplicit,
                            WindingNumberImplicit, evaluate_on_grid, knn_area_weights, load_normalized_mesh,
                            mesh_is_empty, output_path, reconstruct_mesh, sample_point_cloud)

PROFILE_MESH = 'bunny'
NEIGHBOR_COUNTS = [1, 2, 4, 6, 8, 12]


def check_area_weights():
    print('(a) sum(A_i) / true surface area')
    rows = []
    for mesh_name in MESH_NAMES:
        mesh = load_normalized_mesh(mesh_name)
        cloud = sample_point_cloud(mesh, NUM_RECON_POINTS, RECON_SEED)
        true_area = mesh.farea.sum()
        ratios = [knn_area_weights(cloud.point, k)[0].sum() / true_area for k in NEIGHBOR_COUNTS]
        rows.append([mesh_name] + [f'{ratio:.3f}' for ratio in ratios])
        print(f'  {mesh_name:8s}', '  '.join(f'k={k}: {ratio:.2f}' for k, ratio in zip(NEIGHBOR_COUNTS, ratios)))
    with open(output_path('tables', 'area_weight_check.csv'), 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['mesh'] + [f'k={k}' for k in NEIGHBOR_COUNTS])
        writer.writerows(rows)


def profile_across_surface():
    mesh = load_normalized_mesh(PROFILE_MESH)
    cloud = sample_point_cloud(mesh, NUM_RECON_POINTS, RECON_SEED)
    winding = WindingNumberImplicit(cloud, WINDING_NEIGHBORS)
    spacing = winding.mean_neighbor_distances.mean()

    # Cross the true surface at a point that is NOT one of the samples, along
    # the true surface normal (a point from an independent sampling).
    crossing = sample_point_cloud(mesh, 1, TEST_SEED)
    signed_offsets = np.linspace(-0.3, 0.3, 1201)
    line_points = crossing.point[0] + signed_offsets[:, None] * crossing.normal[0]

    variants = [
        ('point-to-plane, handout weight, eps = 0.1h', PointToPlaneImplicit(cloud, 0.1 * spacing, HANDOUT_WEIGHT_POWER), 'tab:blue', '--'),
        ('point-to-plane, handout weight, eps = 1h', PointToPlaneImplicit(cloud, 1.0 * spacing, HANDOUT_WEIGHT_POWER), 'tab:blue', ':'),
        ('point-to-plane, squared weight, eps = 0.1h', PointToPlaneImplicit(cloud, 0.1 * spacing, LOCAL_WEIGHT_POWER), 'tab:green', '-'),
        ('point-to-plane, squared weight, eps = 1h', PointToPlaneImplicit(cloud, 1.0 * spacing, LOCAL_WEIGHT_POWER), 'tab:green', '-.'),
        ('winding number, 0.5 - w', winding, 'tab:orange', '-'),
    ]

    print(f'\n(b) {PROFILE_MESH}, {NUM_RECON_POINTS} samples, mean spacing h = {spacing:.4f}')
    figure, axes = plt.subplots(figsize=(7.5, 4.2))
    axes.plot(signed_offsets, signed_offsets, color='gray', linewidth=1, label='true signed distance')
    for label, implicit, color, style in variants:
        axes.plot(signed_offsets, implicit.evaluate_many(line_points), color=color, linestyle=style, label=label)
        grid_values = evaluate_on_grid(implicit)[0]
        reconstruction = reconstruct_mesh(implicit)[0]
        faces = 0 if mesh_is_empty(reconstruction) else reconstruction.face.shape[0]
        print(f'  {label:45s} grid corners with Phi > 0: {100 * np.mean(grid_values > 0):5.1f}%   '
              f'max Phi on grid: {grid_values.max():+.3f}   marching cubes faces: {faces}')
    axes.axhline(0.0, color='black', linewidth=0.6)
    axes.axvline(0.0, color='black', linewidth=0.6)
    axes.set_ylim(-0.6, 0.6)
    axes.set_xlabel('signed offset from the true surface along its normal (inside < 0)')
    axes.set_ylabel('Phi(x)')
    axes.set_title(f'Implicit functions across the {PROFILE_MESH} surface (h = mean sample spacing)', fontsize=10)
    axes.legend(fontsize=7, loc='upper left')
    figure.tight_layout()
    figure.savefig(output_path('figures', 'fig1_weight_locality.png'), dpi=150)
    plt.close(figure)


if __name__ == '__main__':
    check_area_weights()
    profile_across_surface()
