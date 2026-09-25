"""Experiment 2: choose eps so point-to-plane approximates the winding number.

eps is swept as a multiple of the mean sample spacing h (the same kNN
spacing the winding number uses for A_i), so the chosen value transfers
across meshes and sampling densities. For every eps we reconstruct both
surfaces and measure how far apart they are:
  - symmetric mean surface distance between the two marching-cubes meshes
  - volumetric IoU of the inside regions (Phi < 0) on the grid
Clean point clouds for all meshes, plus noisy versions of the bunny.

Outputs: results/tables/epsilon_sweep.csv, results/figures/fig2_epsilon_sweep.png
Run with --plot-only to redraw the figure from the saved CSV without recomputing.
"""

import csv
import sys

import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from reconstruction import (MESH_NAMES, NOISE_SEED, NUM_RECON_POINTS, RECON_SEED, WINDING_NEIGHBORS,  # noqa: E402
                            PointToPlaneImplicit, WindingNumberImplicit, add_position_noise,
                            load_normalized_mesh, output_path, reconstruct_mesh, sample_point_cloud,
                            symmetric_mean_distance)

EPSILON_RATIOS = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
NOISY_CASES = [('bunny', 0.01), ('bunny', 0.02), ('bunny', 0.04)]


def sweep_case(mesh_name, noise_sigma):
    mesh = load_normalized_mesh(mesh_name)
    cloud = sample_point_cloud(mesh, NUM_RECON_POINTS, RECON_SEED)
    if noise_sigma > 0:
        cloud = add_position_noise(cloud, noise_sigma, NOISE_SEED)
    winding = WindingNumberImplicit(cloud, WINDING_NEIGHBORS)
    spacing = winding.mean_neighbor_distances.mean()
    winding_mesh, winding_grid = reconstruct_mesh(winding)
    winding_inside = winding_grid < 0

    rows = []
    for ratio in EPSILON_RATIOS:
        plane_mesh, plane_grid = reconstruct_mesh(PointToPlaneImplicit(cloud, ratio * spacing))
        plane_inside = plane_grid < 0
        iou = np.sum(plane_inside & winding_inside) / np.sum(plane_inside | winding_inside)
        distance = symmetric_mean_distance(plane_mesh, winding_mesh)
        rows.append(dict(mesh=mesh_name, noise_sigma=noise_sigma, spacing_h=spacing, epsilon_over_h=ratio,
                         epsilon=ratio * spacing, surface_distance=distance, inside_iou=iou))
        print(f'  {mesh_name:7s} noise {noise_sigma:<5} eps = {ratio:<4} h: distance to WN surface {distance:.4f}, IoU {iou:.4f}',
              flush=True)
    return rows


SWEEP_TABLE = ('tables', 'epsilon_sweep.csv')


def run_sweep():
    all_rows = []
    for mesh_name in MESH_NAMES:
        all_rows += sweep_case(mesh_name, 0.0)
    for mesh_name, noise_sigma in NOISY_CASES:
        all_rows += sweep_case(mesh_name, noise_sigma)
    with open(output_path(*SWEEP_TABLE), 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)
    return all_rows


def load_sweep():
    with open(output_path(*SWEEP_TABLE), newline='') as csv_file:
        all_rows = list(csv.DictReader(csv_file))
    for row in all_rows:
        for key in row:
            if key != 'mesh':
                row[key] = float(row[key])
    return all_rows


def case_rows(all_rows, mesh_name, noise_sigma):
    return [row for row in all_rows if row['mesh'] == mesh_name and row['noise_sigma'] == noise_sigma]


def plot_sweep(all_rows):
    figure, (clean_axes, noisy_axes) = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for axes, cases, title in [
            (clean_axes, [(name, 0.0) for name in MESH_NAMES], 'clean samples'),
            (noisy_axes, [('bunny', 0.0)] + NOISY_CASES, 'bunny with position noise')]:
        for mesh_name, noise_sigma in cases:
            rows = case_rows(all_rows, mesh_name, noise_sigma)
            label = mesh_name if axes is clean_axes else f'sigma = {noise_sigma:g}'
            axes.plot([row['epsilon_over_h'] for row in rows], [row['surface_distance'] for row in rows],
                      marker='o', markersize=3, label=label)
        axes.set_xscale('log')
        axes.set_xlabel('eps / h  (h = mean kNN spacing of the samples)')
        axes.set_title(title, fontsize=10)
        axes.grid(True, which='both', alpha=0.3)
        axes.legend(fontsize=8)
    clean_axes.set_ylabel('distance to winding-number surface')
    figure.tight_layout()
    figure.savefig(output_path('figures', 'fig2_epsilon_sweep.png'), dpi=150)
    plt.close(figure)


def print_best(all_rows):
    print('Best eps / h per case (smallest surface distance):')
    for mesh_name, noise_sigma in [(name, 0.0) for name in MESH_NAMES] + NOISY_CASES:
        best = min(case_rows(all_rows, mesh_name, noise_sigma), key=lambda row: row['surface_distance'])
        print(f'  {mesh_name:7s} noise {noise_sigma:<5g}: eps / h = {best["epsilon_over_h"]:<4g} '
              f'(eps = {best["epsilon"]:.4f}, distance {best["surface_distance"]:.4f})')


if __name__ == '__main__':
    sweep_rows = load_sweep() if '--plot-only' in sys.argv else run_sweep()
    plot_sweep(sweep_rows)
    print_best(sweep_rows)
