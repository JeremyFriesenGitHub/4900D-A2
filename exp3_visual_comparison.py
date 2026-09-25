"""Experiment 3: visual comparison of the reconstructions (report figures).

fig3_clean_comparison.png   several meshes, clean samples:
    ground truth | samples | winding number | point-to-plane at three eps
fig4_noise_comparison.png   bunny with increasing position noise:
    noisy samples | winding number | point-to-plane at three eps
fig5_winding_regularization.png   (extension) why the winding-number surface is rough:
    raw kernel vs regularized kernel, next to point-to-plane at the matched eps

Every panel in a row uses the ground-truth extent as its window, so shapes
stay comparable and far-away spurious surfaces are clipped at the panel edge.
Every reconstruction is also saved as .obj under results/meshes/ (open in MeshLab).
"""

from reconstruction import (MATCHED_EPSILON_RATIO, NOISE_SEED, NUM_RECON_POINTS, RECON_SEED,  # noqa: E402
                            REGULARIZATION_RATIO, WINDING_NEIGHBORS, PointToPlaneImplicit, WindingNumberImplicit,
                            add_position_noise, load_normalized_mesh, mesh_is_empty, output_path, reconstruct_mesh,
                            sample_point_cloud)
from render import METHOD_COLORS, MESH_VIEWS, render_grid  # noqa: E402

CLEAN_MESHES = ['bunny', 'camel', 'kid', 'cube']
NOISE_LEVELS = [0.01, 0.02, 0.04]
EPSILON_RATIOS_SHOWN = [MATCHED_EPSILON_RATIO, 0.5, 2.0]


def reconstruct_and_save(implicit, file_stem):
    mesh = reconstruct_mesh(implicit)[0]
    if not mesh_is_empty(mesh):
        mesh.save(output_path('meshes', file_stem + '.obj'))
    return mesh


def comparison_row(mesh_name, cloud, frame, file_prefix, first_panels):
    """Panels for one row: first_panels, then the winding number and point-to-plane at each eps."""
    view = MESH_VIEWS[mesh_name]
    winding = WindingNumberImplicit(cloud, WINDING_NEIGHBORS)
    spacing = winding.mean_neighbor_distances.mean()
    panels = list(first_panels)
    panels.append(dict(kind='mesh', view=view, frame=frame, color=METHOD_COLORS['winding number'],
                       data=reconstruct_and_save(winding, f'{file_prefix}_winding'),
                       title='winding number'))
    for ratio in EPSILON_RATIOS_SHOWN:
        suffix = ' (matched)' if ratio == MATCHED_EPSILON_RATIO else ''
        implicit = PointToPlaneImplicit(cloud, ratio * spacing)
        panels.append(dict(kind='mesh', view=view, frame=frame, color=METHOD_COLORS['point-to-plane'],
                           data=reconstruct_and_save(implicit, f'{file_prefix}_plane_eps{ratio:g}h'),
                           title=f'point-to-plane\neps = {ratio:g}h{suffix}'))
    print(f'  {file_prefix}: h = {spacing:.4f}', flush=True)
    return panels


def clean_figure():
    panels = []
    for mesh_name in CLEAN_MESHES:
        mesh = load_normalized_mesh(mesh_name)
        cloud = sample_point_cloud(mesh, NUM_RECON_POINTS, RECON_SEED)
        view = MESH_VIEWS[mesh_name]
        first_panels = [dict(kind='mesh', data=mesh, view=view, title=f'{mesh_name} (ground truth)'),
                        dict(kind='points', data=cloud, view=view, frame=mesh.vertex,
                             title=f'{NUM_RECON_POINTS} samples')]
        panels += comparison_row(mesh_name, cloud, mesh.vertex, f'{mesh_name}_clean', first_panels)
    render_grid(panels, output_path('figures', 'fig3_clean_comparison.png'), num_columns=6)


def noise_figure():
    mesh = load_normalized_mesh('bunny')
    clean_cloud = sample_point_cloud(mesh, NUM_RECON_POINTS, RECON_SEED)
    panels = []
    for noise_sigma in NOISE_LEVELS:
        cloud = add_position_noise(clean_cloud, noise_sigma, NOISE_SEED)
        first_panels = [dict(kind='points', data=cloud, view=MESH_VIEWS['bunny'], frame=mesh.vertex,
                             title=f'noise sigma = {noise_sigma}')]
        panels += comparison_row('bunny', cloud, mesh.vertex, f'bunny_noise{noise_sigma:g}', first_panels)
    render_grid(panels, output_path('figures', 'fig4_noise_comparison.png'), num_columns=5)


def regularization_figure():
    mesh = load_normalized_mesh('bunny')
    clean_cloud = sample_point_cloud(mesh, NUM_RECON_POINTS, RECON_SEED)
    view = MESH_VIEWS['bunny']
    panels = []
    for noise_sigma in [0.0, 0.02]:
        cloud = add_position_noise(clean_cloud, noise_sigma, NOISE_SEED) if noise_sigma > 0 else clean_cloud
        prefix = f'bunny_noise{noise_sigma:g}'
        winding = WindingNumberImplicit(cloud, WINDING_NEIGHBORS)
        spacing = winding.mean_neighbor_distances.mean()
        regularized = WindingNumberImplicit(cloud, WINDING_NEIGHBORS, REGULARIZATION_RATIO * spacing)
        plane = PointToPlaneImplicit(cloud, MATCHED_EPSILON_RATIO * spacing)
        panels += [
            dict(kind='mesh', view=view, frame=mesh.vertex, color=METHOD_COLORS['winding number'],
                 data=reconstruct_mesh(winding)[0], title=f'winding number\nnoise sigma = {noise_sigma}'),
            dict(kind='mesh', view=view, frame=mesh.vertex, color=METHOD_COLORS['winding number'],
                 data=reconstruct_and_save(regularized, f'{prefix}_winding_regularized'),
                 title=f'winding number, regularized\ndelta = {REGULARIZATION_RATIO:g}h'),
            dict(kind='mesh', view=view, frame=mesh.vertex, color=METHOD_COLORS['point-to-plane'],
                 data=reconstruct_mesh(plane)[0], title=f'point-to-plane\neps = {MATCHED_EPSILON_RATIO:g}h'),
        ]
    render_grid(panels, output_path('figures', 'fig5_winding_regularization.png'), num_columns=3)


if __name__ == '__main__':
    print('Figure 3 (clean)', flush=True)
    clean_figure()
    print('Figure 4 (noise)', flush=True)
    noise_figure()
    print('Figure 5 (regularization)', flush=True)
    regularization_figure()
