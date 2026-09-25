"""Minimal offscreen renderer for report figures (matplotlib only).

Orthographic projection, back-face culling, painter's algorithm and flat
Lambert shading, so faceting and marching-cubes artifacts stay visible.
MeshLab screenshots of the saved .obj files work just as well.
"""

import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

# Per-mesh viewing setup: which model axis points up, and the camera angles.
MESH_VIEWS = {
    'bunny': dict(up_axis='y', azimuth=20, elevation=10),
    'camel': dict(up_axis='y', azimuth=-50, elevation=15),
    'kid': dict(up_axis='z', azimuth=240, elevation=10),
    'sphere': dict(up_axis='y', azimuth=30, elevation=20),
    'cube': dict(up_axis='y', azimuth=35, elevation=25),
}
DEFAULT_VIEW = dict(up_axis='y', azimuth=30, elevation=20)

LIGHT_DIRECTION = np.array([-0.4, 0.6, 1.0]) / np.linalg.norm([-0.4, 0.6, 1.0])
METHOD_COLORS = {
    'ground truth': np.array([0.80, 0.80, 0.80]),
    'winding number': np.array([0.95, 0.62, 0.35]),
    'point-to-plane': np.array([0.45, 0.65, 0.90]),
}


def view_rotation(up_axis='y', azimuth=30, elevation=20):
    """Rotation taking model coordinates to camera coordinates
    (x right, y up, z towards the viewer)."""
    up_to_y = {
        'y': np.eye(3),
        'z': np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]], dtype=float),
        'x': np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]], dtype=float),
    }[up_axis]
    azimuth_rad, elevation_rad = np.radians(azimuth), np.radians(elevation)
    spin = np.array([[np.cos(azimuth_rad), 0, np.sin(azimuth_rad)],
                     [0, 1, 0],
                     [-np.sin(azimuth_rad), 0, np.cos(azimuth_rad)]])
    tilt = np.array([[1, 0, 0],
                     [0, np.cos(elevation_rad), -np.sin(elevation_rad)],
                     [0, np.sin(elevation_rad), np.cos(elevation_rad)]])
    return tilt @ spin @ up_to_y


def draw_mesh(axes, mesh, rotation, base_color, title=None, frame=None):
    if mesh is None or mesh.face.shape[0] == 0:
        axes.text(0.5, 0.5, '(empty surface)', ha='center', va='center', transform=axes.transAxes)
        axes.set_axis_off()
        if title:
            axes.set_title(title, fontsize=9)
        return
    camera_vertices = mesh.vertex @ rotation.T
    triangles = camera_vertices[mesh.face]
    face_normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    normal_lengths = np.linalg.norm(face_normals, axis=1)
    valid = normal_lengths > 1e-14
    face_normals[valid] /= normal_lengths[valid, None]
    facing_viewer = valid & (face_normals[:, 2] > 0)
    triangles, face_normals = triangles[facing_viewer], face_normals[facing_viewer]

    # Painter's algorithm: far triangles first
    draw_order = np.argsort(triangles[:, :, 2].mean(axis=1))
    triangles, face_normals = triangles[draw_order], face_normals[draw_order]
    diffuse = np.clip(face_normals @ LIGHT_DIRECTION, 0.0, 1.0)
    shade = 0.25 + 0.75 * diffuse
    face_colors = np.clip(shade[:, None] * base_color[None, :], 0.0, 1.0)

    collection = PolyCollection(triangles[:, :, :2], facecolors=face_colors, edgecolors=face_colors,
                                linewidths=0.1, antialiaseds=False)
    axes.add_collection(collection)
    finish_axes(axes, camera_vertices if frame is None else frame @ rotation.T, title)


def draw_point_cloud(axes, point_cloud, rotation, title=None, frame=None, marker_size=1.5):
    camera_points = point_cloud.point @ rotation.T
    camera_normals = point_cloud.normal @ rotation.T
    draw_order = np.argsort(camera_points[:, 2])
    shade = 0.3 + 0.7 * np.abs(camera_normals[draw_order] @ LIGHT_DIRECTION)
    axes.scatter(camera_points[draw_order, 0], camera_points[draw_order, 1], s=marker_size,
                 c=shade, cmap='gray', vmin=0, vmax=1.15, linewidths=0)
    finish_axes(axes, camera_points if frame is None else frame @ rotation.T, title)


def finish_axes(axes, camera_points, title):
    center = 0.5 * (camera_points[:, :2].min(axis=0) + camera_points[:, :2].max(axis=0))
    half_extent = 0.58 * (camera_points[:, :2].max(axis=0) - camera_points[:, :2].min(axis=0)).max()
    axes.set_xlim(center[0] - half_extent, center[0] + half_extent)  # content outside is clipped
    axes.set_ylim(center[1] - half_extent, center[1] + half_extent)
    axes.set_aspect('equal')
    axes.set_axis_off()
    if title:
        axes.set_title(title, fontsize=9)


def render_grid(panels, output_path, num_columns, panel_size=2.6, row_labels=None):
    """panels: list of dicts with keys kind ('mesh' | 'points'), data, view, title, color, and
    optionally frame: points whose extent sets the visible window (e.g. the ground-truth
    vertices, so every panel of a row has the same scale and far-away artifacts get clipped)."""
    num_rows = int(np.ceil(len(panels) / num_columns))
    figure, axes_grid = plt.subplots(num_rows, num_columns,
                                     figsize=(panel_size * num_columns, panel_size * num_rows + 0.2),
                                     squeeze=False)
    for panel_index, axes in enumerate(axes_grid.ravel()):
        if panel_index >= len(panels):
            axes.set_axis_off()
            continue
        panel = panels[panel_index]
        rotation = view_rotation(**panel.get('view', DEFAULT_VIEW))
        if panel['kind'] == 'points':
            draw_point_cloud(axes, panel['data'], rotation, panel.get('title'), panel.get('frame'))
        else:
            draw_mesh(axes, panel['data'], rotation, panel.get('color', METHOD_COLORS['ground truth']),
                      panel.get('title'), panel.get('frame'))
    if row_labels:
        for row_index, label in enumerate(row_labels):
            axes_grid[row_index, 0].text(-0.08, 0.5, label, rotation=90, ha='center', va='center',
                                         transform=axes_grid[row_index, 0].transAxes, fontsize=10)
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
