# COMP 4900D · Assignment 2: Surface reconstruction with the winding number

This repository reconstructs surfaces from oriented point clouds with two implicit functions and compares them visually and quantitatively. The first is the point-to-plane blend from class, Φ(x) = Σᵢ φ(‖x − pᵢ‖) nᵢᵀ(x − pᵢ) / Σᵢ φ(‖x − pᵢ‖). The second is the winding number of Barill et al. (SIGGRAPH 2018), w(x) = 1/(4π) Σᵢ Aᵢ nᵢᵀ(pᵢ − x) / ‖pᵢ − x‖³, used as Φ(x) = 0.5 − w(x) with Aᵢ = (mean distance from pᵢ to its k nearest neighbours)². Surfaces are extracted with GeomProc's marching cubes on 64³ cubes over [−1.5, 1.5]³, with every mesh normalized to [−1, 1].

## Setup

```
pip install -r requirements.txt                              # numpy, scipy, matplotlib
git clone https://github.com/ovankaic/GeomProc ~/GeomProc
ln -s ~/GeomProc/geomproc geomproc
```

`geomproc` in this folder is a symlink to the library's package, as in assignment 1, and git ignores it. GeomProc is used unmodified, and every input mesh (sphere, cube, bunny, camel, kid) comes from the `meshes/` folder of the same clone, so neither is included here. Use a copy that includes the marching-cubes fix from September 24, 2026 (`git pull` an older clone).

## Running

```
./run_all.sh        # about 40 minutes on one core
```

| Script | Output (under `results/`) |
|---|---|
| `exp1_weight_locality.py` | `figures/fig1_weight_locality.png`, `tables/area_weight_check.csv` |
| `exp2_epsilon_sweep.py` | `figures/fig2_epsilon_sweep.png`, `tables/epsilon_sweep.csv` |
| `exp3_visual_comparison.py` | `figures/fig3_clean_comparison.png`, `fig4_noise_comparison.png`, `fig5_winding_regularization.png`, `fig6_weight_comparison.png`, `meshes/*.obj` |
| `exp4_quantitative.py` | `tables/quantitative.md`, `tables/quantitative.csv` |

`reconstruction.py` holds both methods and the shared helpers, and `render.py` draws the figures (MeshLab screenshots of `results/meshes/*.obj` work too). Every experiment uses 2000 reconstruction samples, as in GeomProc's RBF demo, and fixed seeds. `python3 exp2_epsilon_sweep.py --plot-only` redraws fig. 2 from the saved table without the 15-minute sweep.

## Findings

Summed over all samples, the handout's weight φ(r) = 1/(r² + ε²) produces no surface (fig. 1), so the point-to-plane blend uses a faster-decaying weight, φ(r) = 1/(r⁶ + ε⁶) with ε = 0.5h (fig. 6 compares it with the squared weight 1/(r² + ε²)²). Near the surface it is the most accurate method: on clean samples the mean distance from the test points to the reconstruction is 0.0015–0.0033, against 0.0040–0.0077 for the winding number (`tables/quantitative.md`), and with noise of σ = 0.02 the two are about equal (0.0083–0.0113 against 0.0082–0.0112). The difference is away from the surface. Point-to-plane only gets the sign right close to the samples, so on the kid, whose feet are sparsely sampled, it still leaves a spurious sheet (mean distance from reconstructed vertices to the true surface 0.28, against 0.006 for the winding number), and on the camel a few fragments (0.023). Slower weights are much worse there: 0.36 and 0.62 for the squared weight, 0.35 and 0.67 for the exponential weight exp(−r²/ε²). The winding number's weakness is a rough surface: its kernel is singular at every sample, which also makes its maximum |Φ| at the test points large (4.7–139 on clean samples) while the median stays at 0.012–0.043. Regularizing the kernel removes the roughness (fig. 5). In the ε sweep, the point-to-plane surface stays about equally close to the winding-number surface up to 0.5h and drifts away beyond it (fig. 2).

## Implementation choices that affect the results

**The point-to-plane weight falls off fast.** With the handout's φ(r) = 1/(r² + ε²) summed over all samples, Φ is negative on both sides of the surface, so marching cubes returns an empty mesh for every ε (`exp1`, fig. 1). Summed over a surface, 1/r² behaves like ∫ dA/r², which keeps growing with the area included, so distant samples dominate the blend; on a closed surface their plane distances average to about −3V/A < 0. Faster-decaying weights make nearby samples dominate. The main method, `SixthPowerPointToPlaneImplicit`, uses φ(r) = 1/(r⁶ + ε⁶). `PointToPlaneImplicit(cloud, epsilon, weight_power, num_neighbors)` covers the handout weight, the squared weight 1/(r² + ε²)² and the 16-nearest-neighbour version, and `GaussianPointToPlaneImplicit` the exponential weight; `exp4` reports all of them for reference. The exponential weight is computed relative to the nearest sample, which leaves Φ unchanged but avoids 0/0 where every exp(−r²/ε²) underflows to zero; computed directly, 28% of the bunny's grid corners come out as NaN at ε = 0.5h.

**k = 6 neighbours for Aᵢ.** For uniformly random samples, Σ Aᵢ matches the surface area at k = 6 (0.95–1.02 on all five meshes, `tables/area_weight_check.csv`). With k = 1 it is about 0.32 of the area, so w stays near 0.32 inside and never crosses 0.5. The neighbours come from GeomProc's KD-tree (`geomproc.KDTree`), and the query skips the point itself. scipy's `cKDTree` is used only where a pure-Python tree would be too slow: the 16-nearest-sample variant, which queries every grid corner, and the distances to meshes in `exp2` and `exp4`.

**ε is a multiple of the sample spacing h**, the mean kNN distance behind Aᵢ, so one setting transfers across meshes. `exp2` sweeps ε/h for the main weight, and `MATCHED_EPSILON_RATIO` records the result (0.5h); the squared-weight reference uses 0.1h, from its own earlier sweep.

**Marching cubes.** `reconstruct_mesh` evaluates Φ on the whole grid in one vectorized pass and gives `geomproc.marching_cubes` a lookup function. The output is identical to passing Φ directly (checked on a 16³ grid), only much faster. GeomProc versions before the September 24, 2026 fix controlled the grid loops with floating-point comparisons and mis-indexed the grid for some sizes. The default 3/64 cell is exact in binary and was never affected (the 64³ bunny reconstruction is identical before and after the fix), and `reconstruct_mesh` raises an error if GeomProc does not evaluate every grid corner exactly once.

**Sampling and noise.** `geomproc.mesh.sample` draws from Python's `random` module, which is what gets seeded. Noise is zero-mean Gaussian on the positions only, with normals kept from the clean mesh, to isolate positional noise; `geomproc.mesh.add_noise` was avoided because it adds uniform noise in [0, scale), which also shifts the shape.

**Error metrics.** As required, the error is |Φ(x)| at 10× more test points sampled independently from the surface, reported as max and mean, with the median added because the winding number's maximum comes from test points that land next to a sample. Since |0.5 − w| is not a length, `exp4` also measures on the reconstructed meshes the distance from test points to the reconstruction and from reconstructed vertices to the true mesh; the second one exposes spurious surfaces.

**Extension, off by default.** `WindingNumberImplicit(cloud, k, kernel_regularization=δ)` replaces ‖p − x‖³ with (‖p − x‖² + δ²)^{3/2}, which removes the spikes at the samples (fig. 5). δ plays the same role for the winding number that ε plays for point-to-plane.

## AI use

As the course policy requires, this records the AI help on this assignment. All code in this repository (`reconstruction.py`, `render.py`, `exp1`–`exp4`, `run_all.sh`) was drafted by Claude (Anthropic's assistant, claude.ai) in September 2026. Claude also ran the experiments that produced `results/` and identified the implementation issues described above: the handout weight, the choice of k, GeomProc's grid constraint and the far-field failures of point-to-plane.

*[Before submitting: add what you reviewed, changed or re-ran yourself, and whether AI helped with the report.]*
