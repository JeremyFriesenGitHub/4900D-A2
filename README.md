# COMP 4900D - A2

## Setup

```
pip install -r requirements.txt                              
git clone https://github.com/ovankaic/GeomProc ~/GeomProc
ln -s ~/GeomProc/geomproc geomproc
```

## Running

```
./run_all.sh        
```


## AI Usage

I used an AI while working on this assignment. I have read and tested all of this code. The breakdown is below:

**Written by the AI**:
  - reconstruction.py, which implements the point-to-plane and winding-number
    implicit functions, the k-nearest-neighbour area weights, sampling and noise,
    the marching-cubes wrapper and the distance metrics.
  - exp1_weight_locality.py, exp2_epsilon_sweep.py, exp3_visual_comparison.py and
    exp4_quantitative.py, which run the experiments and write the figures, tables
    and meshes under results/.
  - render.py, which draws the figures, and run_all.sh, which runs all four
    experiments.

**Written solely by me**:
- I ran experiments that produced results/, and found some
    implementation issues that shaped the method: the handout weight
    1/(r² + ε²) produced no surface, the area weights need k = 6 neighbours,
    GeomProc's marching cubes mis-indexed the grid for some sizes, and
    point-to-plane leaves spurious surfaces far from the samples.

**Report text written with AI help**:
  - TODO: which report sections, if any, the AI helped with.
