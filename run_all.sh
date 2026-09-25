#!/usr/bin/env bash
# Reproduces every figure and table in results/ (about 30 minutes on one core).
# Needs the geomproc symlink in this folder (see README).
set -e
cd "$(dirname "$0")"
python3 exp1_weight_locality.py
python3 exp2_epsilon_sweep.py
python3 exp3_visual_comparison.py
python3 exp4_quantitative.py
