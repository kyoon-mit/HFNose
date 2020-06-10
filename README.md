# HGCNose

Author: K.Yoon

## Setup
1. Put it in your working directory and build using `scram b`.
2. Set environment variables (will put script later). For now, change DIR_DATA and TOP_DIR in `python/ERes_cfg.py` and `scripts/ERes_FitPlot.py`, respectively.

## Analysis
Here is the list of analysis modules.
* plugins/EnergyResolution.cc | Returns energy resolution for a given PID and eta range. Outputs plots in `output` directory. Run using `python/ERes_cfg.py`.

## Plotting
List of plotting modules.
* scripts/ERes_FitPlot.py | Plots energy resolution, mean energy, and detector layer-related plots. PNG images are saved in `plots` directory.

