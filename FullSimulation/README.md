# HGCNose

Author: K.Yoon

## Setup
1. Put it in your working directory and build using `scram b`.
2. Set environment variables: edit `config.sh` to suit your preferences and execute.

## Analysis
Here is the list of analysis modules.
* `plugins/TICLAnalyzer.cc` | Generic module for analyzing TICL-related objects, such as Tracksters. Run using `cmsRun python/TICLAnalyzer.py`. (Add: explanation of variables)
* `plugins/EnergyResolution.cc` | Returns energy resolution for a given PID and eta range. Outputs plots in `output` directory. Run using `cmsRun python/ERes_cfg.py pt=<pt value of photon>`.
* `plugins/MoliereRadius.cc` | (working) Outputs plots related to Moliere radius
*

## Plotting
List of plotting modules.
* `scripts/TICLAnalyzer.py` | Plots TICL related histograms
* `scripts/ERes_FitPlot.py` | Plots energy resolution, mean energy, and detector layer-related plots. PNG images are saved in `plots` directory.

## NanoAODProduction
Flattens EDM-formatted files to NanoAOD Ntuples. Consider this as the data-preprocessing part for ML projects.
