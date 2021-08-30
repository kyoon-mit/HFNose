# HFNose

Author: K.Yoon

## Setup
1. Put it in your working directory and build using `scram b`.
2. Set environment variables: edit `config.sh` to suit your preferences and execute.

## FullSimulation

## Analysis
Here is the list of analysis plugins.
* `Analysis/plugins/TICLAnalyzer.cc` | Generic module for analyzing TICL-related objects, such as Tracksters. Run using `cmsRun python/TICLAnalyzer.py`. (Add: explanation of variables)
* `Analysis/plugins/EnergyResolution.cc` | Returns energy resolution for a given PID and eta range. Outputs plots in `output` directory. Run using `cmsRun python/ERes_cfg.py pt=<pt value of photon>`.
* `Analysis/plugins/PropagatorAnalyzer.cc` | Analysis for propagator

## Plotting
List of plotting modules.
* `Analysis/scripts/TICLAnalyzer.py` | Plots TICL related histograms
* `Analysis/scripts/ERes_FitPlot.py` | Plots energy resolution, mean energy, and detector layer-related plots. PNG images are saved in `plots` directory.

## NanoAODProduction
Flattens EDM-formatted files to NanoAOD Ntuples.

## MLProjects
Workflow:
1. `MLProjects/python/HFNoseNtuplizer_cfg.py` | Creates Ntuples from EDM format.
2. `MLProjects/scripts/Train.py` | Preprocesses the above Ntuples and trains model.
 
