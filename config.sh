#!/usr/bin/env bash

# WARNING: DO NOT MOVE THIS FILE. KEEP IT AT THE TOP-LEVEL DIRECTORY.


###### WARNING: DO NOT CHANGE THESE ENVIRONMENT VARIABLES ######

# This is the top directory
export DIRTOP_HFNOSE=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )

# Analysis directory
export DIRANALYSIS_HFNOSE=${DIRTOP_HFNOSE}/Analysis

# Simulation directory
export DIRSIM_HFNOSE=${DIRTOP_HFNOSE}/FullSimulation

################################################################



###### YOU MAY CHANGE THESE ENVIRONMENT VARIABLES ##############

# This is the directory where your simulation data files are stored.
export DIRDATA_HFNOSE=$PWD

################################################################
