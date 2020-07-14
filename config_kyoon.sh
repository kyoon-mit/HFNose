#!/usr/bin/env bash

# WARNING: DO NOT MOVE THIS FILE. KEEP IT AT THE TOP-LEVEL DIRECTORY.


###### WARNING: DO NOT CHANGE THESE ENVIRONMENT VARIABLES ######

# This is the top directory
export DIRTOP_HGCNOSE=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )

# Analysis directory
export DIRANALYSIS_HGCNOSE=${DIRTOP_HGCNOSE}/Analysis

# Simulation directory
export DIRSIM_HGCNOSE=${DIRTOP_HGCNOSE}/FullSimulation

################################################################



###### YOU MAY CHANGE THESE ENVIRONMENT VARIABLES ##############

# This is the directory where your simulation data files are stored.
export DIRDATA_HGCNOSE=/data/t3home000/kyoon/TICL

################################################################
