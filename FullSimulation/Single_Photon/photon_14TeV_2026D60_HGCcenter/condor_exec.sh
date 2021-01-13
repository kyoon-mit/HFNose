#!/bin/bash

DIRSIM_HGCNOSE=$1
. ${DIRSIM_HGCNOSE}/../config_kyoon.sh
cd ${DIRSIM_HGCNOSE}/Single_Photon/photon_14TeV_2026D60_HGCcenter/
eval `scramv1 runtime -sh`

python run.py # Modify the run file directly for setting parameters
