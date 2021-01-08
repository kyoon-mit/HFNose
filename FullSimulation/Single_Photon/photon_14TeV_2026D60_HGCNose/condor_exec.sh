#!/bin/bash

DIRSIM_HGCNOSE=$1
. ${DIRSIM_HGCNOSE}/../config_kyoon.sh
cd ${DIRSIM_HGCNOSE}/Single_Photon/photon_14TeV_2026D60_HGCNose/
eval `scramv1 runtime -sh`

export DIRDATA_HGCNOSE="${DIRDATA_HGCNOSE}_old_setup_no_highpurity"

python run.py
