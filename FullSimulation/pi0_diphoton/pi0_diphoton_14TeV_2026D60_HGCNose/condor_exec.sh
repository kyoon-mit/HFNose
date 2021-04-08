#!/bin/bash

DIRSIM_HGCNOSE=$1
. ${DIRSIM_HGCNOSE}/../config_kyoon.sh
cd ${DIRSIM_HGCNOSE}/pi0_diphoton/pi0_diphoton_14TeV_2026D60_HGCNose/
eval `scramv1 runtime -sh`

export DIRDATA_HGCNOSE=${DIRDATA_HGCNOSE}

python run.py
