#!/bin/bash

DIRSIM_HGCNOSE=$1
. ${DIRSIM_HGCNOSE}/../config_kyoon.sh
cd ${DIRSIM_HGCNOSE}/Single_Electron/electron_14TeV_2026D60_HGCNose/
eval `scramv1 runtime -sh`

python run.py
