#!/bin/bash

top_dir=$1
. ${top_dir}/../config_kyoon.sh
cd ${DIRSIM_HGCNOSE}/Single_Photon/photon_14TeV_2026D60_HGCNose/
eval `scramv1 runtime -sh`

python run.py
