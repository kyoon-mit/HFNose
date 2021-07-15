#!/bin/bash

top_dir=$1
. ${top_dir}/config_kyoon.sh
cd ${DIRSIM_HFNOSE}/Single_Electron/electron_14TeV_2026D60/
eval `scramv1 runtime -sh`

./run.py -e $2 -y $3 -n $4 -l $5 -d $6 -s 1 2 3
