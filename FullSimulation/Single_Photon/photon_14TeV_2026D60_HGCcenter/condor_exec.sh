#!/bin/bash

this_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. this_dir/../../../config_kyoon.sh
cd ${DIRSIM_HGCNOSE}/Single_Photon/photon_14TeV_2026D60_HGCcenter/
eval `scramv1 runtime -sh`

python run.py # Modify the run file directly for setting parameters
