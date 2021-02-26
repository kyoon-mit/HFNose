#!/bin/bash

cd /afs/cern.ch/user/k/kyoon/CMSSW_11_3_X_2020-12-18-1100/src/HGCNose/FullSimulation/Single_Electron/electron_14TeV_2026D60_HGCNose/
eval `scramv1 runtime -sh`
. /afs/cern.ch/user/k/kyoon/CMSSW_11_3_X_2020-12-18-1100/src/HGCNose/config_kyoon.sh

python run.py
