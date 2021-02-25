#!/usr/bin/env python2.7
# Simply run this python script
# Do python -i run.py to run interactively

from MakeConfigFiles import *
import os

run = RunInstance()
run.set_E(200)
run.set_eta(3.5)
run.set_nevents(500)
run.set_top_save_dir(os.environ.get('DIRDATA_HGCNOSE') + '/photon_2026D60_HGCNose_ageing_3000')
run.makeAllConfigFiles()
run.runSteps(1,2,3)
