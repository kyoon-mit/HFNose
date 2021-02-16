#!/usr/bin/env python2.7
# Simply run this python script
# Do python -i run.py to run interactively

from MakeConfigFiles import *
import os

run = RunInstance()
run.set_E(50, 100, 200, 300)
run.set_eta(2.95)
run.set_nevents(500)
run.set_top_save_dir(os.environ.get('DIRDATA_HGCNOSE'))
run.makeAllConfigFiles()
run.runSteps(1,2,3)
