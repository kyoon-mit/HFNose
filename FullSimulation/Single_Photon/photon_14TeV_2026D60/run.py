#!/usr/bin/env python3

from MakeConfigFiles import *
import os
import argparse
from sys import version_info

if version_info[0] < 3:
    raise Exception("Must be using Python 3")

parser = argparse.ArgumentParser(description='Create and process workflow with user-given parameters.')
parser.add_argument('-e', '--energy', type=int, nargs='+', default=[50])
parser.add_argument('-y', '--eta', type=float, nargs='+', default=[3.5])
parser.add_argument('-n', '--nevents', type=int, nargs='?', default=500)
parser.add_argument('-l', '--lumi', type=int, nargs='?', default=0)
parser.add_argument('-t', '--agetag', type=str, nargs='?', default='')
parser.add_argument('-d', '--savedir', type=str, nargs='?', default=os.environ.get('DIRDATA_HGCNOSE') + '/photon_2026D60_HGCNose')
parser.add_argument('-s', '--steps', type=int, nargs='+', default=[1,2,3])
args = parser.parse_args()

print ("Running with the following parameters:", end=" ")
for arg in vars(args):
    print (arg, getattr(args, arg), end=" ")
print ('\n')
type (args.nevents)

run = RunInstance()
run.set_E(*args.energy)
run.set_eta(*args.eta)
run.set_nevents(args.nevents)
run.set_aging(args.lumi, args.agetag)
run.set_top_save_dir(args.savedir)
run.makeAllConfigFiles()
run.runSteps(*args.steps)
