#!/usr/bin/env python2.7
# Simply run this python script
# Do python -i run.py to run interactively

from MakeConfigFiles import *

set_env()

# Can also set up parser
nevents = 10
pt_string_list = set_pt('4','5','6','7')

makeStep1ConfigFiles (pt_string_list, nevents)
makeStep2ConfigFiles (pt_string_list, nevents)
makeStep3ConfigFiles (pt_string_list, nevents)
makeStep4ConfigFiles (pt_string_list, nevents)

runSteps(pt_string_list, 1,2,3)
