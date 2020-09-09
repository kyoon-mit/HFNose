#!/usr/bin/env python2.7
# Simply run this python script Do python -i run.py to run interactively

from MakeConfigFiles import *

set_env()

# Can also set up parser
nevents = 10
E_string_list = set_E('50', '500')

makeStep1ConfigFiles (E_string_list, nevents)
makeStep2ConfigFiles (E_string_list, nevents)
makeStep3ConfigFiles (E_string_list, nevents)
makeStep4ConfigFiles (E_string_list, nevents)

runSteps(E_string_list, 2,3)
