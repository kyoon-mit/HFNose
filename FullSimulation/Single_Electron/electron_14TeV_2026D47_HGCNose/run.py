#!/usr/bin/env python2.7
# Simply run this python script
# Do python -i run.py to run interactively

from MakeConfigFiles import *

set_env()

# Can also set up parser
nevents = 5000
E_string_list = set_E()

#makeStep1ConfigFiles (E_string_list, nevents)
#makeStep2ConfigFiles (E_string_list, nevents)
#makeStep3ConfigFiles (E_string_list, nevents)
#makeStep4ConfigFiles (E_string_list, nevents)

#runSteps(E_string_list[:4], 1,2,3)
#runSteps(E_string_list[4:8], 1,2,3)
#runSteps(E_string_list[8:11], 1,2,3)
runSteps(E_string_list[11:], 1,2,3)
