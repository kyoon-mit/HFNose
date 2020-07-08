#!/usr/bin/env python2.7
# Simply run this python script
# For interactive session, do: python -i run.py


from MakeConfigFiles import *

set_env()

# Can also set up parser
nevents = 5000
E_string_list = set_E('80' ,'120', '160', '200')

makeStep1ConfigFiles (E_string_list, nevents)
makeStep2ConfigFiles (E_string_list, nevents)
makeStep3ConfigFiles (E_string_list, nevents)
makeStep4ConfigFiles (E_string_list, nevents)

runSteps(E_string_list, 1,2,3)
