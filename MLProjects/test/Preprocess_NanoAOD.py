#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
This module converts a ROOT file which contains the NanoAOD format of the data which was flattened and converted from a CMSSW EDM file.

Usage:
    If running as script, run using the following example.
    If importing as module, refer to Functions.

Command-line example:
    $ python Preproces_NanoAOD.py --infile HFNoseNanoAOD.root --outfile 
    
Functions:
    

Documentations:
    NanoAOD  |  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD
    uproot   |  https://pypi.org/project/uproot/
"""

import uproot4 as uproot
import pkg_resources

pkg_resources.require('uproot4 > 0.1.')





if __name__ == '__main__':
    
