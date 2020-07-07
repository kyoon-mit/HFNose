# HGCNose/FullSimulation

Author: K.Yoon

## Notes
The default setting in CMSSW_11_1_0_pre7 will not generate calotruth in HGCNose. In order to properly do so, follow these steps:
  * `cd $CMSSW_BASE`
  * `git cms-addpkg SimGeneral/MixingModule`
  * `cd src`
  * `mkdir SimGeneral`
  * `mv MixingModule SimGeneral`
  * Uncomment lines 44-51 in SimGeneral/MixingModule/python/caloTruthProducer_cfi.py
  * `cd $CMSSW_BASE/src/SimGeneral`
  * `scram b`
  
You may check that the module has been successfully added by running python and importing the module.
  * `python`
  * `import SimGeneral.MixingModule.python.caloTruthProducer_cfi`
  * `print (SimGeneral.MixingModule.python.caloTruthProducer_cfi.__file__)`
