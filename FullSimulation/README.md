# HGCNose/FullSimulation

Author: K.Yoon

## Notes
* The default setting in CMSSW_11_1_0_pre7 will not generate calotruth in HGCNose. In order to properly do so, follow these steps:
  * cd $CMSSW_BASE/src
  * mkdir SimGeneral
  * cd SimGeneral
  * git cms-addpkg SimGeneral/MixingModule
  * Uncomment lines 44-51 in MixingModule/python/caloTruthProducer_cfi.python
  * cd $CMSSW_BASE/src/SimGeneral
  * scram b
