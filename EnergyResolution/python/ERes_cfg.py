#!/usr/bin/env python2

# Configuration file for the EnergyResolution.h analysis

import FWCore.ParameterSet.Config as cms

process = cms.Process("HFNose")

process.load("FWCore.MessageService.MessageLogger_cfi")
# Important to load geometry configs because they are part of EventSetup
# Not sure why the following works
process.load("Geometry.ForwardCommonData.hfnoseXML_cfi")
process.load("Geometry.ForwardCommonData.hfnoseParametersInitialization_cfi")
process.load("Geometry.ForwardCommonData.hfnoseNumberingInitialization_cfi")
process.load("Geometry.CaloEventSetup.HFNoseTopology_cfi")
process.load("Geometry.ForwardGeometry.HFNoseGeometryESProducer_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:$DIR_DATA/photon_pt2/step3_photon_pt2.root")
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(30)
)

process.myProcess = cms.EDAnalyzer("EnergyResolution",
    TAG_HGCHFNoseRecHits = cms.untracked.InputTag("HGCalRecHit", "HGCHFNoseRecHits")
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("file:/home/kyoon/CMSSW_11_1_0_pre7_RECHIT/src/HFNose/EnergyResolution/output/ERes_test3.root")
)

process.p = cms.Path(process.myProcess)
