#!/usr/bin/env python2

# Configuration file for the EnergyResolution.h analysis

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process('HFNose')

options = VarParsing ('analysis')

pt=1

OUTPUT_DIR = '/home/kyoon/CMSSW_11_1_0_pre7_RECHIT/src/HGCNose/EnergyResolution/output/'
INPUT_DIR = 'file:/data/t3home000/kyoon/gendata/photon_2026D47/'
options.outputFile = OUTPUT_DIR + 'ERes_pt{}.root'.format(pt)
options.inputFiles = INPUT_DIR + 'photon_pt{0}/step3_photon_pt{0}.root'.format(pt)
options.parseArguments()

process.load('FWCore.MessageService.MessageLogger_cfi')
# Important to load geometry configs because they are part of EventSetup
# Not sure why the following works
process.load('Geometry.ForwardCommonData.hfnoseXML_cfi')
process.load('Geometry.ForwardCommonData.hfnoseParametersInitialization_cfi')
process.load('Geometry.ForwardCommonData.hfnoseNumberingInitialization_cfi')
process.load('Geometry.CaloEventSetup.HFNoseTopology_cfi')
process.load('Geometry.ForwardGeometry.HFNoseGeometryESProducer_cfi')

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.analysis = cms.EDAnalyzer('EnergyResolution',
    TAG_HGCHFNoseRecHits = cms.untracked.InputTag('HGCalRecHit', 'HGCHFNoseRecHits')
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string(options.outputFile)
)

process.p = cms.Path(process.analysis)
