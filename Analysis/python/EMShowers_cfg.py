#!/usr/bin/env python2

# Configuration file for the EnergyResolution.h analysis

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.Eras.Era_Phase2C10_cff import Phase2C10

process = cms.Process('AnalysisHGCNoseMoliereRadius', Phase2C10)

options = VarParsing ('analysis')
options.register ('pt', '1', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                    "(type: string) pt value of the photon")
options.parseArguments()

OUTPUT_DIR = '/home/kyoon/CMSSW_11_1_0_pre7_RECHIT/src/HGCNose/Analysis/output/'
INPUT_DIR = 'file:/data/t3home000/kyoon/gendata/photon_2026D47/'
# options.outputFile = OUTPUT_DIR + 'ERes_pt{}.root'.format(pt)
# options.inputFiles = INPUT_DIR + 'photon_pt{0}/step3_photon_pt{0}.root'.format(pt)
outputfile = OUTPUT_DIR + 'EMShowers_pt{}.root'.format(options.pt)
inputfile = INPUT_DIR + 'photon_pt{0}/step3_photon_pt{0}.root'.format(options.pt)

# Process load
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load('Configuration.Geometry.GeometryExtended2026D47_cff')
process.load('Configuration.Geometry.GeometryExtended2026D47Reco_cff')

process.load('Geometry.ForwardCommonData.hfnoseXML_cfi')
process.load('Geometry.ForwardCommonData.hfnoseParametersInitialization_cfi')
process.load('Geometry.ForwardCommonData.hfnoseNumberingInitialization_cfi')
process.load('Geometry.CaloEventSetup.HFNoseTopology_cfi')
process.load('Geometry.ForwardGeometry.HFNoseGeometryESProducer_cfi')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options = cms.untracked.PSet (
    wantSummary = cms.untracked.bool(False),
    numberOfThreads = cms.untracked.uint32(12),
    numberOfStreams = cms.untracked.uint32(12)
)

# Global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

# Input
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(inputfile)
)

# Max events to process
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Process
process.Analysis_EMShowers = cms.EDAnalyzer('EMShowerStudies',
    #TAG_HGCHFNoseRecHits = cms.untracked.InputTag('HGCalRecHit', 'HGCHFNoseRecHits')
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string(outputfile)
)

# Timing
process.Timing = cms.Service("Timing",
  summaryOnly = cms.untracked.bool(True),
  useJobReport = cms.untracked.bool(True)
)

process.p = cms.Path(process.Analysis_EMShowers)
