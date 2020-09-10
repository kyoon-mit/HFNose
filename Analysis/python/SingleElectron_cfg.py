#!/usr/bin/env python2

# Configuration file for the EnergyResolution.h analysis

import os

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.Eras.Era_Phase2C10_cff import Phase2C10

process = cms.Process('AnalysisHGCNoseSingleElectron', Phase2C10)

# Varparsing
options = VarParsing('analysis')
options.register('E', '50', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                    "(type: string) E value of the electron")
options.parseArguments()

# Set output and input paths
OUTPUT_DIR = os.path.abspath(os.environ['DIRANALYSIS_HGCNOSE'] + '/output')
INPUT_DIR = os.path.abspath(os.environ['DIRDATA_HGCNOSE'] + '/electron_2026D60/electron_E{}'.format(options.E))
outputfile = OUTPUT_DIR + '/Single_Electron_E{}.root'.format(options.E)
inputfile = 'file:' + INPUT_DIR + '/step3_electron_E{}.root'.format(options.E)

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Process load
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load('Configuration.Geometry.GeometryExtended2026D60_cff')
process.load('Configuration.Geometry.GeometryExtended2026D60Reco_cff')

process.load('Geometry.ForwardCommonData.hfnoseXML_cfi')
process.load('Geometry.ForwardCommonData.hfnoseParametersInitialization_cfi')
process.load('Geometry.ForwardCommonData.hfnoseNumberingInitialization_cfi')
process.load('Geometry.CaloEventSetup.HFNoseTopology_cfi')
process.load('Geometry.ForwardGeometry.HFNoseGeometryESProducer_cfi')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options = cms.untracked.PSet (
    wantSummary = cms.untracked.bool(False),
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0)
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
process.Analysis_SingleElectron = cms.EDAnalyzer('SingleElectron',
    # TAG_HGCHFNoseRecHits = cms.untracked.InputTag('HGCalRecHit', 'HGCHFNoseRecHits')
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string(outputfile)
)

# Timing
process.Timing = cms.Service("Timing",
  summaryOnly = cms.untracked.bool(True),
  useJobReport = cms.untracked.bool(True)
)

process.p = cms.Path(process.Analysis_SingleElectron)
