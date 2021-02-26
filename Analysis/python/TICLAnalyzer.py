#!/usr/bin/env python2

# Configuration file for the EnergyResolution.h analysis

import os

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.Eras.Era_Phase2C10_cff import Phase2C10

process = cms.Process('TICLAnalyzer', Phase2C10)

# Varparsing
options = VarParsing('analysis')
options.register('inpath', 'photon_2026D60_HGCNose', VarParsing.multiplicity.singleton,
                    VarParsing.varType.string,
                    "last part of input directory path")
options.register('outpath', 'photon_2026D60_HGCNose', VarParsing.multiplicity.singleton,
                    VarParsing.varType.string,
                    "last part of input directory path")
options.register('outsuffix', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                    "outfile suffix to put at the very end before the extension")
options.register('pid', 22, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                    "(type: int) PDG ID of the particle")
options.register('E', 100, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                    "(type: int) E value of the particle")
options.register('eta', 3.5, VarParsing.multiplicity.singleton, VarParsing.varType.float,
                    "(type: float) eta value of the particle")
options.register('deltaR', 0.5, VarParsing.multiplicity.singleton, VarParsing.varType.float,
                    "(type: float) R value of the shower cone")
options.register('itername', 'EMn', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                    "(type: string) itername of the Trackster")
options.parseArguments()

# Set output and input paths
OUTPUT_DIR = os.path.abspath(os.environ['DIRANALYSIS_HGCNOSE'] + '/output/{}/'.format(options.outpath))
INPUT_DIR = os.path.abspath(os.environ['DIRDATA_HGCNOSE'] + '/{}/'.format(options.inpath))
outputfile = OUTPUT_DIR + '/Single_Photon_E{0}_eta{1}{2}.root'.format(options.E, options.eta, options.outsuffix)
inputfile = 'file:' + INPUT_DIR + '/photon_E{0}_eta{1}/step3_photon_E{0}_eta{1}.root'.format(options.E, options.eta)

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Process load
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load('Configuration.Geometry.GeometryExtended2026D60_cff')
process.load('Configuration.Geometry.GeometryExtended2026D60Reco_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options = cms.untracked.PSet (
    wantSummary = cms.untracked.bool(False),
    numberOfThreads = cms.untracked.uint32(4),
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
process.Analysis_TICLAnalyzer = cms.EDAnalyzer('TICLAnalyzer',
    TAG_Trackster = cms.untracked.InputTag('ticlTrackstersHFNoseEM'),
    select_PID = cms.untracked.int32(options.pid),
    select_EtaLow = cms.untracked.double(options.eta - 0.1),
    select_EtaHigh = cms.untracked.double(options.eta + 0.1),
    truth_matching_deltaR = cms.untracked.double(options.deltaR),
    trackster_itername = cms.untracked.string(options.itername)
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string(outputfile)
)

# Timing
process.Timing = cms.Service("Timing",
  summaryOnly = cms.untracked.bool(True),
  useJobReport = cms.untracked.bool(True)
)

process.p = cms.Path(process.Analysis_TICLAnalyzer)
