import os
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C10_cff import Phase2C10

process = cms.Process('HFNoseNtuplizer', Phase2C10)

#-----------------------------------------
# VarParsing
#-----------------------------------------

options = VarParsing('analysis')
options.register('process',
                            'photon',
                            VarParsing.multiplicity.singleton,
                            VarParsing.varType.string,
        		            "(type: string) process name")
options.register('E',
                            100,
                            VarParsing.multiplicity.singleton,
                            VarParsing.varType.int,
                            "(type: int) E value of the particle")
options.register('eta',
                            3.5,
                            VarParsing.multiplicity.singleton,
                            VarParsing.varType.float,
                            "(type: float) eta value of the particle")
options.register('outsuffix',
                            '',
                            VarParsing.multiplicity.singleton,
                            VarParsing.varType.string,
                            "outfile suffix to put at the very end before the file extension")
options.register ('nThreads',
		                    1,
		                    VarParsing.multiplicity.singleton,
		                    VarParsing.varType.int,
		                    "Number of threads (also streams)")
options.parseArguments()

#-----------------------------------------
# Input / Output
#-----------------------------------------

format_dict = dict()
format_dict["process"] = options.process
format_dict["E"] = options.E
format_dict["eta"] = options.eta
format_dict["outsuffix"] = options.outsuffix
format_dict["INPUT_DIR"] = os.path.join("/eos/home-k/kyoon", os.environ["DIRDATA_HFNOSE"])
format_dict["OUTPUT_DIR"] = os.path.join("/eos/home-k/kyoon", os.environ["DIRDATA_HFNOSE"], "HFNoseNtuple")

inputfile = "file:{INPUT_DIR}/{process}_2026D60_HFNose/{process}_E{E}_eta{eta}/step3_{process}_E{E}_eta{eta}.root".format(**format_dict)
outputdir = "{OUTPUT_DIR}/{process}".format(**format_dict)
outputfile = outputdir + "/{process}_E{E}_eta{eta}{outsuffix}.root".format(**format_dict)

if not os.path.exists(outputdir):
    os.makedirs(outputdir)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputfile)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(outputfile)
)

#-----------------------------------------
# Main Configuration
#-----------------------------------------

# Max Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Message Logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet (
    wantSummary = cms.untracked.bool(True),
    #numberOfThreads = cms.untracked.uint32(options.nThreads),
    #numberOfStreams = cms.untracked.uint32(options.nThreads)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Process load
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.Geometry.GeometryExtended2026D60_cff')
process.load('Configuration.Geometry.GeometryExtended2026D60Reco_cff')

# GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

# Process
process.HFNoseNtuple = cms.EDAnalyzer('HFNoseNtuplizer',
    InputTag_Tracksters = cms.untracked.InputTag("ticlTrackstersHFNoseEM"),
    withPU = cms.untracked.bool(False)
)

# Timing
process.Timing = cms.Service("Timing",
    summaryOnly = cms.untracked.bool(True),
    useJobReport = cms.untracked.bool(True)
)

process.p = cms.Path(process.HFNoseNtuple)
