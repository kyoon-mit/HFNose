import os
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.Eras.Era_Phase2C10_cff import Phase2C10

process = cms.Process('NanoAODOutputAging',Phase2C10)

# Varparsing
options = VarParsing('analysis')
options.register('process', 'photon',
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.string,
                    "process name (e.g. muon)")
options.register('geometry', '2026D60',
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.string,
                    "geometry (e.g. 2026D60)")
options.register('detector', 'HGCNose',
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.string,
                    "detector (e.g. HGCNose)")
options.register('suffix', '',
                    VarParsing.multiplicity.singleton, 
                    VarParsing.varType.string,
                    "inpath suffix")
options.register('E', 200,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.int,
                    "(type: int) E value of the particle")
options.register('eta', 3.5,
                    VarParsing.multiplicity.singleton,
                    VarParsing.varType.float,
                    "(type: float) eta value of the particle")
options.parseArguments()

# Set output and input paths
inpath = options.process + '_' + options.geometry + '_' + options.detector + options.suffix
outpath = options.process + '_' + options.geometry

INPUT_DIR = os.path.abspath(os.environ['DIRDATA_HGCNOSE'] + '/' + inpath + '/')
OUTPUT_DIR = os.path.abspath(os.environ['DIRDATA_HGCNOSE'] + '/NanoAOD/' + outpath + '/TICL/')

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

inputfile = 'file:' + INPUT_DIR + '/{0}_E{1}_eta{2}/step3_{0}_E{1}_eta{2}{3}.root'.format(options.process, options.E, options.eta, options.suffix)
outputfile = OUTPUT_DIR + '/{0}_{1}_E{2}_eta{3}{4}.root'.format(options.process, options.detector, options.E, options.eta, options.suffix)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet (
    wantSummary = cms.untracked.bool(False),
    numberOfThreads = cms.untracked.uint32(4),
    numberOfStreams = cms.untracked.uint32(4)
)

process.load('Configuration.Geometry.GeometryExtended2026D60Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D60_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#-----------------------------------------
# INPUT / OUTPUT
#-----------------------------------------

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputfile)
)

process.output = cms.OutputModule("NanoAODOutputModule",
     compressionAlgorithm = cms.untracked.string('LZMA'),
     compressionLevel = cms.untracked.int32(9),
     dataset = cms.untracked.PSet(
         dataTier = cms.untracked.string('NANOAODSIM'),
         filterName = cms.untracked.string('')
     ),
     fileName = cms.untracked.string(outputfile),
     outputCommands = cms.untracked.vstring(
         'drop *',
         "keep nanoaodFlatTable_*Table_*_*",     # event data
     )
 )

#-----------------------------------------
# HGCstuff
#-----------------------------------------

####
process.load("HGCNose.NanoAODProducer.HGCNanoAOD_TICL_cff")

process.hgc = cms.Path(process.HGCNanoAODSequence)

process.finalize = cms.EndPath(process.output)

process.schedule = cms.Schedule(
    process.hgc,
    process.finalize
)
