#!/usr/bin/env python2.7

import os
import subprocess

# Set list of pt values (must be string) to iterate over
def set_pt (*args):
    """ Sets list of pt values in string that will be passed on to the generator.
    
    Parameters
    ----------
    *args
        Variable length list of pt values in string.  If any of the values entered is not string, it will raise an exception.
    
    Returns
    -------
    list(str)
        If no value is provided in args, it will return ['.33', '.66', '1', '2', '3', '4', '5', '6', '7', '8', '9'] by default. Otherwise, it will return args in the order the values were provided.
    
    """
    
    if len(args)==0:
        pt_string_list = ['.33', '.66', '1', '2', '3', '4', '5', '6', '7', '8', '9']
        return pt_string_list
    else:
        for arg in args:
            if type(arg)!=str:
                raise Exception("Please provide string only in args.")
        return args


# Set environment variables
def set_env ():
    """ Exports environment variables that are needed in this program.
    
    DIRDATA_HGCNOSE is the directory where the simulation data will be stored. Can be set in config.sh of the top directory. Exception will be raised if it is not set.
    
    Parameters
    ----------
    (none)
    
    Returns
    -------
    (none)
    
    """
        
    if not 'DIRDATA_HGCNOSE' in os.environ:
        raise Exception("Environment variable DIRDATA_HGCNOSE does not exist! This is the top directory where your generated files will be stored. Please set it by editing DIRDATA_HGCNOSE in config.sh and running it.")
    else:
        print "Generated files will be stored in %s" % (os.environ['DIRDATA_HGCNOSE'])

    return None
        
        
def runSteps (pt_string_list, *steps):
    """ Run the step cfg.py files for the step numbers specified.
    
    All config files in the same steps will run in parallel.
    
    Parameters
    ----------
    pt_string_list: list(str)
        List of string-formatted pt values to run cfg.py files over.
    
    *steps: int or str
        Steps to run. Does not have to be ordered; it will order it for you.
        
    Returns
    -------
    (none)
    
    """
    
    dir_run = os.path.abspath(__file__ + '/../run')
    if not os.path.exists(dir_run):
        os.makedirs(dir_run)
        
    steps = sorted(s for s in steps)
    for step in steps:
        dir_cfg = os.path.abspath(dir_run + "/step{}_config".format(step))
        cmd_list= ["cmsRun {0}/step{1}_config/step{1}_2026D60_14TeV_photon_pt{2}_eta35_cfg.py".format(dir_run, step, pt) for pt in pt_string_list]
        bash_command = " & ".join(cmd_list)
        p = subprocess.Popen(bash_command, shell=True)
        p.wait()
        
    return None
        
        
def purge ():
    """ Clear the run directory of all config files.
    
    Parameters
    ----------
    (none)
    
    Returns
    -------
    (none)
    
    """
    
    dir_run = os.path.abspath(__file__ + '/../run')
    subprocess.run("rm {}/*".format(dir_run))
    
    return None


def makeStep1ConfigFiles (pt_string_list, nevents):
    """ Generates CMSSW cfg python scripts to run step 1 generation.
    
    If there are N pt values given, N files will be generated in the sub-directory run/step1_config.
    
    Parameters
    ----------
    pt_string_list : list(str)
        List of string-formatted pt values. Ideally, this will be given through the function set_pt(*args).
        
    nevents : int
        Number of events to generate
        
    Returns
    -------
    None
    """
    
    # Set output directory to put cfg.py files
    dir_run = os.path.abspath(__file__ + '/../run/')
    if not os.path.exists(dir_run + '/step1_config'):
        os.makedirs(dir_run + '/step1_config')
    dir_step1 = dir_run + '/step1_config'
    
    # Set output directory to put simulation root files
    dir_save = os.path.abspath(os.environ['DIRDATA_HGCNOSE'] + '/photon_2026D60')
    for pt in pt_string_list:
        if not os.path.exists(dir_save + '/photon_pt{}'.format(pt)):
            os.makedirs(dir_save + '/photon_pt{}'.format(pt))
    
    # What to write in file
    filedump_preformatted =\
    """
# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: TTbar_14TeV_TuneCP5_cfi --conditions auto:phase2_realistic_T15 -n 10 --era Phase2C10 --eventcontent FEVTDEBUG --relval 9000,100 -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC14TeV --geometry Extended2026D60 --fileout file:step1.root
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C10_cff import Phase2C10

process = cms.Process('SIM',Phase2C10)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D60Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D60_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
# process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedFlat_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Vertex smearing
process.VtxSmeared.MaxX = cms.double(0.0000)
process.VtxSmeared.MinX = cms.double(-0.0000)
process.VtxSmeared.MaxY = cms.double(0.0000)
process.VtxSmeared.MinY = cms.double(-0.0000)
process.VtxSmeared.MaxZ = cms.double(0.0000)
process.VtxSmeared.MinZ = cms.double(-0.0000)

process.MessageLogger = cms.Service("MessageLogger",
        destinations   = cms.untracked.vstring(
                                                'pt{0}_detailedInfo_step1',
                                                'pt{0}_critical_step1'
        ),
        pt{0}_detailedInfo_step1 = cms.untracked.PSet(
                threshold   = cms.untracked.string('DEBUG'),
                default     = cms.untracked.PSet(
                                limit = cms.untracked.int32(10),
                                timespan = cms.untracked.int32(60)
                ),
                WARNING     = cms.untracked.PSet(
                                limit = cms.untracked.int32(100),
                                timespan = cms.untracked.int32(60)
                ),
                ERROR       = cms.untracked.PSet(
                                limit = cms.untracked.int32(100),
                                timespan = cms.untracked.int32(60)
                )
       ),
       pt{0}_critical_step1     = cms.untracked.PSet(
                threshold   = cms.untracked.string('ERROR')
       )
)

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32({1}),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(4),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1_2026D60_14TeV_photon_pt{0}_eta35 nevts:{1}'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:{2}/photon_pt{0}/step1_photon_pt{0}.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

process.generator = cms.EDProducer("FlatRandomPtGunProducer",
    AddAntiParticle = cms.bool(True),
    PGunParameters = cms.PSet(
        AddAntiParticle = cms.bool(True),
        MaxEta = cms.double(3.5001),
        MinEta = cms.double(3.4999),
        MaxPhi = cms.double(3.14159265359),
        MinPhi = cms.double(-3.14159265359),
        MaxPt = cms.double({0} + 0.001),
        MinPt = cms.double({0} - 0.001),
        PartID = cms.vint32(22)
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('single gamma pt {0}')
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.generator)

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
    """
    
    # iterate over string-formatted pt values
    for pt in pt_string_list:
        outfile = dir_step1 + '/step1_2026D60_14TeV_photon_pt{0}_eta35_cfg.py'.format(pt)
        with open(outfile, 'w') as f:
            f.write(filedump_preformatted.format(pt, nevents, dir_save))
    
    return
        

def makeStep2ConfigFiles (pt_string_list, nevents):
    """ Generates CMSSW cfg python scripts to run step 2 digitization.
    
    If there are N pt values given, N files will be generated in the sub-directory run/step2_config.
    
    Parameters
    ----------
    pt_string_list : list(str)
        List of string-formatted pt values. Ideally, this will be given through the function set_pt(*args).
        
    nevents : int
        Number of events to process
        
    Returns
    -------
    None
    """
    
    # Set output directory to put cfg.py files
    dir_run = os.path.abspath(__file__ + '/../run/')
    if not os.path.exists(dir_run + '/step2_config'):
        os.makedirs(dir_run + '/step2_config')
    dir_step2 = dir_run + '/step2_config'
    
    # Set output directory to put simulation root files
    dir_save = os.path.abspath(os.environ['DIRDATA_HGCNOSE'] + '/photon_2026D60')
    if not os.path.exists(dir_save):
        os.makedirs(dir_save)
    
    # What to write in file
    filedump_preformatted =\
    """
# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --conditions auto:phase2_realistic_T15 -s DIGI:pdigi_valid,L1,L1TrackTrigger,DIGI2RAW,HLT:@fake2 --datatier GEN-SIM-DIGI-RAW -n 10 --geometry Extended2026D60 --era Phase2C10 --eventcontent FEVTDEBUGHLT --filein file:step1.root --fileout file:step2.root
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C10_cff import Phase2C10

process = cms.Process('HLT',Phase2C10)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D60Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_Fake2_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger = cms.Service("MessageLogger",
        destinations   = cms.untracked.vstring(
                                                'pt{0}_detailedInfo_step2',
                                                'pt{0}_critical_step2'
        ),
        pt{0}_detailedInfo_step2 = cms.untracked.PSet(
                threshold   = cms.untracked.string('DEBUG'),
                default     = cms.untracked.PSet(
                                limit = cms.untracked.int32(10),
                                timespan = cms.untracked.int32(60)
                ),
                WARNING     = cms.untracked.PSet(
                                limit = cms.untracked.int32(100),
                                timespan = cms.untracked.int32(60)
                ),
                ERROR       = cms.untracked.PSet(
                                limit = cms.untracked.int32(100),
                                timespan = cms.untracked.int32(60)
                )
        ),
        pt{0}_critical_step2     = cms.untracked.PSet(
                threshold   = cms.untracked.string('ERROR')
        )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32({1}),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('file:{2}/photon_pt{0}/step1_photon_pt{0}.root'),
    inputCommands = cms.untracked.vstring(
        'keep *', 
        'drop *_genParticles_*_*', 
        'drop *_genParticlesForJets_*_*', 
        'drop *_kt4GenJets_*_*', 
        'drop *_kt6GenJets_*_*', 
        'drop *_iterativeCone5GenJets_*_*', 
        'drop *_ak4GenJets_*_*', 
        'drop *_ak7GenJets_*_*', 
        'drop *_ak8GenJets_*_*', 
        'drop *_ak4GenJetsNoNu_*_*', 
        'drop *_ak8GenJetsNoNu_*_*', 
        'drop *_genCandidatesForMET_*_*', 
        'drop *_genParticlesForMETAllVisible_*_*', 
        'drop *_genMetCalo_*_*', 
        'drop *_genMetCaloAndNonPrompt_*_*', 
        'drop *_genMetTrue_*_*', 
        'drop *_genMetIC5GenJs_*_*'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(4),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2_2026D60_14TeV_photon_pt{0}_eta35 nevts:{1}'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:{2}/photon_pt{0}/step2_photon_pt{0}.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1TrackTrigger_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.FEVTDEBUGHLToutput_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
    """
    
    # iterate over string-formatted pt values
    for pt in pt_string_list:
        outfile = dir_step2 + '/step2_2026D60_14TeV_photon_pt{0}_eta35_cfg.py'.format(pt)
        with open(outfile, 'w') as f:
            f.write(filedump_preformatted.format(pt, nevents, dir_save))
    
    return


def makeStep3ConfigFiles (pt_string_list, nevents):
    """ Generates CMSSW cfg python scripts to run step 3 reconstruction.
    
    If there are N pt values given, N files will be generated in the sub-directory run/step3_config.
    
    Parameters
    ----------
    pt_string_list : list(str)
        List of string-formatted pt values. Ideally, this will be given through the function set_pt(*args).
        
    nevents : int
        Number of events to process
        
    Returns
    -------
    None
    """
    
    # Set output directory to put cfg.py files
    dir_run = os.path.abspath(__file__ + '/../run/')
    if not os.path.exists(dir_run + '/step3_config'):
        os.makedirs(dir_run + '/step3_config')
    dir_step3 = dir_run + '/step3_config'
    
    # Set output directory to put simulation root files
    dir_save = os.path.abspath(os.environ['DIRDATA_HGCNOSE'] + '/photon_2026D60')
    if not os.path.exists(dir_save):
        os.makedirs(dir_save)
    
    # What to write in file
    filedump_preformatted =\
    """
# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:phase2_realistic_T15 -n 10 --era Phase2C10 --eventcontent FEVTDEBUGHLT,MINIAODSIM,DQM --runUnscheduled -s RAW2DIGI,L1Reco,RECO,RECOSIM,PAT,VALIDATION:@phase2Validation+@miniAODValidation,DQM:@phase2+@miniAODDQM --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --geometry Extended2026D60 --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C10_cff import Phase2C10

process = cms.Process('RECO',Phase2C10)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D60Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger = cms.Service("MessageLogger",
       destinations   = cms.untracked.vstring(
                                                'pt{0}_detailedInfo_step3',
                                                'pt{0}_critical_step3'
       ),
       pt{0}_detailedInfo_step3 = cms.untracked.PSet(
                threshold   = cms.untracked.string('DEBUG'),
                default     = cms.untracked.PSet(
                                limit = cms.untracked.int32(10),
                                timespan = cms.untracked.int32(60)
                ),
                WARNING     = cms.untracked.PSet(
                                limit = cms.untracked.int32(100),
                                timespan = cms.untracked.int32(60)
                ),
                ERROR       = cms.untracked.PSet(
                                limit = cms.untracked.int32(100),
                                timespan = cms.untracked.int32(60)
                )
       ),
       pt{0}_critical_step3     = cms.untracked.PSet(
                threshold   = cms.untracked.string('ERROR')
       )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32({1}),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:{2}/photon_pt{0}/step2_photon_pt{0}.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(4),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3_2026D60_14TeV_photon_pt{0}_eta35 nevts:{1}'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:{2}/photon_pt{0}/step3_photon_pt{0}.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:{2}/photon_pt{0}/step3_photon_pt{0}_inMINIAODSIM.root'),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:{2}/photon_pt{0}/step3_photon_pt{0}_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

# Path and EndPath definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.L1Reco_step = cms.Path(process.L1Reco)
#process.reconstruction_step = cms.Path(process.reconstruction)
#process.recosim_step = cms.Path(process.recosim)
#process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
#process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
#process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
#process.Flag_trkPOGFilters = cms.Path(~process.logErrorTooManyClusters)
#process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
#process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
#process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
#process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
#process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
#process.Flag_eeBadScFilter = cms.Path()
#process.Flag_METFilters = cms.Path(process.metFilters)
#process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
#process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
#process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
#process.Flag_HBHENoiseIsoFilter = cms.Path()
#process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
#process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
#process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
#process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
#process.Flag_HBHENoiseFilter = cms.Path()
#process.Flag_trkPOG_toomanystripclus53X = cms.Path()
#process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
#process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
#process.Flag_trkPOG_manystripclus53X = cms.Path()
#process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
#process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
#process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
#process.Flag_BadPFMuonDzFilter = cms.Path(process.BadPFMuonDzFilter)
#process.prevalidation_step = cms.Path(process.baseCommonPreValidation)
#process.prevalidation_step1 = cms.Path(process.globalPrevalidationTracking)
#process.prevalidation_step2 = cms.Path(process.globalPrevalidationMuons)
#process.prevalidation_step3 = cms.Path(process.globalPrevalidationJetMETOnly)
#process.prevalidation_step4 = cms.Path(process.prebTagSequenceMC)
#process.prevalidation_step5 = cms.Path(process.produceDenoms)
#process.prevalidation_step6 = cms.Path(process.globalPrevalidationHCAL)
#process.prevalidation_step7 = cms.Path(process.globalPrevalidationHGCal)
#process.prevalidation_step8 = cms.Path(process.prevalidationMiniAOD)
#process.validation_step = cms.EndPath(process.baseCommonValidation)
#process.validation_step1 = cms.EndPath(process.globalValidationTrackingOnly)
#process.validation_step2 = cms.EndPath(process.globalValidationMuons)
#process.validation_step3 = cms.EndPath(process.globalValidationJetMETonly)
#process.validation_step4 = cms.EndPath(process.electronValidationSequence)
#process.validation_step5 = cms.EndPath(process.photonValidationSequence)
#process.validation_step6 = cms.EndPath(process.bTagPlotsMCbcl)
#process.validation_step7 = cms.EndPath((process.TauValNumeratorAndDenominatorQCD+process.TauValNumeratorAndDenominatorRealData+process.TauValNumeratorAndDenominatorRealElectronsData+process.TauValNumeratorAndDenominatorRealMuonsData+process.TauValNumeratorAndDenominatorZEE+process.TauValNumeratorAndDenominatorZMM+process.TauValNumeratorAndDenominatorZTT))
#process.validation_step8 = cms.EndPath(process.globalValidationHCAL)
#process.validation_step9 = cms.EndPath(process.globalValidationHGCal)
#process.validation_step10 = cms.EndPath(process.globalValidationMTD)
#process.validation_step11 = cms.EndPath(process.globalValidationOuterTracker)
#process.validation_step12 = cms.EndPath(process.validationECALPhase2)
#process.validation_step13 = cms.EndPath(process.trackerphase2ValidationSource)
#process.validation_step14 = cms.EndPath(process.validationMiniAOD)
#process.dqmoffline_step = cms.EndPath(process.DQMOfflineTracking)
#process.dqmoffline_1_step = cms.EndPath(process.DQMOuterTracker)
#process.dqmoffline_2_step = cms.EndPath(process.DQMOfflineTrackerPhase2)
#process.dqmoffline_3_step = cms.EndPath(process.DQMOfflineMuon)
#process.dqmoffline_4_step = cms.EndPath(process.DQMOfflineHcal)
#process.dqmoffline_5_step = cms.EndPath(process.DQMOfflineHcal2)
#process.dqmoffline_6_step = cms.EndPath(process.DQMOfflineEGamma)
#process.dqmoffline_7_step = cms.EndPath(process.DQMOfflineMiniAOD)
#process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
#process.dqmofflineOnPAT_1_step = cms.EndPath(process.PostDQMOfflineMiniAOD)
#process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
#process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)
#process.DQMoutput_step = cms.EndPath(process.DQMoutput)

process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path()
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path()
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
process.Flag_HBHENoiseFilter = cms.Path()
process.Flag_trkPOG_toomanystripclus53X = cms.Path()
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path()
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.prevalidation_step = cms.Path(process.baseCommonPreValidation)
process.prevalidation_step1 = cms.Path(process.globalPrevalidationTracking)
process.prevalidation_step2 = cms.Path(process.globalPrevalidationMuons)
process.prevalidation_step3 = cms.Path(process.globalPrevalidationJetMETOnly)
process.prevalidation_step4 = cms.Path(process.prebTagSequenceMC)
process.prevalidation_step5 = cms.Path(process.globalPrevalidationHCAL)
process.prevalidation_step6 = cms.Path(process.prevalidationMiniAOD)
process.validation_step = cms.EndPath(process.baseCommonValidation)
process.validation_step1 = cms.EndPath(process.globalValidationTrackingOnly)
process.validation_step2 = cms.EndPath(process.globalValidationMuons)
process.validation_step3 = cms.EndPath(process.globalValidationJetMETonly)
process.validation_step4 = cms.EndPath(process.bTagPlotsMCbcl)
process.validation_step5 = cms.EndPath(process.globalValidationHCAL)
process.validation_step6 = cms.EndPath(process.globalValidationHGCal)
process.validation_step7 = cms.EndPath(process.globalValidationOuterTracker)
process.validation_step8 = cms.EndPath(process.validationMiniAOD)
process.dqmoffline_step = cms.EndPath(process.DQMOfflineTracking)
process.dqmoffline_1_step = cms.EndPath(process.DQMOuterTracker)
process.dqmoffline_2_step = cms.EndPath(process.DQMOfflineMuon)
process.dqmoffline_3_step = cms.EndPath(process.DQMOfflineHcal)
process.dqmoffline_4_step = cms.EndPath(process.HcalDQMOfflineSequence)
process.dqmoffline_5_step = cms.EndPath(process.DQMOfflineEGamma)
process.dqmoffline_6_step = cms.EndPath(process.DQMOfflineMiniAOD)
process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
process.dqmofflineOnPAT_1_step = cms.EndPath(process.PostDQMOfflineMiniAOD)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadPFMuonDzFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.prevalidation_step,process.prevalidation_step1,process.prevalidation_step2,process.prevalidation_step3,process.prevalidation_step4,process.prevalidation_step5,process.prevalidation_step6,process.prevalidation_step7,process.prevalidation_step8,process.validation_step,process.validation_step1,process.validation_step2,process.validation_step3,process.validation_step4,process.validation_step5,process.validation_step6,process.validation_step7,process.validation_step8,process.validation_step9,process.validation_step10,process.validation_step11,process.validation_step12,process.validation_step13,process.validation_step14,process.dqmoffline_step,process.dqmoffline_1_step,process.dqmoffline_2_step,process.dqmoffline_3_step,process.dqmoffline_4_step,process.dqmoffline_5_step,process.dqmoffline_6_step,process.dqmoffline_7_step,process.dqmofflineOnPAT_step,process.dqmofflineOnPAT_1_step,process.FEVTDEBUGHLToutput_step,process.MINIAODSIMoutput_step,process.DQMoutput_step)
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# TICL
from RecoHGCal.TICL.ticl_iterations import TICL_iterations_withReco,TICL_iterations
# process = TICL_iterations_withReco(process)
# process = TICL_iterations(process)

# End of customisation functions

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
    """
    
    # iterate over string-formatted pt values
    for pt in pt_string_list:
        outfile = dir_step3 + '/step3_2026D60_14TeV_photon_pt{0}_eta35_cfg.py'.format(pt)
        with open(outfile, 'w') as f:
            f.write(filedump_preformatted.format(pt, nevents, dir_save))
    
    return


def makeStep4ConfigFiles (pt_string_list, nevents):
    """ Generates CMSSW cfg python scripts to run step 4 harvesting.
    
    If there are N pt values given, N files will be generated in the sub-directory run/step4_config.
    
    Parameters
    ----------
    pt_string_list : list(str)
        List of string-formatted pt values. Ideally, this will be given through the function set_pt(*args).
        
    nevents : int
        Number of events to process
        
    Returns
    -------
    None
    """
    
    # Set output directory to put cfg.py files
    dir_run = os.path.abspath(__file__ + '/../run/')
    if not os.path.exists(dir_run + '/step4_config'):
        os.makedirs(dir_run + '/step4_config')
    dir_step4 = dir_run + '/step4_config'
    
    # Set output directory to put simulation root files
    dir_save = os.path.abspath(os.environ['DIRDATA_HGCNOSE'] + '/photon_2026D60')
    if not os.path.exists(dir_save):
        os.makedirs(dir_save)
    
    # What to write in file
    filedump_preformatted =\
    """
# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step4 --conditions auto:phase2_realistic_T15 -s HARVESTING:@phase2Validation+@phase2+@miniAODValidation+@miniAODDQM --scenario pp --filetype DQM --geometry Extended2026D60 --era Phase2C10 --mc -n {1} --filein file:step3_inDQM.root --fileout file:step4.root
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C10_cff import Phase2C10

process = cms.Process('HARVESTING',Phase2C10)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D60Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.DQMSaverAtRunEnd_cff')
process.load('Configuration.StandardSequences.Harvesting_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32({1}),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("DQMRootSource",
    fileNames = cms.untracked.vstring('file:{2}/photon_pt{0}/step3_photon_pt{0}_inDQM.root')
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step4_2026D60_14TeV_photon_pt{0}_eta35 nevts:{1}'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

# Path and EndPath definitions
process.genHarvesting = cms.Path(process.postValidation_gen)
process.validationprodHarvesting = cms.Path(process.hltpostvalidation_prod+process.postValidation_gen)
process.validationHarvestingHI = cms.Path(process.postValidationHI)
process.dqmHarvestingExtraHLT = cms.Path(process.DQMOffline_SecondStep_ExtraHLT+process.DQMOffline_Certification)
process.alcaHarvesting = cms.Path()
process.validationHarvestingFS = cms.Path(process.recoMuonPostProcessors+process.postValidationTracking+process.MuIsoValPostProcessor+process.calotowersPostProcessor+process.hcalSimHitsPostProcessor+process.hcaldigisPostProcessor+process.hcalrechitsPostProcessor+process.electronPostValidationSequence+process.photonPostProcessor+process.pfJetClient+process.pfMETClient+process.pfJetResClient+process.pfElectronClient+process.rpcRecHitPostValidation_step+process.makeBetterPlots+process.bTagCollectorSequenceMCbcl+process.METPostProcessor+process.L1GenPostProcessor+process.bdHadronTrackPostProcessor+process.MuonGEMHitsPostProcessors+process.MuonGEMDigisPostProcessors+process.MuonGEMRecHitsPostProcessors+process.hgcalPostProcessor+process.MuonME0DigisPostProcessors+process.MuonME0SegPostProcessors+process.trackerphase2ValidationHarvesting+process.postValidation_gen)
process.validationpreprodHarvesting = cms.Path(process.postValidation_preprod+process.hltpostvalidation_preprod+process.postValidation_gen)
process.validationHarvestingNoHLT = cms.Path(process.postValidation+process.postValidation_gen)
process.validationHarvesting = cms.Path(process.postValidation+process.hltpostvalidation+process.postValidation_gen)
process.validationpreprodHarvestingNoHLT = cms.Path(process.postValidation_preprod+process.postValidation_gen)
process.dqmHarvestingPOGMC = cms.Path(process.DQMOffline_SecondStep_PrePOGMC)
process.dqmHarvestingFakeHLT = cms.Path(process.DQMOffline_SecondStep_FakeHLT+process.DQMOffline_Certification)
process.dqmHarvesting = cms.Path(process.DQMOffline_SecondStep+process.DQMOffline_Certification)
process.postValidation_common_step = cms.Path(process.postValidation_common)
process.postValidationTracking_step = cms.Path(process.postValidationTracking)
process.postValidation_muons_step = cms.Path(process.postValidation_muons)
process.postValidation_JetMET_step = cms.Path(process.postValidation_JetMET)
process.electronPostValidationSequence_step = cms.Path(process.electronPostValidationSequence)
process.photonPostProcessor_step = cms.Path(process.photonPostProcessor)
process.bTagCollectorSequenceMCbcl_step = cms.Path(process.bTagCollectorSequenceMCbcl)
process.runTauEff_step = cms.Path(process.runTauEff)
process.postValidation_HCAL_step = cms.Path(process.postValidation_HCAL)
process.hgcalValidatorPostProcessor_step = cms.Path(process.hgcalValidatorPostProcessor)
process.mtdValidationPostProcessor_step = cms.Path(process.mtdValidationPostProcessor)
process.postValidationOuterTracker_step = cms.Path(process.postValidationOuterTracker)
process.trackerphase2ValidationHarvesting_step = cms.Path(process.trackerphase2ValidationHarvesting)
process.DQMHarvestTracking_step = cms.Path(process.DQMHarvestTracking)
process.DQMHarvestOuterTracker_step = cms.Path(process.DQMHarvestOuterTracker)
process.DQMHarvestTrackerPhase2_step = cms.Path(process.DQMHarvestTrackerPhase2)
process.DQMHarvestMuon_step = cms.Path(process.DQMHarvestMuon)
process.DQMCertMuon_step = cms.Path(process.DQMCertMuon)
process.DQMHarvestHcal_step = cms.Path(process.DQMHarvestHcal)
process.DQMHarvestHcal2_step = cms.Path(process.DQMHarvestHcal2)
process.DQMHarvestEGamma_step = cms.Path(process.DQMHarvestEGamma)
process.DQMCertEGamma_step = cms.Path(process.DQMCertEGamma)
process.DQMHarvestMiniAOD_step = cms.Path(process.DQMHarvestMiniAOD)
process.dqmsave_step = cms.Path(process.DQMSaver)

# Schedule definition
process.schedule = cms.Schedule(process.postValidation_common_step,process.postValidationTracking_step,process.postValidation_muons_step,process.postValidation_JetMET_step,process.electronPostValidationSequence_step,process.photonPostProcessor_step,process.bTagCollectorSequenceMCbcl_step,process.runTauEff_step,process.postValidation_HCAL_step,process.hgcalValidatorPostProcessor_step,process.mtdValidationPostProcessor_step,process.postValidationOuterTracker_step,process.trackerphase2ValidationHarvesting_step,process.DQMHarvestTracking_step,process.DQMHarvestOuterTracker_step,process.DQMHarvestTrackerPhase2_step,process.DQMHarvestMuon_step,process.DQMCertMuon_step,process.DQMHarvestHcal_step,process.DQMHarvestHcal2_step,process.DQMHarvestEGamma_step,process.DQMCertEGamma_step,process.validationHarvestingMiniAOD,process.DQMHarvestMiniAOD_step,process.dqmsave_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
    """
    
    # iterate over string-formatted pt values
    for pt in pt_string_list:
        outfile = dir_step4 + '/step4_2026D60_14TeV_photon_pt{0}_eta35_cfg.py'.format(pt)
        with open(outfile, 'w') as f:
            f.write(filedump_preformatted.format(pt, nevents, dir_save))
    
    return
