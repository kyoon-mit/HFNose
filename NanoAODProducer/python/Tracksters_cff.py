import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var,P3Vars

TracksterEMTable = cms.EDProducer("SimpleTracksterFlatTableProducer",
    src = cms.InputTag("ticlTrackstersEM"),
    cut = cms.string(""), 
    name = cms.string("TICLTrackstersEM"),
    doc  = cms.string("TICLTrackstersEM"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        raw_energy = Var('raw_energy', 'float', precision=14, doc='raw energy'),
        raw_em_energy = Var('raw_em_energy', 'float', precision=14, doc='raw energy'),
        time = Var('time', 'float', precision=14, doc='time'),
        time_error = Var('timeError', 'float', precision=14, doc='time error'),
        sigmas = Var('sigmas', 'float', precision=14, doc='raw energy'),
        # regressed_energy = Var('regressed_energy', 'float', precision=14, doc='regressed energy'),
    )
)

HFNoseTracksterEMTable = cms.EDProducer("SimpleTracksterFlatTableProducer",
    src = cms.InputTag("ticlTrackstersHFNoseEM"),
    cut = cms.string(""), 
    name = cms.string("TICLTrackstersHFNoseEM"),
    doc  = cms.string("TICLTrackstersEM in HFNose"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        raw_energy = Var('raw_energy', 'float', precision=14, doc='raw energy'),
        raw_em_energy = Var('raw_em_energy', 'float', precision=14, doc='raw energy'),
        time = Var('time', 'float', precision=14, doc='time'),
        time_error = Var('timeError', 'float', precision=14, doc='time error'),
        sigmas = Var('sigmas', 'float', precision=14, doc='raw energy'),
        # regressed_energy = Var('regressed_energy', 'float', precision=14, doc='regressed energy'),
    )
)

TracksterSequence = cms.Sequence(
    TracksterEMTable
    +HFNoseTracksterEMTable
)

