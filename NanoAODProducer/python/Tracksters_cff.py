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
        raw_em_energy = Var('raw_em_energy', 'float', precision=14, doc='raw em energy'),
        raw_pt = Var('raw_pt', 'float', precision=14, doc='raw pt'),
        raw_em_pt = Var('raw_em_pt', 'float', precision=14, doc='raw em pt'),
        time = Var('time', 'float', precision=14, doc='time'),
        time_error = Var('timeError', 'float', precision=14, doc='time error'),
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
        raw_em_energy = Var('raw_em_energy', 'float', precision=14, doc='raw em energy'),
        raw_pt = Var('raw_pt', 'float', precision=14, doc='raw pt'),
        raw_em_pt = Var('raw_em_pt', 'float', precision=14, doc='raw em pt'),
        time = Var('time', 'float', precision=14, doc='time'),
        time_error = Var('timeError', 'float', precision=14, doc='time error'),
    )
)

TracksterEMSigmaTable = cms.EDProducer("TracksterExtensionTableProducer",
    src = TracksterEMTable.src,
    cut = TracksterEMTable.cut,
    name = TracksterEMTable.name,
    doc  = TracksterEMTable.doc,
    extension = cms.bool(True)
)

HFNoseTracksterEMSigmaTable = cms.EDProducer("TracksterExtensionTableProducer",
    src = HFNoseTracksterEMTable.src,
    cut = HFNoseTracksterEMTable.cut,
    name = HFNoseTracksterEMTable.name,
    doc  = HFNoseTracksterEMTable.doc,
    extension = cms.bool(True)
)

TracksterSequence = cms.Sequence(
    TracksterEMTable
    +HFNoseTracksterEMTable
    +TracksterEMSigmaTable
    +HFNoseTracksterEMSigmaTable
)

