import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var,P3Vars

HGCEERecHitsTable = cms.EDProducer("SimpleHGCRecHitFlatTableProducer",
    src = cms.InputTag("HGCalRecHit:HGCEERecHits"),
    cut = cms.string(""), 
    name = cms.string("RecHitsHGCEE"),
    doc  = cms.string("RecHits in HGC EE"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        detId = Var('detid().rawId()', 'int', precision=-1, doc='detId'),
        flags = Var('flags', 'int', precision=-1, doc='flag'),
        time = Var('time', 'float', precision=14, doc='hit time'),
        timeError = Var('timeError', 'float', precision=14, doc='time error'),
        energy = Var('energy', 'float', precision=14, doc='energy'),
        isTimeValid = Var('isTimeValid', 'bool', precision=-1, doc='(bool) is time valid'),
        isTimeErrorValid = Var('isTimeErrorValid', 'bool', precision=-1, doc='(bool) is time error valid'),
        outOfTimeEnergy = Var('outOfTimeEnergy', 'float', precision=14, doc='(only for out of time events) out of time energy'),
        chi2 = Var('chi2', 'float', precision=14, doc='chi2'),
        outOfTimeChi2 = Var('outOfTimeChi2', 'float', precision=14, doc='(only for out of time events) out of time chi2'),
        signalOverSigmaNoise = Var('signalOverSigmaNoise', 'float', precision=14, doc='signal over sigma noise')
    )
)

HGCHEFRecHitsTable = cms.EDProducer("SimpleHGCRecHitFlatTableProducer",
    src = cms.InputTag("HGCalRecHit:HGCHEFRecHits"),
    cut = cms.string(""), 
    name = cms.string("RecHitsHGCHEF"),
    doc  = cms.string("RecHits in HGC HEF"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        detId = Var('detid().rawId()', 'int', precision=-1, doc='detId'),
        flags = Var('flags', 'int', precision=-1, doc='flag'),
        time = Var('time', 'float', precision=14, doc='hit time'),
        timeError = Var('timeError', 'float', precision=14, doc='time error'),
        energy = Var('energy', 'float', precision=14, doc='energy'),
        isTimeValid = Var('isTimeValid', 'bool', precision=-1, doc='(bool) is time valid'),
        isTimeErrorValid = Var('isTimeErrorValid', 'bool', precision=-1, doc='(bool) is time error valid'),
        outOfTimeEnergy = Var('outOfTimeEnergy', 'float', precision=14, doc='(only for out of time events) out of time energy'),
        chi2 = Var('chi2', 'float', precision=14, doc='chi2'),
        outOfTimeChi2 = Var('outOfTimeChi2', 'float', precision=14, doc='(only for out of time events) out of time chi2'),
        signalOverSigmaNoise = Var('signalOverSigmaNoise', 'float', precision=14, doc='signal over sigma noise')
    )
)

HGCHEBRecHitsTable = cms.EDProducer("SimpleHGCRecHitFlatTableProducer",
    src = cms.InputTag("HGCalRecHit:HGCHEBRecHits"),
    cut = cms.string(""), 
    name = cms.string("RecHitsHGCHEB"),
    doc  = cms.string("RecHits in HGC HEB"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        detId = Var('detid().rawId()', 'int', precision=-1, doc='detId'),
        flags = Var('flags', 'int', precision=-1, doc='flag'),
        time = Var('time', 'float', precision=14, doc='hit time'),
        timeError = Var('timeError', 'float', precision=14, doc='time error'),
        energy = Var('energy', 'float', precision=14, doc='energy'),
        isTimeValid = Var('isTimeValid', 'bool', precision=-1, doc='(bool) is time valid'),
        isTimeErrorValid = Var('isTimeErrorValid', 'bool', precision=-1, doc='(bool) is time error valid'),
        outOfTimeEnergy = Var('outOfTimeEnergy', 'float', precision=14, doc='(only for out of time events) out of time energy'),
        chi2 = Var('chi2', 'float', precision=14, doc='chi2'),
        outOfTimeChi2 = Var('outOfTimeChi2', 'float', precision=14, doc='(only for out of time events) out of time chi2'),
        signalOverSigmaNoise = Var('signalOverSigmaNoise', 'float', precision=14, doc='signal over sigma noise')
    )
)

HFNoseRecHitsTable = cms.EDProducer("SimpleHGCRecHitFlatTableProducer",
    src = cms.InputTag("HGCalRecHit:HGCHFNoseRecHits"),
    cut = cms.string(""), 
    name = cms.string("RecHitsHFNose"),
    doc  = cms.string("RecHits in HFNose"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        detId = Var('detid().rawId()', 'int', precision=-1, doc='detId'),
        flags = Var('flags', 'int', precision=-1, doc='flag'),
        time = Var('time', 'float', precision=14, doc='hit time'),
        timeError = Var('timeError', 'float', precision=14, doc='time error'),
        energy = Var('energy', 'float', precision=14, doc='energy'),
        isTimeValid = Var('isTimeValid', 'bool', precision=-1, doc='(bool) is time valid'),
        isTimeErrorValid = Var('isTimeErrorValid', 'bool', precision=-1, doc='(bool) is time error valid'),
        outOfTimeEnergy = Var('outOfTimeEnergy', 'float', precision=14, doc='(only for out of time events) out of time energy'),
        chi2 = Var('chi2', 'float', precision=14, doc='chi2'),
        outOfTimeChi2 = Var('outOfTimeChi2', 'float', precision=14, doc='(only for out of time events) out of time chi2'),
        signalOverSigmaNoise = Var('signalOverSigmaNoise', 'float', precision=14, doc='signal over sigma noise')
    )
)

HGCEERecHitsPositionTable = cms.EDProducer("HGCRecHitPositionFromDetIDTableProducer",
    src = HGCEERecHitsTable.src,
    cut = HGCEERecHitsTable.cut,
    name = HGCEERecHitsTable.name,
    doc  = HGCEERecHitsTable.doc,
)

HGCHEFRecHitsPositionTable = cms.EDProducer("HGCRecHitPositionFromDetIDTableProducer",
    src = HGCHEFRecHitsTable.src,
    cut = HGCHEFRecHitsTable.cut, 
    name = HGCHEFRecHitsTable.name,
    doc  = HGCHEFRecHitsTable.doc,
)

HGCHEBRecHitsPositionTable = cms.EDProducer("HGCRecHitPositionFromDetIDTableProducer",
    src = HGCHEBRecHitsTable.src,
    cut = HGCHEBRecHitsTable.cut, 
    name = HGCHEBRecHitsTable.name,
    doc  = HGCHEBRecHitsTable.doc,
)

HFNoseRecHitsPositionTable = cms.EDProducer("HGCRecHitPositionFromDetIDTableProducer",
    src = HFNoseRecHitsTable.src,
    cut = HFNoseRecHitsTable.cut, 
    name = HFNoseRecHitsTable.name,
    doc  = HFNoseRecHitsTable.doc,
)

HGCRecHitsSequence = cms.Sequence(
    HGCEERecHitsTable
    +HGCHEFRecHitsTable
    +HGCHEBRecHitsTable
    +HFNoseRecHitsTable
    +HGCEERecHitsPositionTable
    +HGCHEFRecHitsPositionTable
    +HGCHEBRecHitsPositionTable
    +HFNoseRecHitsPositionTable
)

