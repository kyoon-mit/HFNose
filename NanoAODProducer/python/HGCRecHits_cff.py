import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var,P3Vars

HGCEERecHitsTable = cms.EDProducer("SimpleCaloRecHitFlatTableProducer",
    src = cms.InputTag("HGCalRecHit:HGCEERecHits"),
    cut = cms.string(""), 
    name = cms.string("RecHitsHGCEE"),
    doc  = cms.string("RecHits in HGC EE"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        detId = Var('detid().rawId()', 'int', precision=-1, doc='detId'),
        energy = Var('energy', 'float', precision=14, doc='energy'),
        time = Var('time', 'float', precision=14, doc='hit time'),
    )
)

HGCHEFRecHitsTable = cms.EDProducer("SimpleCaloRecHitFlatTableProducer",
    src = cms.InputTag("HGCalRecHit:HGCHEFRecHits"),
    cut = cms.string(""), 
    name = cms.string("RecHitsHGCHEF"),
    doc  = cms.string("RecHits in HGC HEF"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        detId = Var('detid().rawId()', 'int', precision=-1, doc='detId'),
        energy = Var('energy', 'float', precision=14, doc='energy'),
        time = Var('time', 'float', precision=14, doc='hit time'),
    )
)

HGCHEBRecHitsTable = cms.EDProducer("SimpleCaloRecHitFlatTableProducer",
    src = cms.InputTag("HGCalRecHit:HGCHEBRecHits"),
    cut = cms.string(""), 
    name = cms.string("RecHitsHGCHEB"),
    doc  = cms.string("RecHits in HGC HEB"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        detId = Var('detid().rawId()', 'int', precision=-1, doc='detId'),
        energy = Var('energy', 'float', precision=14, doc='energy'),
        time = Var('time', 'float', precision=14, doc='hit time'),
    )
)

HFNoseRecHitsTable = cms.EDProducer("SimpleCaloRecHitFlatTableProducer",
    src = cms.InputTag("HGCalRecHit:HGCHFNoseRecHits"),
    cut = cms.string(""), 
    name = cms.string("RecHitsHFNose"),
    doc  = cms.string("RecHits in HFNose"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        detId = Var('detid().rawId()', 'int', precision=-1, doc='detId'),
        energy = Var('energy', 'float', precision=14, doc='energy'),
        time = Var('time', 'float', precision=14, doc='hit time'),
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

