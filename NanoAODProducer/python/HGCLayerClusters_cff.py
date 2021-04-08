import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import P3Vars,Var

LayerClusterTable = cms.EDProducer("SimpleCaloClusterFlatTableProducer",
    src = cms.InputTag("hgcalLayerClusters"),
    cut = cms.string(""),
    name = cms.string("LayerClusters"),
    doc  = cms.string("LayerClusters information"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(        
        nHits = Var('hitsAndFractions().size()', 'int', precision=-1, doc='number of hits in a layerCluster'),
        sumHitEnergy = Var('energy', 'float', precision=14, doc='total energy of simhits'),
        eta  = Var("eta()",  float, precision=14),
        phi = Var("phi()", float, precision=14),
    )
)

HFNoseLayerClusterTable = cms.EDProducer("SimpleCaloClusterFlatTableProducer",
    src = cms.InputTag("hgcalLayerClustersHFNose"),
    cut = cms.string(""),
    name = cms.string("LayerClustersHFNose"),
    doc  = cms.string("LayerClusters information"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(        
        nHits = Var('hitsAndFractions().size()', 'int', precision=-1, doc='number of hits in a layerCluster'),
        sumHitEnergy = Var('energy', 'float', precision=14, doc='total energy of simhits'),
        eta  = Var("eta()",  float, precision=14),
        phi = Var("phi()", float, precision=14),
    )
)

LayerClusterLayerNumTable = cms.EDProducer("LayerClusterExtensionProducer",
    src = LayerClusterTable.src,
    cut = LayerClusterTable.cut,
    name = LayerClusterTable.name,
    doc  = LayerClusterTable.doc,
)

HFNoseLayerClusterLayerNumTable = cms.EDProducer("LayerClusterExtensionProducer",
    src = HFNoseLayerClusterTable.src,
    cut = HFNoseLayerClusterTable.cut,
    name = HFNoseLayerClusterTable.name,
    doc  = HFNoseLayerClusterTable.doc,
)

HGCLayerClusterSequence = cms.Sequence(
    LayerClusterTable
    +HFNoseLayerClusterTable
    +LayerClusterLayerNumTable
    +HFNoseLayerClusterLayerNumTable
)
