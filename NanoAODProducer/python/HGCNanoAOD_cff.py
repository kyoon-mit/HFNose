from __future__ import print_function
import FWCore.ParameterSet.Config as cms

###
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import genParticleTable
from PhysicsTools.NanoAOD.genVertex_cff import *
###
from CaloParticles_cff import *
from HGCRecHits_cff import *
from HGCLayerClusters_cff import *
from Tracksters_cff import *

nanoMetadata = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("untagged"),
    )
)

genParticleTable.src = "genParticles"
genParticleTable.variables = cms.PSet(genParticleTable.variables,
    charge = CandVars.charge)

HGCNanoAODSequence = cms.Sequence(nanoMetadata+genVertexTables+genParticleTable
                                  +CaloParticleTable
                                  +HGCRecHitsSequence
                                  +HGCLayerClusterSequence
                                  +TracksterSequence
)
