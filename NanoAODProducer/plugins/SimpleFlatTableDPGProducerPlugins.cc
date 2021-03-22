#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
typedef SimpleFlatTableProducer<CaloParticle> SimpleCaloParticleFlatTableProducer;

#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
typedef SimpleFlatTableProducer<HGCRecHit> SimpleHGCRecHitFlatTableProducer;

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
typedef SimpleFlatTableProducer<reco::CaloCluster> SimpleCaloClusterFlatTableProducer;

#include "DataFormats/HGCalReco/interface/Trackster.h"
typedef SimpleFlatTableProducer<ticl::Trackster> SimpleTracksterFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SimpleCaloParticleFlatTableProducer);
DEFINE_FWK_MODULE(SimpleHGCRecHitFlatTableProducer);
DEFINE_FWK_MODULE(SimpleCaloClusterFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTracksterFlatTableProducer);
