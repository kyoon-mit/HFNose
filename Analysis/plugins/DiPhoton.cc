#include "DiPhoton.h"

#include <iostream>
#include <cmath> // Switch to TMath.h if you need more physics-related functions
#include "DataFormats/Math/interface/deltaR.h"

#include <TH1.h>
// #include <TH2.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Physics Objects
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

// HGCNose
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

// CMS Coordinate System
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

// Detector Geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"


DiPhoton::DiPhoton ( const edm::ParameterSet& iConfig ) :

    TH1_Container_ (),
    // TH2_Container_(),
    
    // (tag name, default value (label, instance, process) -- CHECK SPELLING!!
    // Truth objects
    tag_CaloParticle_MergedCaloTruth_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_MergedCaloTruth", edm::InputTag ("mix", "MergedCaloTruth") ) ), 
    tag_GenParticle_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_GenParticle", edm::InputTag ("genParticles") ) ),

    // Reco level objects
    tag_HGCHFNoseRecHits_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_HGCHFNoseRecHits", edm::InputTag ("HGCalRecHit", "HGCHFNoseRecHits") ) ),
    tag_HGCalLayerClustersHFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_HGCalLayerClusterHFNose", edm::InputTag ("hgcalLayerClustersHFNose") ) ),
    
    // Pre-selection parameters
    select_PID_ ( 22 ),
    select_EtaLow_ ( 3.49 ),
    select_EtaHigh_ ( 3.51 ),
    select_coneR_ ( 0.5 )
    
{
    // consumes: frequent request of additional data | mayConsume: infrequent
    token_CaloParticle_MergedCaloTruth_ = mayConsume<std::vector<CaloParticle>> ( tag_CaloParticle_MergedCaloTruth_ );
    token_GenParticle_ = consumes<reco::GenParticleCollection> ( tag_GenParticle_ );
    
    token_HGCHFNoseRecHits_ = consumes<HGCRecHitCollection> ( tag_HGCHFNoseRecHits_ );
    token_HGCalLayerClustersHFNose_ = consumes<std::vector<reco::CaloCluster>> ( tag_HGCalLayerClustersHFNose_ );
}


DiPhoton::~DiPhoton ()
{
    // Deconstructor
}


void DiPhoton::analyze ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
    // Get CaloParticles: CaloTruth
    edm::Handle<std::vector<CaloParticle>> handle_CaloParticle_MergedCaloTruth_;
    iEvent.getByToken ( token_CaloParticle_MergedCaloTruth_, handle_CaloParticle_MergedCaloTruth_ );

    // Get GenParticles
    edm::Handle<reco::GenParticleCollection> handle_GenParticle;
    iEvent.getByToken ( token_GenParticle_, handle_GenParticle );

    // Get HGCRecHits
    edm::Handle<HGCRecHitCollection> handle_HGCHFNoseRecHits;
    iEvent.getByToken ( token_HGCHFNoseRecHits_, handle_HGCHFNoseRecHits );

    // Get CaloClusters
    edm::Handle<std::vector<reco::CaloCluster>> handle_HGCalLayerClustersHFNose;
    iEvent.getByToken ( token_HGCalLayerClustersHFNose_, handle_HGCalLayerClustersHFNose ); 

    // Get HGCalGeometry
    edm::ESHandle<HGCalGeometry> handle_HGCalGeometry;
    iSetup.get<IdealGeometryRecord>().get( "HGCalHFNoseSensitive", handle_HGCalGeometry );
    // const HGCalGeometry* geom_HGCal = handle_HGCalGeometry.product();
    
    if ( handle_HGCHFNoseRecHits.isValid() && handle_GenParticle.isValid() && handle_HGCalLayerClustersHFNose.isValid() )
    {
//        const HGCRecHitCollection HGCRecHits = *handle_HGCRecHits.product();
//        const HGCalGeometry HGCalGeom = *handle_HGCalGeometry.product();
//        const std::vector<reco::CaloCluster> HFNoseClusters = *handle_HGCalLayerClustersHFNose.product();
    }
    else std::cout << "Handle(s) invalid!" << std::endl;
}


std::vector<math::XYZTLorentzVectorF> DiPhoton::getGenTruthP4 ( const reco::GenParticleCollection & GenParticles )
{
    std::vector<math::XYZTLorentzVectorF> truth_container;
    
    for ( auto const& gen: GenParticles )
    {
        if ( gen.pdgId() == select_PID_
                && abs(gen.eta()) > select_EtaLow_
                && abs(gen.eta()) < select_EtaHigh_  )
        {
            truth_container.emplace_back ( (math::XYZTLorentzVectorF) gen.p4() );
        }
    }
    
    return truth_container;
}


std::vector<math::XYZTLorentzVectorF> DiPhoton::getCaloTruthP4 ( const std::vector<CaloParticle> & CaloParticles )
{
    std::vector<math::XYZTLorentzVectorF> truth_container;
    
    for ( auto const& clp: CaloParticles )
    {
        if ( clp.pdgId() == select_PID_
                && abs(clp.eta()) > select_EtaLow_
                && abs(clp.eta()) < select_EtaHigh_  )
        {
            truth_container.emplace_back ( (math::XYZTLorentzVectorF) clp.p4() );
        }
    }
    
    return truth_container;
}


void DiPhoton::beginJob ()
{
    edm::Service<TFileService> fs;
    
    TH1_Container_["truthEDist"] = fs->make<TH1F>("truthEDist", "Truth Energy Distribution", 200, 0, 200);
}


void DiPhoton::endJob ()
{
    // Job ends
    std::cout << "The job, DiPhoton, has ended. Thank you for your patience." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (DiPhoton);
