#include "EnergyResolution.h"

#include <iostream>
#include <array>

#include <cmath> // Switch to TMath.h if you need more physics-related functions
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Physics objects
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// HFNose (forward + HGCal)
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"

// HF (forward + HCal)
// #include "DataFormats/HCalRecHit/interface/HCalRecHitCollections.h"
// #include "DataFormats/HCalDetId/interface/HCalDetId.h"

// CMS Coordinate System
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

// Detector Geometry
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

// ROOT headers

using namespace edm;

EnergyResolution::EnergyResolution ( const edm::ParameterSet& iConfig ) :

    histContainer_ (),
    
    // (tag name, default value (label, instance, process) -- CHECK SPELLING!!!!!!!
    tag_HGCHFNoseRecHits_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_HGCHFNoseRecHits", edm::InputTag ("HGCalRecHit", "HGCHFNoseRecHits") ) ),
    tag_HGCalLayerClustersHFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_HGCalLayerClusterHFNose", edm::InputTag ("hgcalLayerClustersHFNose") ) ),
    tag_GenParticle_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_GenParticle", edm::InputTag ("genParticles") ) ),
    // tag_CaloParticle_MergedCaloTruth_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_MergedCaloTruth", edm::InputTag ("mix", "MergedCaloTruth") ) ), 
    
    // Pre-selection parameters
    select_PID_ ( 22 ),
    select_EtaLow_ ( 3.49 ),
    select_EtaHigh_ ( 3.51 ),
    select_coneR_ ( 0.5 )
    
{
    // consumes: frequent request of additional data | mayConsume: infrequent
    token_HGCRecHits_ = consumes<HGCRecHitCollection> ( tag_HGCHFNoseRecHits_ );
    token_HGCalLayerClustersHFNose_ = consumes<std::vector<reco::CaloCluster>> ( tag_HGCalLayerClustersHFNose_ );
    token_GenParticle_ = consumes<reco::GenParticleCollection> ( tag_GenParticle_ );
    // token_CaloParticle_MergedCaloTruth_ = mayConsume<std::vector<CaloParticle>> ( tag_CaloParticle_MergedCaloTruth_ );
}


EnergyResolution::~EnergyResolution ()
{
    // Deconstructor
}


void EnergyResolution::analyze ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
    // Get hits
    edm::Handle<HGCRecHitCollection> handle_HGCRecHits;
    iEvent.getByToken ( token_HGCRecHits_, handle_HGCRecHits );

    // Get Geometry
    edm::ESHandle<HGCalGeometry> handle_HGCalGeometry;
    iSetup.get<IdealGeometryRecord>().get( "HGCalHFNoseSensitive", handle_HGCalGeometry );
    // const HGCalGeometry* geom_HGCal = handle_HGCalGeometry.product();
    
    edm::Handle<std::vector<reco::CaloCluster>> handle_HGCalLayerClustersHFNose;
    iEvent.getByToken ( token_HGCalLayerClustersHFNose_, handle_HGCalLayerClustersHFNose ); 

    // Get MC truth (maybe consider using a separate function if code becomes too long)
    // edm::Handle<std::vector<CaloParticle>> handle_CaloParticle_MergedCaloTruth_;
    // iEvent.getByToken ( token_CaloParticle_MergedCaloTruth_, handle_CaloParticle_MergedCaloTruth_ );
    edm::Handle<reco::GenParticleCollection> handle_GenParticle;
    iEvent.getByToken ( token_GenParticle_, handle_GenParticle );
    
    if ( handle_HGCRecHits.isValid() && handle_GenParticle.isValid() && handle_HGCalLayerClustersHFNose.isValid() )
    {
//        const HGCRecHitCollection HGCRecHits = *handle_HGCRecHits.product();
//        const HGCalGeometry HGCalGeom = *handle_HGCalGeometry.product();
//        const std::vector<reco::CaloCluster> HFNoseClusters = *handle_HGCalLayerClustersHFNose.product();
    
        const std::vector<math::XYZTLorentzVectorF> truth_container = getTruthP4 ( *handle_GenParticle.product() );
        
        for ( auto const& truth: truth_container)
        {
            fillHist_HGCalRecHitsEnergy_coneR ( truth, *handle_HGCRecHits.product(), handle_HGCalGeometry.product() );
            fillHist_CaloClustersEnergy_coneR ( truth, *handle_HGCalLayerClustersHFNose.product() );
            histContainer_["truthEDist"]->Fill ( truth.energy() );
        }
    }
    else std::cout << "Handle(s) invalid!" << std::endl;
}


std::vector<math::XYZTLorentzVectorF> EnergyResolution::getTruthP4 ( const reco::GenParticleCollection & GenParticles )
{
    std::vector<math::XYZTLorentzVectorF> container;
    
    for ( auto const& gen: GenParticles )
    {
        if ( gen.pdgId() == select_PID_
                && abs(gen.eta()) > select_EtaLow_
                && abs(gen.eta()) < select_EtaHigh_  )
        {
            container.push_back ( (math::XYZTLorentzVectorF) gen.p4() );
        }
    }
    
    return container;
}


void EnergyResolution::fillHist_HGCalRecHitsEnergy_coneR ( const math::XYZTLorentzVectorF & truth, const HGCRecHitCollection & hits, const HGCalGeometry * geom )
{
// Sum hit energy within coneR of truth and fill histograms
// Also does this for individual silicon layer in HGCal detector

    Float_t sum_E = 0;
    Float_t sum_E_layer[8] = {0};
    
    for ( auto const& hit : hits )
    {
        // Hit-related vars
        uint32_t hit_id = hit.id();
        
        // std::shared_ptr<const CaloCellGeometry> thisCell = geom->getGeometry(hit_DetId);
        
        // Cone reconstruction
        const GlobalPoint & hit_globalPosition = geom->getPosition(hit_id);
        Float_t dR = reco::deltaR ( hit_globalPosition.eta(), hit_globalPosition.phi(), truth.eta(), truth.phi() );
        
        if ( dR < select_coneR_ ) // hit within cone of truth
        {
            Float_t hit_energy = hit.energy();
            HFNoseDetId hit_DetId = HFNoseDetId ( hit_id );

            sum_E += hit_energy;
            sum_E_layer[hit_DetId.layer()-1] += hit_energy;
        }
    }
            
    histContainer_["EDist_hits"]->Fill ( sum_E );
    histContainer_["EDist_hits_layer1"]->Fill ( sum_E_layer[0] );
    histContainer_["EDist_hits_layer2"]->Fill ( sum_E_layer[1] );
    histContainer_["EDist_hits_layer3"]->Fill ( sum_E_layer[2] );
    histContainer_["EDist_hits_layer4"]->Fill ( sum_E_layer[3] );
    histContainer_["EDist_hits_layer5"]->Fill ( sum_E_layer[4] );
    histContainer_["EDist_hits_layer6"]->Fill ( sum_E_layer[5] );
    histContainer_["EDist_hits_layer7"]->Fill ( sum_E_layer[6] );
    histContainer_["EDist_hits_layer8"]->Fill ( sum_E_layer[7] );
    
        // Int_t caloTruth_count_validation = 0; // Must be 1 in this event
    // math::XYZTLorentzVectorF this_truth_p4;
    
    /*
    if ( handle_CaloParticle_MergedCaloTruth_.isValid() )
    {
        for ( auto const& caloTruth : *handle_CaloParticle_MergedCaloTruth_.product() )
        {   
            if ( caloTruth.pdgId() == select_PID_ )
            { // Nested to save a bit of calculation
                if ( ( select_EtaLow_ < std::abs(caloTruth.eta()) ) &&
                     ( std::abs(caloTruth.eta()) < select_EtaHigh_ ) )
                {
                    caloTruth_count_validation++;
                    if ( caloTruth_count_validation > 1 )
                    {
                        std::cout << "CaloTruth count 2; skipping event" << std::endl;
                        return;
                    }
                    
                    this_truth_p4 = caloTruth.p4();
                    
                }
            }
        }
        
    }
    else std::cout << "Handle for CaloParticle:MergeCaloTruth invalid!" << std::endl;
    
    
    if ( caloTruth_count_validation == 0 )
    {
        std::cout << "CaloTruth count 0; skipping event" << std::endl;
        return;
    }
    */

}


void EnergyResolution::fillHist_CaloClustersEnergy_coneR ( const math::XYZTLorentzVectorF & truth, const std::vector<reco::CaloCluster> & Clusters )
{
    Float_t sum_E = 0;
    Float_t sum_E_layer[8] = {0};
    Float_t num_layer[8] = {0};
    //auto sum_E_vector = math::PtEtaPhiELorentzVectorF (0., 0., 0., 0.);
    
    for ( auto const& cl : Clusters )
    {
        // Get Layer # (also, get how many total clusters in layer)
        HFNoseDetId cl_DetId = HFNoseDetId( cl.hitsAndFractions().at(0).first );
        Int_t layer = cl_DetId.layer();
        num_layer[layer-1]++;
        
        // Cone reconstruction
        Float_t dR = reco::deltaR ( cl.eta(), cl.phi(), truth.eta(), truth.phi() );
        
        if ( dR < select_coneR_ ) // clusters within cone of truth
        {
            Float_t cl_energy = cl.energy();

            sum_E += cl_energy;
            sum_E_layer[layer-1] += cl_energy;
            //sum_E_vector += math::PtEtaPhiELorentzVectorF ( cl.eta(), cl.phi(), cl_energy );
        }
    }
            
    histContainer_["EDist_clusters_scalar_sum"]->Fill ( sum_E );
    //histContainer_["EDist_clusters_vector_sum"]->Fill ( sum_E_vector.energy() );
    histContainer_["EDist_clusters_layer1"]->Fill( sum_E_layer[0] );
    histContainer_["EDist_clusters_layer2"]->Fill( sum_E_layer[1] );
    histContainer_["EDist_clusters_layer3"]->Fill( sum_E_layer[2] );
    histContainer_["EDist_clusters_layer4"]->Fill( sum_E_layer[3] );
    histContainer_["EDist_clusters_layer5"]->Fill( sum_E_layer[4] );
    histContainer_["EDist_clusters_layer6"]->Fill( sum_E_layer[5] );
    histContainer_["EDist_clusters_layer7"]->Fill( sum_E_layer[6] );
    histContainer_["EDist_clusters_layer8"]->Fill( sum_E_layer[7] );
    histContainer_["num_clusters_layer1"]->Fill( num_layer[0] );
    histContainer_["num_clusters_layer2"]->Fill( num_layer[1] );
    histContainer_["num_clusters_layer3"]->Fill( num_layer[2] );
    histContainer_["num_clusters_layer4"]->Fill( num_layer[3] );
    histContainer_["num_clusters_layer5"]->Fill( num_layer[4] );
    histContainer_["num_clusters_layer6"]->Fill( num_layer[5] );
    histContainer_["num_clusters_layer7"]->Fill( num_layer[6] );
    histContainer_["num_clusters_layer8"]->Fill( num_layer[7] );
}


void EnergyResolution::beginJob ()
{
    edm::Service<TFileService> fs;
    
    histContainer_["truthEDist"] = fs->make<TH1F>("truthEDist", "Truth Energy Distribution", 200, 0, 200);
    
    histContainer_["EDist_hits"] = fs->make<TH1F>("EDist_hits", "Energy Distribution (hits)", 200, 0, 200);
    histContainer_["EDist_hits_layer1"] = fs->make<TH1F>("EDist_hits_layer1", "Energy Distribution Layer 1 (hits)", 200, 0, 200);
    histContainer_["EDist_hits_layer2"] = fs->make<TH1F>("EDist_hits_layer2", "Energy Distribution Layer 2 (hits)", 200, 0, 200);
    histContainer_["EDist_hits_layer3"] = fs->make<TH1F>("EDist_hits_layer3", "Energy Distribution Layer 3 (hits)", 200, 0, 200);
    histContainer_["EDist_hits_layer4"] = fs->make<TH1F>("EDist_hits_layer4", "Energy Distribution Layer 4 (hits)", 200, 0, 200);
    histContainer_["EDist_hits_layer5"] = fs->make<TH1F>("EDist_hits_layer5", "Energy Distribution Layer 5 (hits)", 200, 0, 200);
    histContainer_["EDist_hits_layer6"] = fs->make<TH1F>("EDist_hits_layer6", "Energy Distribution Layer 6 (hits)", 200, 0, 200);
    histContainer_["EDist_hits_layer7"] = fs->make<TH1F>("EDist_hits_layer7", "Energy Distribution Layer 7 (hits)", 200, 0, 200);
    histContainer_["EDist_hits_layer8"] = fs->make<TH1F>("EDist_hits_layer8", "Energy Distribution Layer 8 (hits)", 200, 0, 200);

    histContainer_["EDist_clusters_scalar_sum"] = fs->make<TH1F>("EDist_clusters_scalar_sum", "Energy Distribution (clusters, scaler sum)", 200, 0, 200);
    //histContainer_["EDist_clusters_vector_sum"] = fs->make<TH1F>("EDist_clusters_vector_sum", "Energy Distribution (clusters--vector sum)", 200, 0, 200);
    histContainer_["EDist_clusters_layer1"] = fs->make<TH1F>("EDist_clusters_layer1", "Energy Distribution Layer 1 (clusters)", 200, 0, 200);
    histContainer_["EDist_clusters_layer2"] = fs->make<TH1F>("EDist_clusters_layer2", "Energy Distribution Layer 2 (clusters)", 200, 0, 200);
    histContainer_["EDist_clusters_layer3"] = fs->make<TH1F>("EDist_clusters_layer3", "Energy Distribution Layer 3 (clusters)", 200, 0, 200);
    histContainer_["EDist_clusters_layer4"] = fs->make<TH1F>("EDist_clusters_layer4", "Energy Distribution Layer 4 (clusters)", 200, 0, 200);
    histContainer_["EDist_clusters_layer5"] = fs->make<TH1F>("EDist_clusters_layer5", "Energy Distribution Layer 5 (clusters)", 200, 0, 200);
    histContainer_["EDist_clusters_layer6"] = fs->make<TH1F>("EDist_clusters_layer6", "Energy Distribution Layer 6 (clusters)", 200, 0, 200);
    histContainer_["EDist_clusters_layer7"] = fs->make<TH1F>("EDist_clusters_layer7", "Energy Distribution Layer 7 (clusters)", 200, 0, 200);
    histContainer_["EDist_clusters_layer8"] = fs->make<TH1F>("EDist_clusters_layer8", "Energy Distribution Layer 8 (clusters)", 200, 0, 200);
    
    histContainer_["num_clusters_layer1"] = fs->make<TH1F>("num_clusters_layer1", "Number of Clusters Layer 1", 30, 0, 30);
    histContainer_["num_clusters_layer2"] = fs->make<TH1F>("num_clusters_layer2", "Number of Clusters Layer 2", 30, 0, 30);
    histContainer_["num_clusters_layer3"] = fs->make<TH1F>("num_clusters_layer3", "Number of Clusters Layer 3", 30, 0, 30);
    histContainer_["num_clusters_layer4"] = fs->make<TH1F>("num_clusters_layer4", "Number of Clusters Layer 4", 30, 0, 30);
    histContainer_["num_clusters_layer5"] = fs->make<TH1F>("num_clusters_layer5", "Number of Clusters Layer 5", 30, 0, 30);
    histContainer_["num_clusters_layer6"] = fs->make<TH1F>("num_clusters_layer6", "Number of Clusters Layer 6", 30, 0, 30);
    histContainer_["num_clusters_layer7"] = fs->make<TH1F>("num_clusters_layer7", "Number of Clusters Layer 7", 30, 0, 30);
    histContainer_["num_clusters_layer8"] = fs->make<TH1F>("num_clusters_layer8", "Number of Clusters Layer 8", 30, 0, 30);
}

void EnergyResolution::endJob ()
{
    // THIS SECTION USED TO WORK JUST FINE. WHY IS IT NOW GIVING A SEG FAULT?????
    /*
    const Float_t det_E_Mean        = histContainer_["EDist_hits`"]->GetMean();
    const Float_t det_E_StdDev      = histContainer_["EDist_hits"]->GetStdDev();
    const Float_t truth_E_Mean      = histContainer_["truthEDist"]->GetMean();
    const Float_t truth_E_StdDev    = histContainer_["truthEDist"]->GetStdDev();
    const Float_t E_Res             = det_E_StdDev / det_E_Mean;
    
    // Control
    std::cout << "Truth mean, std dev: " << truth_E_Mean << ", " << truth_E_StdDev << std::endl;
    std::cout << "Measured mean, std dev: " << det_E_Mean << ", " << det_E_StdDev << std::endl;
    std::cout << "ENERGY RESOLUTION: " << E_Res * 100. << " %" << std::endl;
    */

    // Job ends
    std::cout << "The job, EnergyResolution, has ended. Thank you for your patience." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (EnergyResolution);
