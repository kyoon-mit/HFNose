#include "TICLAnalyzer.h"

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

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
// #include "DataFormats/HGCalReco/interface/TICLCandidate.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

// ROOT headers

using namespace edm;

TICLAnalyzer::TICLAnalyzer ( const edm::ParameterSet& iConfig ) :
    histContainer_ (),
    
    tag_CaloParticle_MergedCaloTruth_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_MergedCaloTruth", edm::InputTag ("mix", "MergedCaloTruth") ) ),
    
//    tag_SimHits_HFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_g4SimHitsHFNoseHits", edm::InputTag("g4SimHits","HFNoseHits") ) ),
    tag_RecHits_HFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_HGCHFNoseRecHits", edm::InputTag("HGCalRecHit:HGCHFNoseRecHits") ) ),
    
    tag_LayerClusters_HFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_LayerClustersHFNose", edm::InputTag("hgcalLayerClustersHFNose") ) ),
    
    tag_Trackster_HFNoseEM_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_TracksterHFNoseEM", edm::InputTag ("ticlTrackstersHFNoseEM") ) ),
    tag_Trackster_HFNoseTrkEM_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_TracksterHFNoseTrk", edm::InputTag ("ticlTrackstersHFNoseTrkEM") ) ),
    
    // Custom parameters
    select_PID_ ( 22 ),
    select_EtaLow_ ( 3.49 ),
    select_EtaHigh_ ( 3.51 ),
    truth_matching_deltaR_ ( 0.2 )
{
    // consumes: frequent request of additional data | mayConsume: infrequent
    token_CaloParticle_MergedCaloTruth_ = mayConsume<std::vector<CaloParticle>> ( tag_CaloParticle_MergedCaloTruth_ );
    
//    token_SimHits_HFNose_ = mayConsume<edm::PCaloHitContainer> ( tag_SimHits_HFNose_ );
    token_RecHits_HFNose_ = consumes<HGCRecHitCollection> ( tag_RecHits_HFNose_ );
    
    token_LayerClusters_HFNose_ = consumes<std::vector<reco::CaloCluster>>( tag_LayerClusters_HFNose_ );

    token_Trackster_HFNoseEM_ = consumes<std::vector<ticl::Trackster>> ( tag_Trackster_HFNoseEM_ );
    token_Trackster_HFNoseTrkEM_ = consumes<std::vector<ticl::Trackster>> ( tag_Trackster_HFNoseTrkEM_ );
}


TICLAnalyzer::~TICLAnalyzer ()
{
    // Deconstructor
}


void TICLAnalyzer::analyze ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
    // Get Geometry
    edm::ESHandle<HGCalGeometry> handle_HGCalGeometry;
    iSetup.get<IdealGeometryRecord>().get( "HGCalHFNoseSensitive", handle_HGCalGeometry );

    // Get CaloTruth
    edm::Handle<std::vector<CaloParticle>> handle_CaloParticle_MergedCaloTruth;
    iEvent.getByToken ( token_CaloParticle_MergedCaloTruth_, handle_CaloParticle_MergedCaloTruth );
    
    // Get SimHits
//    edm::Handle<edm::PCaloHitContainer> handle_SimHits_HFNose;
//    iEvent.getByToken( token_SimHits_HFNose_, handle_SimHits_HFNose );
//    
    // Get RecHits
    edm::Handle<HGCRecHitCollection> handle_RecHits_HFNose;
    iEvent.getByToken( token_RecHits_HFNose_, handle_RecHits_HFNose );
    
    // Get LayerClusters
    edm::Handle<std::vector<reco::CaloCluster>> handle_LayerClusters_HFNose;
    iEvent.getByToken( token_LayerClusters_HFNose_, handle_LayerClusters_HFNose );
    
    // Get ticlTrackstersHFNoseEM
    edm::Handle<std::vector<ticl::Trackster>> handle_Trackster_HFNoseEM;
    iEvent.getByToken ( token_Trackster_HFNoseEM_, handle_Trackster_HFNoseEM );
    
    // Get ticlTrackstersHFNoseTrkEM
    edm::Handle<std::vector<ticl::Trackster>> handle_Trackster_HFNoseTrkEM;
    iEvent.getByToken ( token_Trackster_HFNoseTrkEM_, handle_Trackster_HFNoseTrkEM );
    
    if ( handle_CaloParticle_MergedCaloTruth.isValid() &&
//         handle_SimHits_HFNose.isValid() &&
         handle_RecHits_HFNose.isValid() &&
         handle_LayerClusters_HFNose.isValid() &&
         handle_Trackster_HFNoseEM.isValid() //&&
//         handle_Trackster_HFNoseTrkEM.isValid()
       )
    {

        const std::vector<CaloParticle> & caloTruthParticles = *handle_CaloParticle_MergedCaloTruth.product();
        const HGCalGeometry* geom = handle_HGCalGeometry.product();

        analyzeTICLTrackster ( caloTruthParticles, *handle_Trackster_HFNoseEM.product(), "EMn" );
//        analyzeTICLTrackster ( caloTruthParticles, *handle_Trackster_HFNoseTrkEM.product(), "TrkEMn" );
        analyzeRecHits ( caloTruthParticles, *handle_RecHits_HFNose.product(), geom );
        analyzeLayerClusters ( caloTruthParticles, *handle_LayerClusters_HFNose.product() );
        
    }
    else std::cout << "Handle(s) invalid!" << std::endl;
}


std::vector<math::XYZTLorentzVectorF> TICLAnalyzer::getTruthP4 ( const std::vector<CaloParticle> & caloTruthParticles )
{
    std::vector<math::XYZTLorentzVectorF> container;
    
    for ( auto const& ct: caloTruthParticles )
    {
        if ( abs(ct.pdgId()) == select_PID_
                && abs(ct.eta()) > select_EtaLow_
                && abs(ct.eta()) < select_EtaHigh_ )
        {
            container.push_back ( (math::XYZTLorentzVectorF) ct.p4() );
        }
    }
    
    return container;
}


void TICLAnalyzer::analyzeTICLTrackster ( const std::vector<CaloParticle> & caloTruthParticles, const std::vector<ticl::Trackster> & tracksters, std::string tag )
{

    const std::vector<math::XYZTLorentzVectorF> selected_calotruths = getTruthP4 ( caloTruthParticles );

    for ( auto const& trs: tracksters )
    {
        Float_t trackster_raw_energy = trs.raw_energy();
        Float_t trackster_eta = trs.barycenter().eta();
        Float_t trackster_phi = trs.barycenter().phi();
        
        for ( auto const& truth: selected_calotruths )
        {
            histContainer_["truthE"]->Fill( truth.E() );
            histContainer_["truthEta"]->Fill( truth.Eta() );
            
            Float_t dR = reco::deltaR ( trackster_eta, trackster_phi, truth.eta(), truth.phi() );

            if ( dR < truth_matching_deltaR_ )
            {
                /* Add Trackster dissecting code here */
                
                // 1. Number of LayerClusters in Trackster, cumulative and per layer
                int count_total_layerClusters = 0;
                int count_perLayer_layerClusters[8] = 0;
                
                // 2. Cosine of angle between doublets (related to the min_cos_theta cut)
                
                // 3. Cosine of pointing angle of outer doublet (related to min_cos_pointing cut)
                
                // 4. Number of missing layers
                
                // 5. Timing
            
            
            
                if ( tag == "EMn" )
                {
                    histContainer_["EDist_tracksterHFNoseEM"]->Fill( trackster_raw_energy );
                    histContainer_["EScale_tracksterHFNoseEM"]->Fill( truth.E() - trackster_raw_energy );
                    histContainer_["DeltaR_tracksterHFNoseEM"]->Fill( dR );
                }
                else if ( tag == "TrkEMn")
                {
                    histContainer_["EDist_tracksterHFNoseTrkEM"]->Fill( trackster_raw_energy );
                    histContainer_["EScale_tracksterHFNoseTrkEM"]->Fill( truth.E() - trackster_raw_energy );
                    histContainer_["DeltaR_tracksterHFNoseTrkEM"]->Fill( dR );
                 }   
                    // case "EM":
                    
                    // case "TrkEM":                
            }
        }
    }
}


//void TICLAnalyzer::analyzeSimHits ( const std::vector<CaloParticle> & caloTruthParticles, const edm::PCaloHitContainer & simHits, const HGCalGeometry* geom )
//{

//// NOT FINISHED

//    const std::vector<math::XYZTLorentzVectorF> selected_calotruths = getTruthP4 ( caloTruthParticles );

//    for ( auto const& truth: selected_calotruths )
//    {
//    
//        Float_t truth_eta = truth.eta();
//        Float_t truth_phi = truth.phi();
//        
//        Float_t hit_sum_energy = 0;
//        Float_t hit_sum_energy_layer[8] = {0};
//    
//        for ( auto const& hit: simHits )
//        {
//            uint32_t hit_id = hit.id();
//            HFNoseDetId detId = HFNoseDetId(hit_id);
//            int hit_layer = detId.layer();
//            
//            const GlobalPoint & hit_globalPosition = geom->getPosition(hit_id);
//            Float_t hit_eta = hit_globalPosition.eta();
//            Float_t hit_phi = hit_globalPosition.phi();
//    
//            Float_t dR = reco::deltaR ( hit_eta, hit_phi, truth_eta, truth_phi );
//            // https://github.com/mariadalfonso/HGCnoseUtils/blob/master/Analyzer/plugins/GenAnalyzer.cc#L515-L623
//            
//            if ( dR < truth_matching_deltaR_ )
//            {
//                Float_t hit_raw_energy = hit.energy();
//                
//            }
//        }
//    }
//}


void TICLAnalyzer::analyzeRecHits ( const std::vector<CaloParticle> & caloTruthParticles, const HGCRecHitCollection & recHits, const HGCalGeometry * geom )
{

    const std::vector<math::XYZTLorentzVectorF> selected_calotruths = getTruthP4 ( caloTruthParticles );

    for ( auto const& truth: selected_calotruths )
    {
        Float_t truth_eta = truth.eta();
        Float_t truth_phi = truth.phi();
        
        Float_t hit_sum_energy = 0;
        Float_t hit_sum_energy_layer[8] = {0};
        
        for ( auto const& hit: recHits )
        {
            uint32_t hit_id = hit.id();
        
            const GlobalPoint & hit_globalPosition = geom->getPosition(hit_id);
            Float_t hit_eta = hit_globalPosition.eta();
            Float_t hit_phi = hit_globalPosition.phi();

            Float_t dR = reco::deltaR ( hit_eta, hit_phi, truth_eta, truth_phi );
            
            if ( dR < truth_matching_deltaR_ )
            {
                int hit_layer = HFNoseDetId(hit_id).layer();
                
                Float_t hit_raw_energy = hit.energy();

                hit_sum_energy += hit_raw_energy;
                hit_sum_energy_layer[hit_layer-1] += hit_raw_energy;
            }
        }
        
        histContainer_["EDist_recHits"]->Fill ( hit_sum_energy );
        
        histContainer_["EDist_recHits_layer1"]->Fill ( hit_sum_energy_layer[0] );
        histContainer_["EDist_recHits_layer2"]->Fill ( hit_sum_energy_layer[1] );
        histContainer_["EDist_recHits_layer3"]->Fill ( hit_sum_energy_layer[2] );
        histContainer_["EDist_recHits_layer4"]->Fill ( hit_sum_energy_layer[3] );
        histContainer_["EDist_recHits_layer5"]->Fill ( hit_sum_energy_layer[4] );
        histContainer_["EDist_recHits_layer6"]->Fill ( hit_sum_energy_layer[5] );
        histContainer_["EDist_recHits_layer7"]->Fill ( hit_sum_energy_layer[6] );
        histContainer_["EDist_recHits_layer8"]->Fill ( hit_sum_energy_layer[7] );
        
    }
}


void TICLAnalyzer::analyzeLayerClusters ( const std::vector<CaloParticle> & caloTruthParticles, const std::vector<reco::CaloCluster> & clusters )
{

    const std::vector<math::XYZTLorentzVectorF> selected_calotruths = getTruthP4 ( caloTruthParticles );

    for ( auto const& truth: selected_calotruths )
    {
        Float_t truth_eta = truth.eta();
        Float_t truth_phi = truth.phi();
        
        Float_t cluster_sum_energy = 0;
        Float_t cluster_sum_energy_layer[8] = {0};
        Float_t cluster_count_layer[8] = {0};
        
        for ( auto const& cl: clusters )
        {
            Float_t dR = reco::deltaR ( cl.eta(), cl.phi(), truth_eta, truth_phi );
            
            if ( dR < truth_matching_deltaR_ ) // clusters within cone of truth
            {
                int cluster_layer = HFNoseDetId(cl.hitsAndFractions().at(0).first).layer();
                
                Float_t cluster_raw_energy = cl.energy();
                
                cluster_sum_energy += cluster_raw_energy;
                cluster_sum_energy_layer[cluster_layer-1] += cluster_raw_energy;
                cluster_count_layer[cluster_layer-1] ++;
            }
        }
        
        histContainer_["EDist_layerClusters_scalar_sum"]->Fill ( cluster_sum_energy );
        
        histContainer_["EDist_layerClusters_layer1"]->Fill( cluster_sum_energy_layer[0] );
        histContainer_["EDist_layerClusters_layer2"]->Fill( cluster_sum_energy_layer[1] );
        histContainer_["EDist_layerClusters_layer3"]->Fill( cluster_sum_energy_layer[2] );
        histContainer_["EDist_layerClusters_layer4"]->Fill( cluster_sum_energy_layer[3] );
        histContainer_["EDist_layerClusters_layer5"]->Fill( cluster_sum_energy_layer[4] );
        histContainer_["EDist_layerClusters_layer6"]->Fill( cluster_sum_energy_layer[5] );
        histContainer_["EDist_layerClusters_layer7"]->Fill( cluster_sum_energy_layer[6] );
        histContainer_["EDist_layerClusters_layer8"]->Fill( cluster_sum_energy_layer[7] );
        
        histContainer_["Count_layerClusters_layer1"]->Fill( cluster_count_layer[0] );
        histContainer_["Count_layerClusters_layer2"]->Fill( cluster_count_layer[1] );
        histContainer_["Count_layerClusters_layer3"]->Fill( cluster_count_layer[2] );
        histContainer_["Count_layerClusters_layer4"]->Fill( cluster_count_layer[3] );
        histContainer_["Count_layerClusters_layer5"]->Fill( cluster_count_layer[4] );
        histContainer_["Count_layerClusters_layer6"]->Fill( cluster_count_layer[5] );
        histContainer_["Count_layerClusters_layer7"]->Fill( cluster_count_layer[6] );
        histContainer_["Count_layerClusters_layer8"]->Fill( cluster_count_layer[7] );
    }

}


void TICLAnalyzer::beginJob ()
{
    edm::Service<TFileService> fs;
    
    // CaloTruth
    histContainer_["truthE"] = fs->make<TH1F>("truthE", "Truth Energy Distribution", 100, 0, 1000);
    histContainer_["truthEta"] = fs->make<TH1F>("truthEta", "Truth Eta Distribution", 40, 3.0, 4.0);
    
    // Tracksters
    histContainer_["EDist_tracksterHFNoseEM"] = fs->make<TH1F>("EDist_tracksterHFNoseEM", "TracksterHFNoseEM  Energy Distribution", 100, 0, 1000);
    histContainer_["EScale_tracksterHFNoseEM"] = fs->make<TH1F>("EScale_tracksterHFNoseEM", "TracksterHFNoseEM Energy Scale", 100, 0, 1000);
    histContainer_["DeltaR_tracksterHFNoseEM"] = fs->make<TH1F>("DeltaR_tracksterHFNoseEM", "TracksterHFNoseEM #Delta R", 20, 0, 0.5);
    
    histContainer_["EDist_tracksterHFNoseEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    histContainer_["EScale_tracksterHFNoseEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    histContainer_["DeltaR_tracksterHFNoseEM"]->GetXaxis()->SetTitle("|#Delta R_{trackster - caloParticle}|");
    
    histContainer_["EDist_tracksterHFNoseTrkEM"] = fs->make<TH1F>("EDist_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM  Energy Distribution", 100, 0, 1000);
    histContainer_["EScale_tracksterHFNoseTrkEM"] = fs->make<TH1F>("EScale_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM Energy Scale", 100, 0, 1000);
    histContainer_["DeltaR_tracksterHFNoseTrkEM"] = fs->make<TH1F>("DeltaR_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM #Delta R", 20, 0, 0.5);
    
    histContainer_["EDist_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    histContainer_["EScale_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    histContainer_["DeltaR_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("|#Delta R_{trackster - caloParticle}|");
    
    // Clusters
    histContainer_["EDist_layerClusters_scalar_sum"] = fs->make<TH1F>("EDist_layerClusters_scalar_sum", "Energy Distribution (clusters, scaler sum)", 100, 0, 1000);
    
    histContainer_["EDist_layerClusters_layer1"] = fs->make<TH1F>("EDist_layerClusters_layer1", "Energy Distribution Layer 1 (clusters)", 100, 0, 1000);
    histContainer_["EDist_layerClusters_layer2"] = fs->make<TH1F>("EDist_layerClusters_layer2", "Energy Distribution Layer 2 (clusters)", 100, 0, 1000);
    histContainer_["EDist_layerClusters_layer3"] = fs->make<TH1F>("EDist_layerClusters_layer3", "Energy Distribution Layer 3 (clusters)", 100, 0, 1000);
    histContainer_["EDist_layerClusters_layer4"] = fs->make<TH1F>("EDist_layerClusters_layer4", "Energy Distribution Layer 4 (clusters)", 100, 0, 1000);
    histContainer_["EDist_layerClusters_layer5"] = fs->make<TH1F>("EDist_layerClusters_layer5", "Energy Distribution Layer 5 (clusters)", 100, 0, 1000);
    histContainer_["EDist_layerClusters_layer6"] = fs->make<TH1F>("EDist_layerClusters_layer6", "Energy Distribution Layer 6 (clusters)", 100, 0, 1000);
    histContainer_["EDist_layerClusters_layer7"] = fs->make<TH1F>("EDist_layerClusters_layer7", "Energy Distribution Layer 7 (clusters)", 100, 0, 1000);
    histContainer_["EDist_layerClusters_layer8"] = fs->make<TH1F>("EDist_layerClusters_layer8", "Energy Distribution Layer 8 (clusters)", 100, 0, 1000);
    
    histContainer_["Count_layerClusters_layer1"] = fs->make<TH1F>("Count_layerClusters_layer1", "Number of Clusters Layer 1", 50, 0, 50);
    histContainer_["Count_layerClusters_layer2"] = fs->make<TH1F>("Count_layerClusters_layer2", "Number of Clusters Layer 2", 50, 0, 50);
    histContainer_["Count_layerClusters_layer3"] = fs->make<TH1F>("Count_layerClusters_layer3", "Number of Clusters Layer 3", 50, 0, 50);
    histContainer_["Count_layerClusters_layer4"] = fs->make<TH1F>("Count_layerClusters_layer4", "Number of Clusters Layer 4", 50, 0, 50);
    histContainer_["Count_layerClusters_layer5"] = fs->make<TH1F>("Count_layerClusters_layer5", "Number of Clusters Layer 5", 50, 0, 50);
    histContainer_["Count_layerClusters_layer6"] = fs->make<TH1F>("Count_layerClusters_layer6", "Number of Clusters Layer 6", 50, 0, 50);
    histContainer_["Count_layerClusters_layer7"] = fs->make<TH1F>("Count_layerClusters_layer7", "Number of Clusters Layer 7", 50, 0, 50);
    histContainer_["Count_layerClusters_layer8"] = fs->make<TH1F>("Count_layerClusters_layer8", "Number of Clusters Layer 8", 50, 0, 50);
    
    // RecHits
    histContainer_["EDist_recHits"] = fs->make<TH1F>("EDist_recHits", "Energy Distribution (recHits)", 100, 0, 1000);
    
    histContainer_["EDist_recHits_layer1"] = fs->make<TH1F>("EDist_recHits_layer1", "Energy Distribution Layer 1 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer2"] = fs->make<TH1F>("EDist_recHits_layer2", "Energy Distribution Layer 2 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer3"] = fs->make<TH1F>("EDist_recHits_layer3", "Energy Distribution Layer 3 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer4"] = fs->make<TH1F>("EDist_recHits_layer4", "Energy Distribution Layer 4 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer5"] = fs->make<TH1F>("EDist_recHits_layer5", "Energy Distribution Layer 5 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer6"] = fs->make<TH1F>("EDist_recHits_layer6", "Energy Distribution Layer 6 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer7"] = fs->make<TH1F>("EDist_recHits_layer7", "Energy Distribution Layer 7 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer8"] = fs->make<TH1F>("EDist_recHits_layer8", "Energy Distribution Layer 8 (recHits)", 100, 0, 1000);
    
    // SimHits
}


void TICLAnalyzer::endJob ()
{
    // Job ends
    std::cout << "The job, TICLAnalyzer, has ended. Thank you for your patience." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (TICLAnalyzer);
