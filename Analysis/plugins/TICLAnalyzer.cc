#include "TICLAnalyzer.h"

#include <iostream>
#include <array>
#include <algorithm>

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
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
// #include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

// TODO: since it is computationally expensive to iterate over calotruths multiple times,
// do it once and pass on the resulting vector into other functions

// ROOT headers

using namespace edm;

TICLAnalyzer::TICLAnalyzer ( const edm::ParameterSet& iConfig ) :
    histContainer_ (),
    
    tag_CaloParticle_MergedCaloTruth_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_MergedCaloTruth", edm::InputTag ("mix", "MergedCaloTruth", "HLT") ) ),
    
    tag_Tracks_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_Tracks", edm::InputTag("generalTracks") ) ),

    // tag_SimHits_HFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_g4SimHitsHFNoseHits", edm::InputTag("g4SimHits","HFNoseHits") ) ),
    tag_RecHits_HFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_HGCHFNoseRecHits", edm::InputTag("HGCalRecHit:HGCHFNoseRecHits") ) ),
    //tag_RecHits_EE_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_HGCEERecHits", edm::InputTag("HGCalRecHit::HGCEERecHits") ) ),
    
    tag_LayerClusters_HFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_LayerClustersHFNose", edm::InputTag("hgcalLayerClusters") ) ),
    //tag_LayerClusters_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_LayerClustersHFNose", edm::InputTag("hgcalLayerClusters") ) ),
    
    tag_Trackster_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_Trackster", edm::InputTag ("ticlTrackstersHFNoseEM") ) ),
    
    // Custom parameters
    select_PID_ ( iConfig.getUntrackedParameter<int> ("select_PID", 22 ) ),
    select_EtaLow_ ( iConfig.getUntrackedParameter<double> ("select_EtaLow", 3.4 ) ),
    select_EtaHigh_ ( iConfig.getUntrackedParameter<double> ("select_EtaHigh", 3.6 ) ),
    truth_matching_deltaR_ ( iConfig.getUntrackedParameter<double> ("truth_matching_deltaR_", 0.5 ) ),
    trackster_itername_ ( iConfig.getUntrackedParameter<std::string> ("trackster_itername", "EMn") )
    
{
    // consumes: frequent request of additional data | mayConsume: infrequent
    token_CaloParticle_MergedCaloTruth_ = mayConsume<std::vector<CaloParticle>> ( tag_CaloParticle_MergedCaloTruth_ );
//    token_SimHits_HFNose_ = mayConsume<edm::PCaloHitContainer> ( tag_SimHits_HFNose_ );

    token_Tracks_ = consumes<std::vector<reco::Track>>( tag_Tracks_);

    token_RecHits_HFNose_ = consumes<HGCRecHitCollection> ( tag_RecHits_HFNose_ );
    //token_RecHits_EE_ = consumes<HGCRecHitCollection> ( tag_RecHits_EE_ );
    
    token_LayerClusters_HFNose_ = consumes<std::vector<reco::CaloCluster>>( tag_LayerClusters_HFNose_ );
    //token_LayerClusters_ = consumes<std::vector<reco::CaloCluster>>( tag_LayerClusters_ );    
    
    token_Trackster_ = consumes<std::vector<ticl::Trackster>> ( tag_Trackster_ );
}


TICLAnalyzer::~TICLAnalyzer ()
{
    // Deconstructor
}


void TICLAnalyzer::analyze ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
    // Get Geometry
    edm::ESHandle<HGCalGeometry> handle_HGCalGeometry_HFNose ;
    iSetup.get<IdealGeometryRecord>().get( "HGCalHFNoseSensitive", handle_HGCalGeometry_HFNose );
    
    //edm::ESHandle<HGCalGeometry> handle_HGCalGeometry_EE ;
    //iSetup.get<IdealGeometryRecord>().get( "HGCalEESensitive", handle_HGCalGeometry_EE );
    
    // Get Tracks
    edm::Handle<std::vector<reco::Track>> handle_Tracks;
    iEvent.getByToken ( token_Tracks_, handle_Tracks );

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
    
    //edm::Handle<HGCRecHitCollection> handle_RecHits_EE;
    //iEvent.getByToken( token_RecHits_EE_, handle_RecHits_EE );
    
    // Get LayerClusters
    edm::Handle<std::vector<reco::CaloCluster>> handle_LayerClusters_HFNose;
    iEvent.getByToken( token_LayerClusters_HFNose_, handle_LayerClusters_HFNose );
    
    //edm::Handle<std::vector<reco::CaloCluster>> handle_LayerClusters;
    //iEvent.getByToken( token_LayerClusters_, handle_LayerClusters );
    
    // Get ticlTracksters
    edm::Handle<std::vector<ticl::Trackster>> handle_Trackster;
    iEvent.getByToken ( token_Trackster_, handle_Trackster );
    
    if ( handle_CaloParticle_MergedCaloTruth.isValid() &&
         handle_Tracks.isValid() &&
         //handle_SimHits_HFNose.isValid() &&
         handle_RecHits_HFNose.isValid() &&
         //handle_RecHits_EE.isValid() &&
         handle_LayerClusters_HFNose.isValid() &&
         //handle_LayerClusters.isValid() &&
         handle_Trackster.isValid()
       )
    {

        const std::vector<CaloParticle> & caloTruthParticles = *handle_CaloParticle_MergedCaloTruth.product();
        const HGCalGeometry* nose_geom = handle_HGCalGeometry_HFNose.product();
        //const HGCalGeometry* EE_geom = handle_HGCalGeometry_EE.product();

        fillTruthHistograms ( caloTruthParticles );
        analyzeTICLTrackster ( caloTruthParticles, *handle_Trackster.product(), trackster_itername_ );
        analyzeTrackPosition ( caloTruthParticles, *handle_Tracks.product() );
        analyzeRecHits ( caloTruthParticles, *handle_RecHits_HFNose.product(), nose_geom );
        //analyzeRecHits ( caloTruthParticles, *handle_RecHits_EE.product(), EE_geom );
        analyzeLayerClusters ( caloTruthParticles, *handle_LayerClusters_HFNose.product() );
        //analyzeLayerClusters ( caloTruthParticles, *handle_LayerClusters.product() );
        
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


void TICLAnalyzer::fillTruthHistograms ( const std::vector<CaloParticle> & caloTruthParticles )
{

    const std::vector<math::XYZTLorentzVectorF> selected_calotruths = getTruthP4 ( caloTruthParticles );
    
    for ( auto const& truth: selected_calotruths )
    {
        // Fill truth histograms
        histContainer_["truthE"]->Fill( truth.energy() );
        histContainer_["truthEta"]->Fill( truth.eta() );   
    }

}


void TICLAnalyzer::analyzeTICLTrackster ( const std::vector<CaloParticle> & caloTruthParticles, const std::vector<ticl::Trackster> & tracksters, std::string tag )
{

    const std::vector<math::XYZTLorentzVectorF> selected_calotruths = getTruthP4 ( caloTruthParticles );

    for ( auto const& trs: tracksters )
    {
        Float_t trackster_raw_energy = trs.raw_energy();
        Float_t trackster_eta = trs.barycenter().eta();
        Float_t trackster_phi = trs.barycenter().phi();
        //std::array<Float_t, 8> trackster_id_probs = trs.id_probabilities();
        // Float_t trackster_sigmas = trs.sigmas();
        // Float_t trackster_sigmasPCA = trs.sigmasPCA();
                
        for ( auto const& truth: selected_calotruths )
        {        
            Float_t dR = reco::deltaR ( trackster_eta, trackster_phi, truth.eta(), truth.phi() );

            if ( dR < (truth_matching_deltaR_))
            {
            
                if ( tag == "EMn" )
                {
                    histContainer_["RawEDist_tracksterHFNoseEM"]->Fill( trackster_raw_energy );
                    histContainer_["RawEScale_tracksterHFNoseEM"]->Fill( truth.E() - trackster_raw_energy );
                    histContainer_["DeltaR_tracksterHFNoseEM"]->Fill( dR );
                }
                else if ( tag == "TrkEMn")
                {
                    histContainer_["RawEDist_tracksterHFNoseTrkEM"]->Fill( trackster_raw_energy );
                    histContainer_["RawEScale_tracksterHFNoseTrkEM"]->Fill( truth.E() - trackster_raw_energy );
                    histContainer_["DeltaR_tracksterHFNoseTrkEM"]->Fill( dR );
                 }   
                 else if ( tag == "EM" )
                 {
                    histContainer_["RawEDist_tracksterEM"]->Fill( trackster_raw_energy );
                    histContainer_["RawEScale_tracksterEM"]->Fill( truth.E() - trackster_raw_energy );
                    histContainer_["DeltaR_tracksterEM"]->Fill( dR );
                 }
                    
            }
        }
    }
}


//void TICLAnalyzer::analyzeSimHits ( const std::vector<CaloParticle> & caloTruthParticles, const edm::PCaloHitContainer & simHits, const HGCalGeometry* geom )
//{

//// NOT FINISHED

//}


void TICLAnalyzer::analyzeRecHits ( const std::vector<CaloParticle> & caloTruthParticles, const HGCRecHitCollection & hits, const HGCalGeometry * geom )
{

    const std::vector<math::XYZTLorentzVectorF> selected_calotruths = getTruthP4 ( caloTruthParticles );

    for ( auto const& truth: selected_calotruths )
    {
        Float_t truth_eta = truth.eta();
        Float_t truth_phi = truth.phi();
        Float_t truth_energy = truth.energy();
        
        Float_t hit_sum_energy = 0;
        Float_t hit_sum_energy_layer[8] = {0};
        
        for ( auto const& hit: hits )
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

        histContainer_["ERatio_recHits_to_truth"]->Fill ( hit_sum_energy / truth_energy );
        
        histContainer_["ERatio_recHits_layer1"]->Fill ( hit_sum_energy_layer[0] / truth_energy );
        histContainer_["ERatio_recHits_layer2"]->Fill ( hit_sum_energy_layer[1] / truth_energy );
        histContainer_["ERatio_recHits_layer3"]->Fill ( hit_sum_energy_layer[2] / truth_energy );
        histContainer_["ERatio_recHits_layer4"]->Fill ( hit_sum_energy_layer[3] / truth_energy );
        histContainer_["ERatio_recHits_layer5"]->Fill ( hit_sum_energy_layer[4] / truth_energy );
        histContainer_["ERatio_recHits_layer6"]->Fill ( hit_sum_energy_layer[5] / truth_energy );
        histContainer_["ERatio_recHits_layer7"]->Fill ( hit_sum_energy_layer[6] / truth_energy );
        histContainer_["ERatio_recHits_layer8"]->Fill ( hit_sum_energy_layer[7] / truth_energy );
        
    }
}


void TICLAnalyzer::analyzeLayerClusters ( const std::vector<CaloParticle> & caloTruthParticles, const std::vector<reco::CaloCluster> & clusters )
{

    const std::vector<math::XYZTLorentzVectorF> selected_calotruths = getTruthP4 ( caloTruthParticles );

    for ( auto const& truth: selected_calotruths )
    {
        Float_t truth_eta = truth.eta();
        Float_t truth_phi = truth.phi();
        Float_t truth_energy = truth.energy();
        
        Float_t cluster_sum_energy = 0;
        Float_t cluster_sum_energy_layer[8] = {0};
        Float_t cluster_count_layer[8] = {0};
        
        for ( auto const& cl: clusters )
        {
            Float_t dR = reco::deltaR ( cl.eta(), cl.phi(), truth_eta, truth_phi );
            
            if ( dR < truth_matching_deltaR_ ) // clusters within cone of truth
            {
                int cluster_layer = rhtools_.getLayer(cl.hitsAndFractions().at(0).first);
                
                // Sanity check
//                for ( int i = 0; int i < size(cl.hitsAndFractions()); ++i )
//                {
//                    int cluster_layer = smth ()  cl.hitsAndFractions().at(i).first // layer of
//                    if ( i == 0 ) first_layer = layer;
//                    // check if they are all the same
//                }
                
                Float_t cluster_raw_energy = cl.energy();
                cluster_sum_energy += cluster_raw_energy;

                if ( cluster_layer <= 8 )
                {
                    cluster_sum_energy_layer[cluster_layer-1] += cluster_raw_energy;
                    cluster_count_layer[cluster_layer-1] ++;
                }
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
        
        histContainer_["ERatio_layerClusters_to_truth"]->Fill ( cluster_sum_energy / truth_energy );
        
        histContainer_["ERatio_layerClusters_layer1"]->Fill ( cluster_sum_energy_layer[0] / truth_energy );
        histContainer_["ERatio_layerClusters_layer2"]->Fill ( cluster_sum_energy_layer[1] / truth_energy );
        histContainer_["ERatio_layerClusters_layer3"]->Fill ( cluster_sum_energy_layer[2] / truth_energy );
        histContainer_["ERatio_layerClusters_layer4"]->Fill ( cluster_sum_energy_layer[3] / truth_energy );
        histContainer_["ERatio_layerClusters_layer5"]->Fill ( cluster_sum_energy_layer[4] / truth_energy );
        histContainer_["ERatio_layerClusters_layer6"]->Fill ( cluster_sum_energy_layer[5] / truth_energy );
        histContainer_["ERatio_layerClusters_layer7"]->Fill ( cluster_sum_energy_layer[6] / truth_energy );
        histContainer_["ERatio_layerClusters_layer8"]->Fill ( cluster_sum_energy_layer[7] / truth_energy );
    }

}


void TICLAnalyzer::analyzeTrackPosition ( const std::vector<CaloParticle> & caloTruthParticles, const std::vector<reco::Track> & generalTracks )
{

    const std::vector<math::XYZTLorentzVectorF> selected_calotruths = getTruthP4 ( caloTruthParticles );

    for ( auto const& truth: selected_calotruths )
    {
        Float_t truth_eta = truth.eta();
        Float_t truth_phi = truth.phi();
        
        Float_t min_R = 9999.;
        Float_t truth_track_dEta = -1.;
        Float_t truth_track_dPhi = -1.;
        reco::Track min_R_track;
        
        for ( auto const& track: generalTracks )
        {
            Float_t track_outerEta = track.outerEta();
            Float_t track_outerPhi = track.outerPhi();
            
            Float_t dR = reco::deltaR ( track_outerEta, track_outerPhi, truth_eta, truth_phi );
            
            if ( dR < truth_matching_deltaR_ && dR < min_R ) // clusters within cone of truth
            {
                min_R = dR;
                min_R_track = track;
                truth_track_dEta = abs(truth_eta - track_outerEta);
                truth_track_dPhi = abs(truth_phi - track_outerPhi);
            }
        }
        
        if ( min_R < truth_matching_deltaR_ )
        {
            histContainer_["dEta_Truth_Track"]->Fill(truth_track_dEta);
            histContainer_["dPhi_Truth_Track"]->Fill(truth_track_dPhi);
        }
    }

}


void TICLAnalyzer::beginJob ()
{
    edm::Service<TFileService> fs;
    
    // CaloTruth
    histContainer_["truthE"] = fs->make<TH1F>("truthE", "Truth Energy Distribution", 100, 0, 500);
    histContainer_["truthEta"] = fs->make<TH1F>("truthEta", "Truth Eta Distribution", 40, (select_EtaLow_ + select_EtaHigh_)/2 - 0.5, (select_EtaLow_ + select_EtaHigh_)/2 + 0.5);
    
    // Tracksters
    histContainer_["RawEDist_tracksterHFNoseEM"] = fs->make<TH1F>("EDist_tracksterHFNoseEM", "TracksterHFNoseEM Raw Energy Distribution", 100, 0, 500);
    histContainer_["RawEScale_tracksterHFNoseEM"] = fs->make<TH1F>("EScale_tracksterHFNoseEM", "TracksterHFNoseEM Raw Energy Scale", 100, 0, 500);
    histContainer_["DeltaR_tracksterHFNoseEM"] = fs->make<TH1F>("DeltaR_tracksterHFNoseEM", "TracksterHFNoseEM #Delta R", 20, 0, 0.5);
    
    histContainer_["RawEDist_tracksterHFNoseEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    histContainer_["RawEScale_tracksterHFNoseEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    histContainer_["DeltaR_tracksterHFNoseEM"]->GetXaxis()->SetTitle("|R_{trackster - caloParticle}|");
    
    histContainer_["RawEDist_tracksterHFNoseTrkEM"] = fs->make<TH1F>("EDist_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM Raw Energy Distribution", 100, 0, 500);
    histContainer_["RawEScale_tracksterHFNoseTrkEM"] = fs->make<TH1F>("EScale_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM Raw Energy Scale", 100, 0, 500);
    histContainer_["DeltaR_tracksterHFNoseTrkEM"] = fs->make<TH1F>("DeltaR_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM #Delta R", 20, 0, 0.5);
    
    histContainer_["RawEDist_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    histContainer_["RawEScale_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    histContainer_["DeltaR_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("|R_{trackster} - R_{caloParticle}|");
    
    histContainer_["RawEDist_tracksterEM"] = fs->make<TH1F>("EDist_tracksterEM", "TracksterEM Raw Energy Distribution", 500, 0, 500);
    histContainer_["RawEScale_tracksterEM"] = fs->make<TH1F>("EScale_tracksterEM", "TracksterEM Raw Energy Scale", 100, 0, 500);
    histContainer_["DeltaR_tracksterEM"] = fs->make<TH1F>("DeltaR_tracksterEM", "TracksterEM #Delta R", 20, 0, 0.5);
    
    histContainer_["RawEDist_tracksterEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    histContainer_["RawEScale_tracksterEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    histContainer_["DeltaR_tracksterEM"]->GetXaxis()->SetTitle("|R_{trackster} - R_{caloParticle}|");
    
    // Clusters
    histContainer_["EDist_layerClusters_scalar_sum"] = fs->make<TH1F>("EDist_layerClusters_scalar_sum", "Energy Distribution (clusters, scalar sum)", 100, 0, 1000);
    
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
    
    // Ratios
    histContainer_["ERatio_recHits_to_truth"] = fs->make<TH1F>("ERatio_recHits_to_truth", "E_{recHits} / E_{truth}", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_to_truth"] = fs->make<TH1F>("ERatio_layerClusters_to_truth", "E_{layerClusters} / E_{truth}", 100, 0, 1.0);
    
    histContainer_["ERatio_recHits_to_truth"]->GetXaxis()->SetTitle("E_{recHits} / E_{truth}");
    histContainer_["ERatio_layerClusters_to_truth"]->GetXaxis()->SetTitle("E_{layerClusters} / E_{truth}");
    
    histContainer_["ERatio_recHits_layer1"] = fs->make<TH1F>("ERatio_recHits_layer1", "E_{recHits} / E_{truth} Layer 1", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer2"] = fs->make<TH1F>("ERatio_recHits_layer2", "E_{recHits} / E_{truth} Layer 2", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer3"] = fs->make<TH1F>("ERatio_recHits_layer3", "E_{recHits} / E_{truth} Layer 3", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer4"] = fs->make<TH1F>("ERatio_recHits_layer4", "E_{recHits} / E_{truth} Layer 4", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer5"] = fs->make<TH1F>("ERatio_recHits_layer5", "E_{recHits} / E_{truth} Layer 5", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer6"] = fs->make<TH1F>("ERatio_recHits_layer6", "E_{recHits} / E_{truth} Layer 6", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer7"] = fs->make<TH1F>("ERatio_recHits_layer7", "E_{recHits} / E_{truth} Layer 7", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer8"] = fs->make<TH1F>("ERatio_recHits_layer8", "E_{recHits} / E_{truth} Layer 8", 100, 0, 1.0);
    
    histContainer_["ERatio_layerClusters_layer1"] = fs->make<TH1F>("ERatio_layerClusters_layer1", "E_{layerClusters} / E_{truth} Layer 1", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer2"] = fs->make<TH1F>("ERatio_layerClusters_layer2", "E_{layerClusters} / E_{truth} Layer 2", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer3"] = fs->make<TH1F>("ERatio_layerClusters_layer3", "E_{layerClusters} / E_{truth} Layer 3", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer4"] = fs->make<TH1F>("ERatio_layerClusters_layer4", "E_{layerClusters} / E_{truth} Layer 4", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer5"] = fs->make<TH1F>("ERatio_layerClusters_layer5", "E_{layerClusters} / E_{truth} Layer 5", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer6"] = fs->make<TH1F>("ERatio_layerClusters_layer6", "E_{layerClusters} / E_{truth} Layer 6", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer7"] = fs->make<TH1F>("ERatio_layerClusters_layer7", "E_{layerClusters} / E_{truth} Layer 7", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer8"] = fs->make<TH1F>("ERatio_layerClusters_layer8", "E_{layerClusters} / E_{truth} Layer 8", 100, 0, 1.0);
    
    // SimHits
    
    // Tracks
    histContainer_["dEta_Truth_Track"] = fs->make<TH1F>("dEta_Truth_Track", "#eta difference between track and calotruth at HGCal interface", 50, 0, 0.5);
    histContainer_["dPhi_Truth_Track"] = fs->make<TH1F>("dPhi_Truth_Track", "#phi difference between track and calotruth at HGCal interface", 50, 0, 0.5);
    
    histContainer_["dEta_Truth_Track"]->GetXaxis()->SetTitle("|#eta_{caloTruth} - #eta_{generalTrack}|");
    histContainer_["dPhi_Truth_Track"]->GetXaxis()->SetTitle("|#phi_{caloTruth} - #phi_{generalTrack}|");
}


void TICLAnalyzer::endJob ()
{
    // Job ends
    std::cout << "The job, TICLAnalyzer, has ended. Thank you for your patience." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (TICLAnalyzer);
