#include "TICLAnalyzer.h"

#include <iostream>

#include "TH1.h"
#include "TH2.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

using namespace edm;

TICLAnalyzer::TICLAnalyzer ( const edm::ParameterSet& iConfig ) :
    TH1Container_ (),
    TH2Container_ (),
    tag_CaloParticle_MergedCaloTruth_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_MergedCaloTruth", edm::InputTag ("mix", "MergedCaloTruth", "HLT") ) ),
    tag_Tracks_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_Tracks", edm::InputTag("generalTracks") ) ),
    tag_RecHits_HFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_HGCHFNoseRecHits", edm::InputTag("HGCalRecHit:HGCHFNoseRecHits") ) ),
    tag_LayerClusters_HFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_LayerClustersHFNose", edm::InputTag("hgcalLayerClustersHFNose") ) ),
    tag_Trackster_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_Trackster", edm::InputTag ("ticlTrackstersHFNoseEM") ) ),
    
    // Custom parameters
    select_PID_ ( iConfig.getUntrackedParameter<int> ("select_PID", 22) ),
    select_EtaLow_ ( iConfig.getUntrackedParameter<double> ("select_EtaLow", 3.48) ),
    select_EtaHigh_ ( iConfig.getUntrackedParameter<double> ("select_EtaHigh", 3.52) ),
    truth_matching_deltaR_ ( iConfig.getUntrackedParameter<double> ("truth_matching_deltaR_", 0.5) ),
    trackster_itername_ ( iConfig.getUntrackedParameter<std::string> ("trackster_itername", "EMn") ),
    cutTk_ ( iConfig.getUntrackedParameter<std::string> ("cutTk"))
    
{
    token_CaloParticle_MergedCaloTruth_ = mayConsume<std::vector<CaloParticle>> ( tag_CaloParticle_MergedCaloTruth_ );
    token_Tracks_ = consumes<std::vector<reco::Track>>( tag_Tracks_);
    token_RecHits_HFNose_ = consumes<HGCRecHitCollection> ( tag_RecHits_HFNose_ );
    token_LayerClusters_HFNose_ = consumes<std::vector<reco::CaloCluster>>( tag_LayerClusters_HFNose_ );
    
    token_TracksterHFNoseEM_ = consumes<std::vector<ticl::Trackster>> ( edm::InputTag ("ticlTrackstersHFNoseEM") );
    token_TracksterHFNoseTrkEM_ = consumes<std::vector<ticl::Trackster>> ( edm::InputTag ("ticlTrackstersHFNoseTrkEM") );
    token_TracksterHFNoseHAD_ = consumes<std::vector<ticl::Trackster>> ( edm::InputTag ("ticlTrackstersHFNoseHAD") );
    token_TracksterHFNoseMIP_ = consumes<std::vector<ticl::Trackster>> ( edm::InputTag ("ticlTrackstersHFNoseMIP") );
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
    
    // Get BField
    edm::ESHandle<MagneticField> handle_MagneticField;
    iSetup.get<IdealMagneticFieldRecord>().get( handle_MagneticField );
    
    // Get HGCalDDDConstants
    edm::ESHandle<HGCalDDDConstants> handle_HGCalDDDConstants;
    iSetup.get<IdealGeometryRecord>().get( "HGCalHFNoseSensitive", handle_HGCalDDDConstants );
    
    // Get Propagator
    edm::ESHandle<Propagator> handle_Propagator1;
    iSetup.get<TrackingComponentsRecord>().get( "RungeKuttaTrackerPropagator", handle_Propagator1 );
    
    edm::ESHandle<Propagator> handle_Propagator2;
    iSetup.get<TrackingComponentsRecord>().get( "SteppingHelixPropagatorAny", handle_Propagator2 );
    
    // Get Tracks
    edm::Handle<std::vector<reco::Track>> handle_Tracks;
    iEvent.getByToken ( token_Tracks_, handle_Tracks );

    // Get CaloTruth
    edm::Handle<std::vector<CaloParticle>> handle_CaloParticle_MergedCaloTruth;
    iEvent.getByToken ( token_CaloParticle_MergedCaloTruth_, handle_CaloParticle_MergedCaloTruth );
    
    // Get RecHits
    edm::Handle<HGCRecHitCollection> handle_RecHits_HFNose;
    iEvent.getByToken( token_RecHits_HFNose_, handle_RecHits_HFNose );
    
    // Get LayerClusters
    edm::Handle<std::vector<reco::CaloCluster>> handle_LayerClusters_HFNose;
    iEvent.getByToken( token_LayerClusters_HFNose_, handle_LayerClusters_HFNose );
    
    // Get ticlTracksters
    edm::Handle<std::vector<ticl::Trackster>> handle_TracksterHFNoseEM;
    edm::Handle<std::vector<ticl::Trackster>> handle_TracksterHFNoseTrkEM;
    edm::Handle<std::vector<ticl::Trackster>> handle_TracksterHFNoseHAD;
    edm::Handle<std::vector<ticl::Trackster>> handle_TracksterHFNoseMIP;
    iEvent.getByToken ( token_TracksterHFNoseEM_, handle_TracksterHFNoseEM );
    iEvent.getByToken ( token_TracksterHFNoseTrkEM_, handle_TracksterHFNoseTrkEM );
    iEvent.getByToken ( token_TracksterHFNoseHAD_, handle_TracksterHFNoseHAD );
    iEvent.getByToken ( token_TracksterHFNoseMIP_, handle_TracksterHFNoseMIP );

    if ( handle_CaloParticle_MergedCaloTruth.isValid() &&
         handle_Tracks.isValid() &&
         handle_RecHits_HFNose.isValid() &&
         handle_LayerClusters_HFNose.isValid() &&
         handle_TracksterHFNoseEM.isValid() &&
         handle_TracksterHFNoseTrkEM.isValid() &&
         handle_TracksterHFNoseHAD.isValid() &&
         handle_TracksterHFNoseMIP.isValid()
       )
    {
        const HGCalGeometry* nose_geom = handle_HGCalGeometry_HFNose.product();
        const MagneticField* b_field = handle_MagneticField.product();
        const HGCalDDDConstants* hgcons = handle_HGCalDDDConstants.product();
        
        const std::vector<math::XYZTLorentzVectorF> selected_calotruths = getTruthP4 ( *handle_CaloParticle_MergedCaloTruth.product() );
        const std::vector<std::pair<GlobalPoint, Int_t>> simhits_globalpoints = getSimHitGlobalPoint ( *handle_CaloParticle_MergedCaloTruth.product(), nose_geom );
        fillTruthHistograms ( selected_calotruths, simhits_globalpoints );
        
        buildFirstLayers ( hgcons );
        analyzeTrackPosition ( simhits_globalpoints, *handle_Tracks.product(), b_field, *handle_Propagator1 );
        propagatorSpeedTest ( *handle_Tracks.product(), b_field, *handle_Propagator1, *handle_Propagator2 );
        analyzeRecHits ( selected_calotruths, *handle_RecHits_HFNose.product(), nose_geom );
        analyzeLayerClusters ( selected_calotruths, *handle_LayerClusters_HFNose.product() );
        analyzeTICLTrackster ( selected_calotruths, *handle_TracksterHFNoseEM.product(), "EMn" );
        analyzeTICLTrackster ( selected_calotruths, *handle_TracksterHFNoseTrkEM.product(), "TrkEMn" );
        analyzeTICLTrackster ( selected_calotruths, *handle_TracksterHFNoseHAD.product(), "HADn" );
        analyzeTICLTrackster ( selected_calotruths, *handle_TracksterHFNoseMIP.product(), "MIPn" );
    }
    else std::cout << "Handle(s) invalid! Please investigate." << std::endl;
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


std::vector<std::pair<GlobalPoint, Int_t>> TICLAnalyzer::getSimHitGlobalPoint ( const std::vector<CaloParticle> & caloTruthParticles, const HGCalGeometry * geom )
{ // SimHit GlobalPoint, number of simclusters

    std::vector<std::pair<GlobalPoint, Int_t>> container;
    
    for ( auto const& ct: caloTruthParticles )
    {
        if ( abs(ct.pdgId()) == select_PID_
                && abs(ct.eta()) > select_EtaLow_
                && abs(ct.eta()) < select_EtaHigh_ )
        {
        
            const SimClusterRefVector & ct_sim_clusters = ct.simClusters();
            
            Int_t    num_sim_clusters = 0;
            Double_t simhit_max_energy = 0.;
            uint32_t simhit_detid = 0;
            Int_t    num_sim_clusters_first_layer = 0; // # sim clusters in layer 1
            
            for ( auto const& it_sc: ct_sim_clusters )
            {
                const SimCluster sc = *it_sc;
                //Int_t first_layer = 99;
                num_sim_clusters++;
                Int_t simcluster_seed_layer = 99; // # layer where it begins

                for ( auto const & he: sc.hits_and_energies() )
                {
                    Int_t layer = rhtools_.getLayer(he.first);
                    if ( layer < simcluster_seed_layer )
                    {
                        simcluster_seed_layer = rhtools_.getLayer(he.first);
                    }
                    if ( layer == 1 && simhit_max_energy < he.second ) 
                    {
                        simhit_detid = he.first;
                        simhit_max_energy = he.second;
                    }
                }
                
                if ( simcluster_seed_layer < 99 )
                {
                    TH1Container_["SimCluster_Seed_Layer"]->Fill( simcluster_seed_layer );
                    if ( simcluster_seed_layer == 1 ) num_sim_clusters_first_layer++;
                }
            }

            if ( num_sim_clusters > 1 ) continue;
            
            try
            {
                const GlobalPoint & hit_globalPosition = geom->getPosition(simhit_detid);
                if ( abs(hit_globalPosition.eta()) > select_EtaLow_
                        && abs(hit_globalPosition.eta()) < select_EtaHigh_ )
                {
                    container.push_back ( std::make_pair(hit_globalPosition, num_sim_clusters) );
                    TH1Container_["nSimClusters_layer1"]->Fill( num_sim_clusters_first_layer );
                }
            }
            catch ( cms::Exception const& e )
            {
                std::cout << e.what() << std::endl;
            }
        }
    }
    return container;
}


void TICLAnalyzer::fillTruthHistograms ( const std::vector<math::XYZTLorentzVectorF> & caloTruthP4, const std::vector<std::pair<GlobalPoint, Int_t>> & simHitsGlobalPoints )
{
    for ( auto const& truth: caloTruthP4 )
    {
        // Fill truth histograms
        TH1Container_["truthE"]->Fill( truth.energy() );
        TH1Container_["truthEta"]->Fill( truth.eta() );   
    }
    for ( auto const& hit_point: simHitsGlobalPoints )
    {
        TH1Container_["first_layer_SimHit_z"]->Fill( abs(hit_point.first.z()) );
        TH1Container_["first_layer_SimHit_eta"]->Fill( abs(hit_point.first.eta()) );
        TH2Container_["nSimCluster_SimHitEta"]->Fill( hit_point.second, abs(hit_point.first.eta()) );
    }
}


void TICLAnalyzer::buildFirstLayers ( const HGCalDDDConstants* hgcons )
// Copied from RecoHGCal/TICL/plugins/SeedingRegionByTracks.cc
{
    float zVal = hgcons->waferZ(1, true);
    std::pair<double, double> rMinMax = hgcons->rangeR(zVal, true);

    for (int iSide = 0; iSide < 2; ++iSide)
    {
        float zSide = (iSide == 0) ? (-1. * zVal) : zVal;
        firstDisk_[iSide] =
        std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
                                              Disk::RotationType(),
                                              SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5, zSide + 0.5)).get());
    }
}


void TICLAnalyzer::analyzeTrackPosition ( const std::vector<math::XYZTLorentzVectorF> & caloTruthP4, const std::vector<reco::Track> & generalTracks, const MagneticField * BField, const Propagator & prop )
{

    for ( auto const& truth: caloTruthP4 )
    {
        Float_t truth_eta = truth.eta();
        Float_t truth_phi = truth.phi();
        
        Float_t min_R = 999.;
        Float_t truth_prop_dEta = 999.;
        Float_t truth_prop_dPhi = 999.;
        
        Float_t track_prop_dEta = 999.;
        Float_t track_prop_dPhi = 999.;
        
        for ( auto const& track: generalTracks )
        {
            if ( !cutTk_(track) )
            {
                continue;
            }
            
            FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(track, BField);
            int iSide = int(track.eta() > 0);
            TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
            
            if (tsos.isValid())
            {            
                Float_t tsos_eta = tsos.globalPosition().eta();
                Float_t tsos_phi = tsos.globalPosition().phi();
                
                Float_t dR = reco::deltaR ( tsos_eta, tsos_phi, truth_eta, truth_phi );

                if ( dR < truth_matching_deltaR_ && dR < min_R ) // within cone of truth
                {
                    min_R = dR;
                    truth_prop_dEta = truth_eta - tsos_eta;
                    truth_prop_dPhi = truth_phi - tsos_phi;
                    track_prop_dEta = track.outerEta() - tsos_eta;
                    track_prop_dPhi = track.outerPhi() - tsos_phi;
                }
            }
            
        }
        
        if ( min_R < truth_matching_deltaR_ )
        {
            // simhit & propagator histograms
            TH1Container_["dEta_Truth_Prop"]->Fill(truth_prop_dEta);
            TH1Container_["dPhi_Truth_Prop"]->Fill(truth_prop_dPhi);
            TH1Container_["dR_Truth_Prop"]->Fill(min_R);
            
            // simhit & track histograms
            TH1Container_["dEta_Track_Prop"]->Fill(abs(track_prop_dEta));
            TH1Container_["dPhi_Track_Prop"]->Fill(abs(track_prop_dPhi));
        }
    }

}


void TICLAnalyzer::propagatorSpeedTest ( const std::vector<reco::Track> & generalTracks, const MagneticField * BField, const Propagator & prop1, const Propagator & prop2 )
{

    for ( auto const& track: generalTracks )
    {
        if ( !cutTk_(track) )
        {
            continue;
        }
        
        FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(track, BField);
        int iSide = int(track.eta() > 0);
        TrajectoryStateOnSurface tsos1 = prop1.propagate(fts, firstDisk_[iSide]->surface());
        TrajectoryStateOnSurface tsos2 = prop2.propagate(fts, firstDisk_[iSide]->surface());
    }

}


void TICLAnalyzer::analyzeTrackPosition ( const std::vector<std::pair<GlobalPoint, Int_t>> & simHitsGlobalPoints, const std::vector<reco::Track> & generalTracks, const MagneticField * BField, const Propagator & prop )
{

    for ( auto const& hit_point: simHitsGlobalPoints )
    {
        Float_t truth_eta = hit_point.first.eta();
        Float_t truth_phi = hit_point.first.phi();
        
        Float_t min_R = 999.;
        Float_t truth_prop_dEta = 999.;
        Float_t truth_prop_dPhi = 999.;
        
        Float_t track_prop_dEta = 999.;
        Float_t track_prop_dPhi = 999.;
        
        Int_t num_sim_clusters = 0;
        
        for ( auto const& track: generalTracks )
        {
            if ( !cutTk_(track) )
            {
                continue;
            }
            
            FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(track, BField);
            int iSide = int(track.eta() > 0);
            TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
            
            if (tsos.isValid())
            {            
                Float_t tsos_eta = tsos.globalPosition().eta();
                Float_t tsos_phi = tsos.globalPosition().phi();
                
                Float_t dR = reco::deltaR ( tsos_eta, tsos_phi, truth_eta, truth_phi );

                if ( dR < truth_matching_deltaR_ && dR < min_R ) // within cone of truth
                {
                    min_R = dR;
                    truth_prop_dEta = truth_eta - tsos_eta;
                    truth_prop_dPhi = truth_phi - tsos_phi;
                    track_prop_dEta = track.outerEta() - tsos_eta;
                    track_prop_dPhi = track.outerPhi() - tsos_phi;
                    num_sim_clusters = hit_point.second;
                }
            }
            
        }
        
        if ( min_R < truth_matching_deltaR_ )
        {
            // simhit & propagator histograms
            TH1Container_["dEta_Truth_Prop"]->Fill(truth_prop_dEta);
            TH1Container_["dPhi_Truth_Prop"]->Fill(truth_prop_dPhi);
            TH1Container_["dR_Truth_Prop"]->Fill(min_R);
            TH2Container_["nSimCluster_dPhi_Truth_Prop"]->Fill(num_sim_clusters, truth_prop_dPhi);
            TH2Container_["nSimCluster_dEta_Truth_Prop"]->Fill(num_sim_clusters, truth_prop_dEta);
            TH2Container_["nSimCluster_dR_Truth_Prop"]->Fill(num_sim_clusters, min_R);
            
            // simhit & track histograms
            TH1Container_["dEta_Track_Prop"]->Fill(abs(track_prop_dEta));
            TH1Container_["dPhi_Track_Prop"]->Fill(abs(track_prop_dPhi));
        }
    }

}


void TICLAnalyzer::analyzeRecHits ( const std::vector<math::XYZTLorentzVectorF> & caloTruthP4, const HGCRecHitCollection & hits, const HGCalGeometry * geom )
{

    for ( auto const& truth: caloTruthP4 )
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
                int hit_layer = rhtools_.getLayer(hit_id);
                
                Float_t hit_raw_energy = hit.energy();

                hit_sum_energy += hit_raw_energy;
                hit_sum_energy_layer[hit_layer-1] += hit_raw_energy;
            }
        }
        
        TH1Container_["EDist_recHits"]->Fill ( hit_sum_energy );
        
        TH1Container_["EDist_recHits_layer1"]->Fill ( hit_sum_energy_layer[0] );
        TH1Container_["EDist_recHits_layer2"]->Fill ( hit_sum_energy_layer[1] );
        TH1Container_["EDist_recHits_layer3"]->Fill ( hit_sum_energy_layer[2] );
        TH1Container_["EDist_recHits_layer4"]->Fill ( hit_sum_energy_layer[3] );
        TH1Container_["EDist_recHits_layer5"]->Fill ( hit_sum_energy_layer[4] );
        TH1Container_["EDist_recHits_layer6"]->Fill ( hit_sum_energy_layer[5] );
        TH1Container_["EDist_recHits_layer7"]->Fill ( hit_sum_energy_layer[6] );
        TH1Container_["EDist_recHits_layer8"]->Fill ( hit_sum_energy_layer[7] );

        TH1Container_["ERatio_recHits_to_truth"]->Fill ( hit_sum_energy / truth_energy );
        
        TH1Container_["ERatio_recHits_layer1"]->Fill ( hit_sum_energy_layer[0] / truth_energy );
        TH1Container_["ERatio_recHits_layer2"]->Fill ( hit_sum_energy_layer[1] / truth_energy );
        TH1Container_["ERatio_recHits_layer3"]->Fill ( hit_sum_energy_layer[2] / truth_energy );
        TH1Container_["ERatio_recHits_layer4"]->Fill ( hit_sum_energy_layer[3] / truth_energy );
        TH1Container_["ERatio_recHits_layer5"]->Fill ( hit_sum_energy_layer[4] / truth_energy );
        TH1Container_["ERatio_recHits_layer6"]->Fill ( hit_sum_energy_layer[5] / truth_energy );
        TH1Container_["ERatio_recHits_layer7"]->Fill ( hit_sum_energy_layer[6] / truth_energy );
        TH1Container_["ERatio_recHits_layer8"]->Fill ( hit_sum_energy_layer[7] / truth_energy );
        
    }
}


void TICLAnalyzer::analyzeLayerClusters ( const std::vector<math::XYZTLorentzVectorF> & caloTruthP4, const std::vector<reco::CaloCluster> & clusters )
{

    for ( auto const& truth: caloTruthP4 )
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
            
            if ( dR < truth_matching_deltaR_ )
            {
//                for ( auto const& cluster: cl.hitsAndFractions() )
//                {
//                    std::cout << "lc layer: " << rhtools_.getLayer(cluster.first) << std::endl;
//                }
//                std::cout << "----" << std::endl;
                int cluster_layer = rhtools_.getLayer(cl.hitsAndFractions().at(0).first);
                
                Float_t cluster_raw_energy = cl.energy();
                cluster_sum_energy += cluster_raw_energy;

                if ( cluster_layer <= 8 )
                {
                    cluster_sum_energy_layer[cluster_layer-1] += cluster_raw_energy;
                    cluster_count_layer[cluster_layer-1] ++;
                }
            }
        }
        
        TH1Container_["EDist_layerClusters_sum"]->Fill ( cluster_sum_energy );
        
        TH1Container_["EDist_layerClusters_layer1"]->Fill( cluster_sum_energy_layer[0] );
        TH1Container_["EDist_layerClusters_layer2"]->Fill( cluster_sum_energy_layer[1] );
        TH1Container_["EDist_layerClusters_layer3"]->Fill( cluster_sum_energy_layer[2] );
        TH1Container_["EDist_layerClusters_layer4"]->Fill( cluster_sum_energy_layer[3] );
        TH1Container_["EDist_layerClusters_layer5"]->Fill( cluster_sum_energy_layer[4] );
        TH1Container_["EDist_layerClusters_layer6"]->Fill( cluster_sum_energy_layer[5] );
        TH1Container_["EDist_layerClusters_layer7"]->Fill( cluster_sum_energy_layer[6] );
        TH1Container_["EDist_layerClusters_layer8"]->Fill( cluster_sum_energy_layer[7] );
        
        TH1Container_["Count_layerClusters_layer1"]->Fill( cluster_count_layer[0] );
        TH1Container_["Count_layerClusters_layer2"]->Fill( cluster_count_layer[1] );
        TH1Container_["Count_layerClusters_layer3"]->Fill( cluster_count_layer[2] );
        TH1Container_["Count_layerClusters_layer4"]->Fill( cluster_count_layer[3] );
        TH1Container_["Count_layerClusters_layer5"]->Fill( cluster_count_layer[4] );
        TH1Container_["Count_layerClusters_layer6"]->Fill( cluster_count_layer[5] );
        TH1Container_["Count_layerClusters_layer7"]->Fill( cluster_count_layer[6] );
        TH1Container_["Count_layerClusters_layer8"]->Fill( cluster_count_layer[7] );
        
        TH1Container_["ERatio_layerClusters_to_truth"]->Fill ( cluster_sum_energy / truth_energy );
        
        TH1Container_["ERatio_layerClusters_layer1"]->Fill ( cluster_sum_energy_layer[0] / truth_energy );
        TH1Container_["ERatio_layerClusters_layer2"]->Fill ( cluster_sum_energy_layer[1] / truth_energy );
        TH1Container_["ERatio_layerClusters_layer3"]->Fill ( cluster_sum_energy_layer[2] / truth_energy );
        TH1Container_["ERatio_layerClusters_layer4"]->Fill ( cluster_sum_energy_layer[3] / truth_energy );
        TH1Container_["ERatio_layerClusters_layer5"]->Fill ( cluster_sum_energy_layer[4] / truth_energy );
        TH1Container_["ERatio_layerClusters_layer6"]->Fill ( cluster_sum_energy_layer[5] / truth_energy );
        TH1Container_["ERatio_layerClusters_layer7"]->Fill ( cluster_sum_energy_layer[6] / truth_energy );
        TH1Container_["ERatio_layerClusters_layer8"]->Fill ( cluster_sum_energy_layer[7] / truth_energy );
    }

}


void TICLAnalyzer::analyzeTICLTrackster ( const std::vector<math::XYZTLorentzVectorF> & caloTruthP4, const std::vector<ticl::Trackster> & tracksters, std::string tag )
{
    for ( auto const& trs: tracksters )
    {
        Float_t trackster_raw_energy = trs.raw_energy();
        Float_t trackster_eta = trs.barycenter().eta();
        Float_t trackster_phi = trs.barycenter().phi();
        //std::array<Float_t, 8> trackster_id_probs = trs.id_probabilities();
        // Float_t trackster_sigmas = trs.sigmas();
        // Float_t trackster_sigmasPCA = trs.sigmasPCA();
                
        for ( auto const& truth: caloTruthP4 )
        {
            Float_t dR = reco::deltaR ( trackster_eta, trackster_phi, truth.eta(), truth.phi() );

            if ( dR < (truth_matching_deltaR_))
            {
            
                if ( tag == "EMn" )
                {
                    TH1Container_["RawEDist_tracksterHFNoseEM"]->Fill( trackster_raw_energy );
                    TH1Container_["RawEScale_tracksterHFNoseEM"]->Fill( truth.E() - trackster_raw_energy );
                    TH1Container_["DeltaR_tracksterHFNoseEM"]->Fill( dR );
                }
                else if ( tag == "TrkEMn")
                {
                    TH1Container_["RawEDist_tracksterHFNoseTrkEM"]->Fill( trackster_raw_energy );
                    TH1Container_["RawEScale_tracksterHFNoseTrkEM"]->Fill( truth.E() - trackster_raw_energy );
                    TH1Container_["DeltaR_tracksterHFNoseTrkEM"]->Fill( dR );
                 }   
                 else if ( tag == "EM" )
                 {
                    TH1Container_["RawEDist_tracksterEM"]->Fill( trackster_raw_energy );
                    TH1Container_["RawEScale_tracksterEM"]->Fill( truth.E() - trackster_raw_energy );
                    TH1Container_["DeltaR_tracksterEM"]->Fill( dR );
                 }
                 else if ( tag == "HADn" )
                 {
                    TH1Container_["RawEDist_tracksterHFNoseHAD"]->Fill( trackster_raw_energy );
                    TH1Container_["RawEScale_tracksterHFNoseHAD"]->Fill( truth.E() - trackster_raw_energy );
                    TH1Container_["DeltaR_tracksterHFNoseHAD"]->Fill( dR );
                 }
                 else if ( tag == "MIPn" )
                 {
                    TH1Container_["RawEDist_tracksterHFNoseMIP"]->Fill( trackster_raw_energy );
                    TH1Container_["RawEScale_tracksterHFNoseMIP"]->Fill( truth.E() - trackster_raw_energy );
                    TH1Container_["DeltaR_tracksterHFNoseMIP"]->Fill( dR );
                 }
                    
            }
        }
    }
}


void TICLAnalyzer::beginJob ()
{
    edm::Service<TFileService> fs;
    
    // ---- CaloTruth ------
    TH1Container_["truthE"] = fs->make<TH1F>("truthE", "Truth Energy Distribution", 100, 0, 500);
    TH1Container_["truthEta"] = fs->make<TH1F>("truthEta", "Truth Eta Distribution", 1000, (select_EtaLow_ + select_EtaHigh_)/2 - 0.5, (select_EtaLow_ + select_EtaHigh_)/2 + 0.5);
    
    // ---- SimHits ------
    TH1Container_["first_layer_SimHit_z"] = fs->make<TH1F>("first_layer_SimHit_z", "|z| of highest energy simhit in first layer", 150, 0, 1500);
    TH1Container_["first_layer_SimHit_eta"] = fs->make<TH1F>("first_layer_SimHit_eta", "|eta| of highest energy simhit in first layer", 150, 3.0, 4.5);
    
    TH2Container_["nSimCluster_SimHitEta"] = fs->make<TH2F>("nSimCluster_SimHitEta", "# SimClusters vs |#eta_{simhit}|, per generated particle", 9, 1, 10, 150, 3.0, 4.5);
    
    TH2Container_["nSimCluster_SimHitEta"]->GetXaxis()->SetTitle("# SimClusters");
    TH2Container_["nSimCluster_SimHitEta"]->GetYaxis()->SetTitle("#eta_{simhit}");

    // ---- SimClusters ------
    TH1Container_["SimCluster_Seed_Layer"] = fs->make<TH1F>("SimCluster_First_Layer", "starting position layer of simcluster", 8, 1, 9);
    TH1Container_["nSimClusters_layer1"] = fs->make<TH1F>("nSimClusters_layer1", "number of simclusters in layer 1", 15, 1, 15);
    
    // ---- Propagator & SimHit ------
    TH1Container_["dEta_Truth_Prop"] = fs->make<TH1F>("dEta_Truth_Prop", "#eta difference between reco and sim at HGCal interface", 600, -0.3, 0.3);
    TH1Container_["dPhi_Truth_Prop"] = fs->make<TH1F>("dPhi_Truth_Prop", "#phi difference between reco and sim at HGCal interface", 600, -0.3, 0.3);
    TH1Container_["dR_Truth_Prop"] = fs->make<TH1F>("dR_Truth_Prop", "#Delta R between reco and sim at HGCal interface", 300, 0.0, 0.3);
    
    TH2Container_["nSimCluster_dEta_Truth_Prop"] = fs->make<TH2F>("nSimCluster_dEta_Truth_Prop", "# SimClusters vs (#eta_{simhit} - #eta_{prop}), per generated particle", 9, 1, 10, 60, -0.3, 0.3);
    TH2Container_["nSimCluster_dPhi_Truth_Prop"] = fs->make<TH2F>("nSimCluster_dPhi_Truth_Prop", "# SimClusters vs (#phi_{simhit} - #phi_{prop}), per generated particle", 9, 1, 10, 60, -0.3, 0.3);
    TH2Container_["nSimCluster_dR_Truth_Prop"] = fs->make<TH2F>("nSimCluster_dR_Truth_Prop", "# SimClusters vs #Delta{R}_{simhit, prop}, per generated particle", 9, 1, 10, 30, 0, 0.3);
    
    TH1Container_["dEta_Truth_Prop"]->GetXaxis()->SetTitle("#eta_{simhit} - #eta_{prop}");
    TH1Container_["dPhi_Truth_Prop"]->GetXaxis()->SetTitle("phi_{simhit} - #phi_{prop}");
    TH1Container_["dR_Truth_Prop"]->GetXaxis()->SetTitle("#Delta{R}_{simhit, prop}");
     
    TH2Container_["nSimCluster_dEta_Truth_Prop"]->GetXaxis()->SetTitle("# SimClusters");
    TH2Container_["nSimCluster_dEta_Truth_Prop"]->GetYaxis()->SetTitle("#eta_{simhit} - #eta_{prop}");
    
    TH2Container_["nSimCluster_dPhi_Truth_Prop"]->GetXaxis()->SetTitle("# SimClusters");
    TH2Container_["nSimCluster_dPhi_Truth_Prop"]->GetYaxis()->SetTitle("#phi_{simhit} - #phi_{prop}");
     
    TH2Container_["nSimCluster_dR_Truth_Prop"]->GetXaxis()->SetTitle("# SimClusters");
    TH2Container_["nSimCluster_dR_Truth_Prop"]->GetYaxis()->SetTitle("#Delta{R}_{simhit, prop}");
    
    // ---- Propagator & Track ------
    
    TH1Container_["dEta_Track_Prop"] = fs->make<TH1F>("dEta_Track_Prop", "#eta difference between #eta_{prop} and #eta_{track, outer}", 300, 0., 0.3);
    TH1Container_["dPhi_Track_Prop"] = fs->make<TH1F>("dPhi_Track_Prop", "#phi difference between #phi_{prop} and #phi_{track, outer}", 300, 0., 0.3);
    
    TH1Container_["dEta_Truth_Prop"]->GetXaxis()->SetTitle("#eta_{simhit} - #eta_{prop}");
    TH1Container_["dPhi_Truth_Prop"]->GetXaxis()->SetTitle("phi_{simhit} - #phi_{prop}");   
    
    // ---- RecHits ------
    TH1Container_["EDist_recHits"] = fs->make<TH1F>("EDist_recHits", "Energy Distribution (recHits)", 100, 0, 1000);
    
    TH1Container_["EDist_recHits_layer1"] = fs->make<TH1F>("EDist_recHits_layer1", "Energy Distribution Layer 1 (recHits)", 100, 0, 1000);
    TH1Container_["EDist_recHits_layer2"] = fs->make<TH1F>("EDist_recHits_layer2", "Energy Distribution Layer 2 (recHits)", 100, 0, 1000);
    TH1Container_["EDist_recHits_layer3"] = fs->make<TH1F>("EDist_recHits_layer3", "Energy Distribution Layer 3 (recHits)", 100, 0, 1000);
    TH1Container_["EDist_recHits_layer4"] = fs->make<TH1F>("EDist_recHits_layer4", "Energy Distribution Layer 4 (recHits)", 100, 0, 1000);
    TH1Container_["EDist_recHits_layer5"] = fs->make<TH1F>("EDist_recHits_layer5", "Energy Distribution Layer 5 (recHits)", 100, 0, 1000);
    TH1Container_["EDist_recHits_layer6"] = fs->make<TH1F>("EDist_recHits_layer6", "Energy Distribution Layer 6 (recHits)", 100, 0, 1000);
    TH1Container_["EDist_recHits_layer7"] = fs->make<TH1F>("EDist_recHits_layer7", "Energy Distribution Layer 7 (recHits)", 100, 0, 1000);
    TH1Container_["EDist_recHits_layer8"] = fs->make<TH1F>("EDist_recHits_layer8", "Energy Distribution Layer 8 (recHits)", 100, 0, 1000);
    
    TH1Container_["ERatio_recHits_to_truth"] = fs->make<TH1F>("ERatio_recHits_to_truth", "E_{recHits} / E_{truth}", 100, 0, 1.0);
    TH1Container_["ERatio_recHits_to_truth"]->GetXaxis()->SetTitle("E_{recHits} / E_{truth}");

    TH1Container_["ERatio_recHits_layer1"] = fs->make<TH1F>("ERatio_recHits_layer1", "E_{recHits} / E_{truth} Layer 1", 100, 0, 1.0);
    TH1Container_["ERatio_recHits_layer2"] = fs->make<TH1F>("ERatio_recHits_layer2", "E_{recHits} / E_{truth} Layer 2", 100, 0, 1.0);
    TH1Container_["ERatio_recHits_layer3"] = fs->make<TH1F>("ERatio_recHits_layer3", "E_{recHits} / E_{truth} Layer 3", 100, 0, 1.0);
    TH1Container_["ERatio_recHits_layer4"] = fs->make<TH1F>("ERatio_recHits_layer4", "E_{recHits} / E_{truth} Layer 4", 100, 0, 1.0);
    TH1Container_["ERatio_recHits_layer5"] = fs->make<TH1F>("ERatio_recHits_layer5", "E_{recHits} / E_{truth} Layer 5", 100, 0, 1.0);
    TH1Container_["ERatio_recHits_layer6"] = fs->make<TH1F>("ERatio_recHits_layer6", "E_{recHits} / E_{truth} Layer 6", 100, 0, 1.0);
    TH1Container_["ERatio_recHits_layer7"] = fs->make<TH1F>("ERatio_recHits_layer7", "E_{recHits} / E_{truth} Layer 7", 100, 0, 1.0);
    TH1Container_["ERatio_recHits_layer8"] = fs->make<TH1F>("ERatio_recHits_layer8", "E_{recHits} / E_{truth} Layer 8", 100, 0, 1.0);
    
    // ---- Clusters ------
    TH1Container_["EDist_layerClusters_sum"] = fs->make<TH1F>("EDist_layerClusters_sum", "Energy Distribution (clusters, summed all layers)", 100, 0, 1000);
    
    TH1Container_["ERatio_layerClusters_to_truth"] = fs->make<TH1F>("ERatio_layerClusters_to_truth", "E_{layerClusters} / E_{truth}", 100, 0, 1.0);
    TH1Container_["ERatio_layerClusters_to_truth"]->GetXaxis()->SetTitle("E_{layerClusters} / E_{truth}");
    
    TH1Container_["EDist_layerClusters_layer1"] = fs->make<TH1F>("EDist_layerClusters_layer1", "Energy Distribution Layer 1 (clusters)", 100, 0, 1000);
    TH1Container_["EDist_layerClusters_layer2"] = fs->make<TH1F>("EDist_layerClusters_layer2", "Energy Distribution Layer 2 (clusters)", 100, 0, 1000);
    TH1Container_["EDist_layerClusters_layer3"] = fs->make<TH1F>("EDist_layerClusters_layer3", "Energy Distribution Layer 3 (clusters)", 100, 0, 1000);
    TH1Container_["EDist_layerClusters_layer4"] = fs->make<TH1F>("EDist_layerClusters_layer4", "Energy Distribution Layer 4 (clusters)", 100, 0, 1000);
    TH1Container_["EDist_layerClusters_layer5"] = fs->make<TH1F>("EDist_layerClusters_layer5", "Energy Distribution Layer 5 (clusters)", 100, 0, 1000);
    TH1Container_["EDist_layerClusters_layer6"] = fs->make<TH1F>("EDist_layerClusters_layer6", "Energy Distribution Layer 6 (clusters)", 100, 0, 1000);
    TH1Container_["EDist_layerClusters_layer7"] = fs->make<TH1F>("EDist_layerClusters_layer7", "Energy Distribution Layer 7 (clusters)", 100, 0, 1000);
    TH1Container_["EDist_layerClusters_layer8"] = fs->make<TH1F>("EDist_layerClusters_layer8", "Energy Distribution Layer 8 (clusters)", 100, 0, 1000);
    
    TH1Container_["Count_layerClusters_layer1"] = fs->make<TH1F>("Count_layerClusters_layer1", "Number of Clusters Layer 1", 50, 0, 50);
    TH1Container_["Count_layerClusters_layer2"] = fs->make<TH1F>("Count_layerClusters_layer2", "Number of Clusters Layer 2", 50, 0, 50);
    TH1Container_["Count_layerClusters_layer3"] = fs->make<TH1F>("Count_layerClusters_layer3", "Number of Clusters Layer 3", 50, 0, 50);
    TH1Container_["Count_layerClusters_layer4"] = fs->make<TH1F>("Count_layerClusters_layer4", "Number of Clusters Layer 4", 50, 0, 50);
    TH1Container_["Count_layerClusters_layer5"] = fs->make<TH1F>("Count_layerClusters_layer5", "Number of Clusters Layer 5", 50, 0, 50);
    TH1Container_["Count_layerClusters_layer6"] = fs->make<TH1F>("Count_layerClusters_layer6", "Number of Clusters Layer 6", 50, 0, 50);
    TH1Container_["Count_layerClusters_layer7"] = fs->make<TH1F>("Count_layerClusters_layer7", "Number of Clusters Layer 7", 50, 0, 50);
    TH1Container_["Count_layerClusters_layer8"] = fs->make<TH1F>("Count_layerClusters_layer8", "Number of Clusters Layer 8", 50, 0, 50);
    
    TH1Container_["ERatio_layerClusters_layer1"] = fs->make<TH1F>("ERatio_layerClusters_layer1", "E_{layerClusters} / E_{truth} Layer 1", 100, 0, 1.0);
    TH1Container_["ERatio_layerClusters_layer2"] = fs->make<TH1F>("ERatio_layerClusters_layer2", "E_{layerClusters} / E_{truth} Layer 2", 100, 0, 1.0);
    TH1Container_["ERatio_layerClusters_layer3"] = fs->make<TH1F>("ERatio_layerClusters_layer3", "E_{layerClusters} / E_{truth} Layer 3", 100, 0, 1.0);
    TH1Container_["ERatio_layerClusters_layer4"] = fs->make<TH1F>("ERatio_layerClusters_layer4", "E_{layerClusters} / E_{truth} Layer 4", 100, 0, 1.0);
    TH1Container_["ERatio_layerClusters_layer5"] = fs->make<TH1F>("ERatio_layerClusters_layer5", "E_{layerClusters} / E_{truth} Layer 5", 100, 0, 1.0);
    TH1Container_["ERatio_layerClusters_layer6"] = fs->make<TH1F>("ERatio_layerClusters_layer6", "E_{layerClusters} / E_{truth} Layer 6", 100, 0, 1.0);
    TH1Container_["ERatio_layerClusters_layer7"] = fs->make<TH1F>("ERatio_layerClusters_layer7", "E_{layerClusters} / E_{truth} Layer 7", 100, 0, 1.0);
    TH1Container_["ERatio_layerClusters_layer8"] = fs->make<TH1F>("ERatio_layerClusters_layer8", "E_{layerClusters} / E_{truth} Layer 8", 100, 0, 1.0);
    
    // ---- Tracksters ------
    TH1Container_["RawEDist_tracksterHFNoseEM"] = fs->make<TH1F>("EDist_tracksterHFNoseEM", "TracksterHFNoseEM Raw Energy Distribution", 100, 0, 500);
    TH1Container_["RawEScale_tracksterHFNoseEM"] = fs->make<TH1F>("EScale_tracksterHFNoseEM", "TracksterHFNoseEM Raw Energy Scale", 100, 0, 500);
    TH1Container_["DeltaR_tracksterHFNoseEM"] = fs->make<TH1F>("DeltaR_tracksterHFNoseEM", "TracksterHFNoseEM #Delta R", 20, 0, 0.5);
    
    TH1Container_["RawEDist_tracksterHFNoseEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    TH1Container_["RawEScale_tracksterHFNoseEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    TH1Container_["DeltaR_tracksterHFNoseEM"]->GetXaxis()->SetTitle("#Delta{R}_{trackster, caloParticle}");
    
    //
    TH1Container_["RawEDist_tracksterHFNoseTrkEM"] = fs->make<TH1F>("EDist_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM Raw Energy Distribution", 100, 0, 500);
    TH1Container_["RawEScale_tracksterHFNoseTrkEM"] = fs->make<TH1F>("EScale_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM Raw Energy Scale", 100, 0, 500);
    TH1Container_["DeltaR_tracksterHFNoseTrkEM"] = fs->make<TH1F>("DeltaR_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM #Delta R", 20, 0, 0.5);
    
    TH1Container_["RawEDist_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    TH1Container_["RawEScale_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    TH1Container_["DeltaR_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("#Delta{R}_{trackster, caloParticle}");
    
    //
    TH1Container_["RawEDist_tracksterEM"] = fs->make<TH1F>("EDist_tracksterEM", "TracksterEM Raw Energy Distribution", 500, 0, 500);
    TH1Container_["RawEScale_tracksterEM"] = fs->make<TH1F>("EScale_tracksterEM", "TracksterEM Raw Energy Scale", 100, 0, 500);
    TH1Container_["DeltaR_tracksterEM"] = fs->make<TH1F>("DeltaR_tracksterEM", "TracksterEM #Delta R", 20, 0, 0.5);
    
    TH1Container_["RawEDist_tracksterEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    TH1Container_["RawEScale_tracksterEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    TH1Container_["DeltaR_tracksterEM"]->GetXaxis()->SetTitle("#Delta{R}_{trackster, caloParticle}");
    
    //
    TH1Container_["RawEDist_tracksterHFNoseHAD"] = fs->make<TH1F>("EDist_tracksterHFNoseHAD", "TracksterHFNoseHAD Raw Energy Distribution", 500, 0, 500);
    TH1Container_["RawEScale_tracksterHFNoseHAD"] = fs->make<TH1F>("EScale_tracksterHFNoseHAD", "TracksterHFNoseHAD Raw Energy Scale", 100, 0, 500);
    TH1Container_["DeltaR_tracksterHFNoseHAD"] = fs->make<TH1F>("DeltaR_tracksterHFNoseHAD", "TracksterHFNoseHAD #Delta R", 20, 0, 0.5);
    
    TH1Container_["RawEDist_tracksterHFNoseHAD"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    TH1Container_["RawEScale_tracksterHFNoseHAD"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    TH1Container_["DeltaR_tracksterHFNoseHAD"]->GetXaxis()->SetTitle("#Delta{R}_{trackster, caloParticle}");
    
    //
    TH1Container_["RawEDist_tracksterHFNoseMIP"] = fs->make<TH1F>("EDist_tracksterHFNoseMIP", "TracksterHFNoseMIP Raw Energy Distribution", 500, 0, 500);
    TH1Container_["RawEScale_tracksterHFNoseMIP"] = fs->make<TH1F>("EScale_tracksterHFNoseMIP", "TracksterHFNoseMIP Raw Energy Scale", 100, 0, 500);
    TH1Container_["DeltaR_tracksterHFNoseMIP"] = fs->make<TH1F>("DeltaR_tracksterHFNoseMIP", "TracksterHFNoseMIP #Delta R", 20, 0, 0.5);
    
    TH1Container_["RawEDist_tracksterHFNoseMIP"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    TH1Container_["RawEScale_tracksterHFNoseMIP"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    TH1Container_["DeltaR_tracksterHFNoseMIP"]->GetXaxis()->SetTitle("#Delta{R}_{trackster, caloParticle}");

}


void TICLAnalyzer::endJob ()
{
    // Job ends
    std::cout << "The job, TICLAnalyzer, has ended. Thank you for your patience." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (TICLAnalyzer);
