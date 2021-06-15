#include "TICLAnalyzer.h"

#include <iostream>
#include <array>
#include <algorithm>
#include <string>

#include <cmath>

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
    histContainer_ (),
    tag_CaloParticle_MergedCaloTruth_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_MergedCaloTruth", edm::InputTag ("mix", "MergedCaloTruth", "HLT") ) ),
    tag_Tracks_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_Tracks", edm::InputTag("generalTracks") ) ),
    tag_RecHits_HFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_HGCHFNoseRecHits", edm::InputTag("HGCalRecHit:HGCHFNoseRecHits") ) ),
    tag_LayerClusters_HFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_LayerClustersHFNose", edm::InputTag("hgcalLayerClustersHFNose") ) ),
    tag_Trackster_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_Trackster", edm::InputTag ("ticlTrackstersHFNoseEM") ) ),
    
    // Custom parameters
    select_PID_ ( iConfig.getUntrackedParameter<int> ("select_PID", 22) ),
    select_EtaLow_ ( iConfig.getUntrackedParameter<double> ("select_EtaLow", 3.4) ),
    select_EtaHigh_ ( iConfig.getUntrackedParameter<double> ("select_EtaHigh", 3.6) ),
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
    edm::ESHandle<Propagator> handle_Propagator;
    iSetup.get<TrackingComponentsRecord>().get( "RungeKuttaTrackerPropagator", handle_Propagator );
    
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
        const std::vector<GlobalPoint> simhits_globalpoints = getSimHitGlobalPoint ( *handle_CaloParticle_MergedCaloTruth.product(), nose_geom );
        fillTruthHistograms ( selected_calotruths, simhits_globalpoints );
        
        buildFirstLayers ( hgcons );
        analyzeTrackPosition ( simhits_globalpoints, *handle_Tracks.product(), b_field, *handle_Propagator );
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


std::vector<GlobalPoint> TICLAnalyzer::getSimHitGlobalPoint ( const std::vector<CaloParticle> & caloTruthParticles, const HGCalGeometry * geom )
{

    std::vector<GlobalPoint> container;
    
    for ( auto const& ct: caloTruthParticles )
    {
        if ( abs(ct.pdgId()) == select_PID_
                && abs(ct.eta()) > select_EtaLow_
                && abs(ct.eta()) < select_EtaHigh_ )
        {
            const SimClusterRefVector & ct_sim_clusters = ct.simClusters();
            
            Double_t simhit_max_energy = 0.;
            DetId simhit_detid;
            
//            Double_t simhit_second_max_energy = 0.;
//            DetId simhit_second_detid;
//            
//            Double_t simhit_third_max_energy = 0.;
//            DetId simhit_third_detid;
            
            for ( auto const& it_sc: ct_sim_clusters )
            {
                const SimCluster sc = *it_sc;
                for ( auto const & he: sc.hits_and_energies() )
                {
                    if ( rhtools_.getLayer(he.first) != 1 ) continue;
                    if ( simhit_max_energy < he.second ) 
                    {
                        simhit_max_energy = he.second;
                        simhit_detid = he.first;
                    }
//                    else if ( simhit_second_max_energy < he.second )
//                    {
//                        simhit_second_max_energy = he.second;
//                        simhit_second_detid = he.first;
//                    }
//                    else if ( simhit_third_max_energy < he.second )
//                    {
//                        simhit_third_max_energy = he.second;
//                        simhit_third_detid = he.first;
//                    }
                }
            }
            const GlobalPoint & hit_globalPosition = geom->getPosition(simhit_detid);
//            const GlobalPoint & hit_second_globalPosition = geom->getPosition(simhit_second_detid);
//            const GlobalPoint & hit_third_globalPosition = geom->getPosition(simhit_third_detid);
//            std::cout << "simhit max E: " << simhit_max_energy << std::endl;
//            std::cout << "simhit z: " << hit_globalPosition.z() << std::endl;
//            std::cout << "simhit eta: " << hit_globalPosition.eta() << std::endl;
//            std::cout << "simhit phi: " << hit_globalPosition.phi() << std::endl;
//            std::cout << "simhit second max E: " << simhit_second_max_energy << std::endl;
//            std::cout << "simhit second z: " << hit_second_globalPosition.z() << std::endl;
//            std::cout << "simhit second eta: " << hit_second_globalPosition.eta() << std::endl;
//            std::cout << "simhit second phi: " << hit_second_globalPosition.phi() << std::endl;
//            std::cout << "simhit third max E: " << simhit_third_max_energy << std::endl;
//            std::cout << "simhit third z: " << hit_third_globalPosition.z() << std::endl;
//            std::cout << "simhit third eta: " << hit_third_globalPosition.eta() << std::endl;
//            std::cout << "simhit third phi: " << hit_third_globalPosition.phi() << std::endl;
            container.push_back ( hit_globalPosition );
        }
    }
    return container;
}


void TICLAnalyzer::fillTruthHistograms ( const std::vector<math::XYZTLorentzVectorF> & caloTruthP4, const std::vector<GlobalPoint> & simHitsGlobalPoints )
{
    for ( auto const& truth: caloTruthP4 )
    {
        // Fill truth histograms
        histContainer_["truthE"]->Fill( truth.energy() );
        histContainer_["truthEta"]->Fill( truth.eta() );   
    }
    for ( auto const& hit_point: simHitsGlobalPoints )
    {
        histContainer_["first_layer_SimHit_z"]->Fill( hit_point.z() );
        histContainer_["first_layer_SimHit_eta"]->Fill( hit_point.eta() );
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
        
//        std::cout << truth.z() << std::endl;
        
        Float_t min_R = 9999.;
        Float_t truth_track_dEta = -1.;
        Float_t truth_track_dPhi = -1.;
        
        for ( auto const& track: generalTracks )
        {
            if ( !cutTk_(track) )
            {
                continue;
            }
        
            FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(track, BField);
            int iSide = int(track.eta() > 0);
            auto begin = std::chrono::high_resolution_clock::now();
            TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
            auto end = std::chrono::high_resolution_clock::now();
            
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            
            std::cout << tsos.isValid() << std::endl;
            if (tsos.isValid())
            {
//                std::cout.precision(6);
//                std::cout << "surface z: " << firstDisk_[iSide]->surface().position().z() << "\n"
//                        << "track outer eta: " << track.outerEta() << "\n"
//                        << "track outer phi: " << track.outerPhi() << "\n"
//                        << "tsos eta: "   << tsos.globalPosition().eta() << "\n"
//                        << "tsos phi: " << tsos.globalPosition().phi() << "\n"
//                        << "tsos r: "   << tsos.globalPosition().perp() << "\n"
//                        << "tsos z: "   << tsos.globalPosition().z() << "\n"
//                        << "B Field nominalValue: " << BField->nominalValue() << "\n"
//                        << "-------" << std::endl;
            
                Float_t tsos_eta = tsos.globalPosition().eta();
                Float_t tsos_phi = tsos.globalPosition().phi();
                
                Float_t dR = reco::deltaR ( tsos_eta, tsos_phi, truth_eta, truth_phi );

                if ( dR < truth_matching_deltaR_ && dR < min_R ) // within cone of truth
                {
                    min_R = dR;
                    truth_track_dEta = abs(truth_eta - tsos_eta);
                    truth_track_dPhi = abs(truth_phi - tsos_phi);
                }
            }
            
        }
        
        if ( min_R < truth_matching_deltaR_ )
        {
            histContainer_["dEta_Truth_Track"]->Fill(truth_track_dEta);
            histContainer_["dPhi_Truth_Track"]->Fill(truth_track_dPhi);
        }
    }

}


void TICLAnalyzer::analyzeTrackPosition ( const std::vector<GlobalPoint> & simHitsGlobalPoints, const std::vector<reco::Track> & generalTracks, const MagneticField * BField, const Propagator & prop )
{

    for ( auto const& hit_point: simHitsGlobalPoints )
    {
        Float_t truth_eta = hit_point.eta();
        Float_t truth_phi = hit_point.phi();
        
        Float_t min_R = 9999.;
        Float_t truth_track_dEta = -1.;
        Float_t truth_track_dPhi = -1.;
        
        for ( auto const& track: generalTracks )
        {
            if ( !cutTk_(track) )
            {
                continue;
            }
        
            FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(track, BField);
            int iSide = int(track.eta() > 0);
            TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
            
            std::cout << tsos.isValid() << std::endl;
            if (tsos.isValid())
            {            
                Float_t tsos_eta = tsos.globalPosition().eta();
                Float_t tsos_phi = tsos.globalPosition().phi();
                
                Float_t dR = reco::deltaR ( tsos_eta, tsos_phi, truth_eta, truth_phi );

                if ( dR < truth_matching_deltaR_ && dR < min_R ) // within cone of truth
                {
                    min_R = dR;
                    truth_track_dEta = abs(truth_eta - tsos_eta);
                    truth_track_dPhi = abs(truth_phi - tsos_phi);
                }
            }
            
        }
        
        if ( min_R < truth_matching_deltaR_ )
        {
            histContainer_["dEta_Truth_Track"]->Fill(truth_track_dEta);
            histContainer_["dPhi_Truth_Track"]->Fill(truth_track_dPhi);
            histContainer_["dR_Truth_Track"]->Fill(min_R);
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
        
        histContainer_["EDist_layerClusters_sum"]->Fill ( cluster_sum_energy );
        
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
                 else if ( tag == "HADn" )
                 {
                    histContainer_["RawEDist_tracksterHFNoseHAD"]->Fill( trackster_raw_energy );
                    histContainer_["RawEScale_tracksterHFNoseHAD"]->Fill( truth.E() - trackster_raw_energy );
                    histContainer_["DeltaR_tracksterHFNoseHAD"]->Fill( dR );
                 }
                 else if ( tag == "MIPn" )
                 {
                    histContainer_["RawEDist_tracksterHFNoseMIP"]->Fill( trackster_raw_energy );
                    histContainer_["RawEScale_tracksterHFNoseMIP"]->Fill( truth.E() - trackster_raw_energy );
                    histContainer_["DeltaR_tracksterHFNoseMIP"]->Fill( dR );
                 }
                    
            }
        }
    }
}


void TICLAnalyzer::beginJob ()
{
    edm::Service<TFileService> fs;
    
    // ---- CaloTruth ------
    histContainer_["truthE"] = fs->make<TH1F>("truthE", "Truth Energy Distribution", 100, 0, 500);
    histContainer_["truthEta"] = fs->make<TH1F>("truthEta", "Truth Eta Distribution", 40, (select_EtaLow_ + select_EtaHigh_)/2 - 0.5, (select_EtaLow_ + select_EtaHigh_)/2 + 0.5);
    
    // ---- SimHits ------
     histContainer_["first_layer_SimHit_z"] = fs->make<TH1F>("first_layer_SimHit_z", "z of highest energy simhit in first layer", 0, 1500, 150);
      histContainer_["first_layer_SimHit_eta"] = fs->make<TH1F>("first_layer_SimHit_eta", "eta of highest energy simhit in first layer", 3.0, 4.5, 150);
    
    // ---- Tracks -------
    histContainer_["dEta_Truth_Track"] = fs->make<TH1F>("dEta_Truth_Track", "#eta difference between reco and sim at HGCal interface", 30, 0, 0.3);
    histContainer_["dPhi_Truth_Track"] = fs->make<TH1F>("dPhi_Truth_Track", "#phi difference between reco and sim at HGCal interface", 30, 0, 0.3);
    histContainer_["dR_Truth_Track"] = fs->make<TH1F>("dR_Truth_Track", "#Delta R between reco and sim at HGCal interface", 30, 0, 0.3);
    
    histContainer_["dEta_Truth_Track"]->GetXaxis()->SetTitle("|#eta_{sim} - #eta_{reco}|");
    histContainer_["dPhi_Truth_Track"]->GetXaxis()->SetTitle("|#phi_{sim} - #phi_{reco}|");
    histContainer_["dR_Truth_Track"]->GetXaxis()->SetTitle("|R_{sim} - R_{reco}|");
    
    // ---- RecHits ------
    histContainer_["EDist_recHits"] = fs->make<TH1F>("EDist_recHits", "Energy Distribution (recHits)", 100, 0, 1000);
    
    histContainer_["EDist_recHits_layer1"] = fs->make<TH1F>("EDist_recHits_layer1", "Energy Distribution Layer 1 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer2"] = fs->make<TH1F>("EDist_recHits_layer2", "Energy Distribution Layer 2 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer3"] = fs->make<TH1F>("EDist_recHits_layer3", "Energy Distribution Layer 3 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer4"] = fs->make<TH1F>("EDist_recHits_layer4", "Energy Distribution Layer 4 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer5"] = fs->make<TH1F>("EDist_recHits_layer5", "Energy Distribution Layer 5 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer6"] = fs->make<TH1F>("EDist_recHits_layer6", "Energy Distribution Layer 6 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer7"] = fs->make<TH1F>("EDist_recHits_layer7", "Energy Distribution Layer 7 (recHits)", 100, 0, 1000);
    histContainer_["EDist_recHits_layer8"] = fs->make<TH1F>("EDist_recHits_layer8", "Energy Distribution Layer 8 (recHits)", 100, 0, 1000);
    
    histContainer_["ERatio_recHits_to_truth"] = fs->make<TH1F>("ERatio_recHits_to_truth", "E_{recHits} / E_{truth}", 100, 0, 1.0);
    histContainer_["ERatio_recHits_to_truth"]->GetXaxis()->SetTitle("E_{recHits} / E_{truth}");

    histContainer_["ERatio_recHits_layer1"] = fs->make<TH1F>("ERatio_recHits_layer1", "E_{recHits} / E_{truth} Layer 1", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer2"] = fs->make<TH1F>("ERatio_recHits_layer2", "E_{recHits} / E_{truth} Layer 2", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer3"] = fs->make<TH1F>("ERatio_recHits_layer3", "E_{recHits} / E_{truth} Layer 3", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer4"] = fs->make<TH1F>("ERatio_recHits_layer4", "E_{recHits} / E_{truth} Layer 4", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer5"] = fs->make<TH1F>("ERatio_recHits_layer5", "E_{recHits} / E_{truth} Layer 5", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer6"] = fs->make<TH1F>("ERatio_recHits_layer6", "E_{recHits} / E_{truth} Layer 6", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer7"] = fs->make<TH1F>("ERatio_recHits_layer7", "E_{recHits} / E_{truth} Layer 7", 100, 0, 1.0);
    histContainer_["ERatio_recHits_layer8"] = fs->make<TH1F>("ERatio_recHits_layer8", "E_{recHits} / E_{truth} Layer 8", 100, 0, 1.0);
    
    // ---- Clusters ------
    histContainer_["EDist_layerClusters_sum"] = fs->make<TH1F>("EDist_layerClusters_sum", "Energy Distribution (clusters, summed all layers)", 100, 0, 1000);
    
    histContainer_["ERatio_layerClusters_to_truth"] = fs->make<TH1F>("ERatio_layerClusters_to_truth", "E_{layerClusters} / E_{truth}", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_to_truth"]->GetXaxis()->SetTitle("E_{layerClusters} / E_{truth}");
    
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
    
    histContainer_["ERatio_layerClusters_layer1"] = fs->make<TH1F>("ERatio_layerClusters_layer1", "E_{layerClusters} / E_{truth} Layer 1", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer2"] = fs->make<TH1F>("ERatio_layerClusters_layer2", "E_{layerClusters} / E_{truth} Layer 2", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer3"] = fs->make<TH1F>("ERatio_layerClusters_layer3", "E_{layerClusters} / E_{truth} Layer 3", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer4"] = fs->make<TH1F>("ERatio_layerClusters_layer4", "E_{layerClusters} / E_{truth} Layer 4", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer5"] = fs->make<TH1F>("ERatio_layerClusters_layer5", "E_{layerClusters} / E_{truth} Layer 5", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer6"] = fs->make<TH1F>("ERatio_layerClusters_layer6", "E_{layerClusters} / E_{truth} Layer 6", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer7"] = fs->make<TH1F>("ERatio_layerClusters_layer7", "E_{layerClusters} / E_{truth} Layer 7", 100, 0, 1.0);
    histContainer_["ERatio_layerClusters_layer8"] = fs->make<TH1F>("ERatio_layerClusters_layer8", "E_{layerClusters} / E_{truth} Layer 8", 100, 0, 1.0);
    
    // ---- Tracksters ------
    histContainer_["RawEDist_tracksterHFNoseEM"] = fs->make<TH1F>("EDist_tracksterHFNoseEM", "TracksterHFNoseEM Raw Energy Distribution", 100, 0, 500);
    histContainer_["RawEScale_tracksterHFNoseEM"] = fs->make<TH1F>("EScale_tracksterHFNoseEM", "TracksterHFNoseEM Raw Energy Scale", 100, 0, 500);
    histContainer_["DeltaR_tracksterHFNoseEM"] = fs->make<TH1F>("DeltaR_tracksterHFNoseEM", "TracksterHFNoseEM #Delta R", 20, 0, 0.5);
    
    histContainer_["RawEDist_tracksterHFNoseEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    histContainer_["RawEScale_tracksterHFNoseEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    histContainer_["DeltaR_tracksterHFNoseEM"]->GetXaxis()->SetTitle("|R_{trackster - caloParticle}|");
    
    //
    histContainer_["RawEDist_tracksterHFNoseTrkEM"] = fs->make<TH1F>("EDist_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM Raw Energy Distribution", 100, 0, 500);
    histContainer_["RawEScale_tracksterHFNoseTrkEM"] = fs->make<TH1F>("EScale_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM Raw Energy Scale", 100, 0, 500);
    histContainer_["DeltaR_tracksterHFNoseTrkEM"] = fs->make<TH1F>("DeltaR_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM #Delta R", 20, 0, 0.5);
    
    histContainer_["RawEDist_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    histContainer_["RawEScale_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    histContainer_["DeltaR_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("|R_{trackster} - R_{caloParticle}|");
    
    //
    histContainer_["RawEDist_tracksterEM"] = fs->make<TH1F>("EDist_tracksterEM", "TracksterEM Raw Energy Distribution", 500, 0, 500);
    histContainer_["RawEScale_tracksterEM"] = fs->make<TH1F>("EScale_tracksterEM", "TracksterEM Raw Energy Scale", 100, 0, 500);
    histContainer_["DeltaR_tracksterEM"] = fs->make<TH1F>("DeltaR_tracksterEM", "TracksterEM #Delta R", 20, 0, 0.5);
    
    histContainer_["RawEDist_tracksterEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    histContainer_["RawEScale_tracksterEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    histContainer_["DeltaR_tracksterEM"]->GetXaxis()->SetTitle("|R_{trackster} - R_{caloParticle}|");
    
    //
    histContainer_["RawEDist_tracksterHFNoseHAD"] = fs->make<TH1F>("EDist_tracksterHFNoseHAD", "TracksterHFNoseHAD Raw Energy Distribution", 500, 0, 500);
    histContainer_["RawEScale_tracksterHFNoseHAD"] = fs->make<TH1F>("EScale_tracksterHFNoseHAD", "TracksterHFNoseHAD Raw Energy Scale", 100, 0, 500);
    histContainer_["DeltaR_tracksterHFNoseHAD"] = fs->make<TH1F>("DeltaR_tracksterHFNoseHAD", "TracksterHFNoseHAD #Delta R", 20, 0, 0.5);
    
    histContainer_["RawEDist_tracksterHFNoseHAD"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    histContainer_["RawEScale_tracksterHFNoseHAD"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    histContainer_["DeltaR_tracksterHFNoseHAD"]->GetXaxis()->SetTitle("|R_{trackster} - R_{caloParticle}|");
    
    //
    histContainer_["RawEDist_tracksterHFNoseMIP"] = fs->make<TH1F>("EDist_tracksterHFNoseMIP", "TracksterHFNoseMIP Raw Energy Distribution", 500, 0, 500);
    histContainer_["RawEScale_tracksterHFNoseMIP"] = fs->make<TH1F>("EScale_tracksterHFNoseMIP", "TracksterHFNoseMIP Raw Energy Scale", 100, 0, 500);
    histContainer_["DeltaR_tracksterHFNoseMIP"] = fs->make<TH1F>("DeltaR_tracksterHFNoseMIP", "TracksterHFNoseMIP #Delta R", 20, 0, 0.5);
    
    histContainer_["RawEDist_tracksterHFNoseMIP"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    histContainer_["RawEScale_tracksterHFNoseMIP"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    histContainer_["DeltaR_tracksterHFNoseMIP"]->GetXaxis()->SetTitle("|R_{trackster} - R_{caloParticle}|");

}


void TICLAnalyzer::endJob ()
{
    // Job ends
    std::cout << "The job, TICLAnalyzer, has ended. Thank you for your patience." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (TICLAnalyzer);
