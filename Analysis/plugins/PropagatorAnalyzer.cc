#include "PropagatorAnalyzer.h"

#include <iostream>
#include <string>

#include <cmath>

#include "TH1F.h"
#include "TH2F.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

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

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

using namespace edm;

PropagatorAnalyzer::PropagatorAnalyzer ( const edm::ParameterSet& iConfig ) :
    TH1Container_ (),
    TH2Container_ (),
    tag_CaloParticle_MergedCaloTruth_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_MergedCaloTruth", edm::InputTag ("mix", "MergedCaloTruth", "HLT") ) ),
    tag_Tracks_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_Tracks", edm::InputTag("generalTracks") ) ),
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
    
    token_TracksterHFNoseEM_ = consumes<std::vector<ticl::Trackster>> ( edm::InputTag ("ticlTrackstersHFNoseEM") );
    token_TracksterHFNoseTrkEM_ = consumes<std::vector<ticl::Trackster>> ( edm::InputTag ("ticlTrackstersHFNoseTrkEM") );
    token_TracksterHFNoseHAD_ = consumes<std::vector<ticl::Trackster>> ( edm::InputTag ("ticlTrackstersHFNoseHAD") );
    token_TracksterHFNoseMIP_ = consumes<std::vector<ticl::Trackster>> ( edm::InputTag ("ticlTrackstersHFNoseMIP") );
}


PropagatorAnalyzer::~PropagatorAnalyzer ()
{
    // Deconstructor
}


void PropagatorAnalyzer::analyze ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
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
    
    // Get Propagators
    edm::ESHandle<Propagator> handle_RKPropagator;
    iSetup.get<TrackingComponentsRecord>().get( "RungeKuttaTrackerPropagator", handle_RKPropagator );
    
    edm::ESHandle<Propagator> handle_SteppingHelixPropagator;
    iSetup.get<TrackingComponentsRecord>().get( "SteppingHelixPropagatorAny", handle_SteppingHelixPropagator );
    
    edm::ESHandle<Propagator> handle_AnalyticalPropagator;
    iSetup.get<TrackingComponentsRecord>().get( "AnalyticalPropagator", handle_AnalyticalPropagator );
    
    // Get Tracks
    edm::Handle<std::vector<reco::Track>> handle_Tracks;
    iEvent.getByToken ( token_Tracks_, handle_Tracks );

    // Get CaloTruth
    edm::Handle<std::vector<CaloParticle>> handle_CaloParticle_MergedCaloTruth;
    iEvent.getByToken ( token_CaloParticle_MergedCaloTruth_, handle_CaloParticle_MergedCaloTruth );
    
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
         handle_TracksterHFNoseEM.isValid() &&
         handle_TracksterHFNoseTrkEM.isValid() &&
         handle_TracksterHFNoseHAD.isValid() &&
         handle_TracksterHFNoseMIP.isValid() &&
         handle_AnalyticalPropagator.isValid() &&
         handle_RKPropagator.isValid() &&
         handle_SteppingHelixPropagator.isValid()
       )
    {
        const HGCalGeometry* nose_geom = handle_HGCalGeometry_HFNose.product();
        const MagneticField* b_field = handle_MagneticField.product();
        const HGCalDDDConstants* hgcons = handle_HGCalDDDConstants.product();
        
        const std::vector<reco::Track> generalTracks = *handle_Tracks.product();
        const std::vector<CaloParticle> caloTruthParticles = *handle_CaloParticle_MergedCaloTruth.product();
        
        const std::vector<math::XYZTLorentzVectorF> selected_calotruths = getTruthP4 ( caloTruthParticles );
        const std::vector<GlobalPoint> simhits_globalpoints = getSimHitGlobalPoint ( caloTruthParticles, nose_geom );
        fillTruthHistograms ( selected_calotruths, simhits_globalpoints );
        
        buildFirstLayersAndSteps ( hgcons );
        analyzeTrackPosition ( simhits_globalpoints, generalTracks, b_field, *handle_AnalyticalPropagator );
        analyzeTICLTrackster ( selected_calotruths, *handle_TracksterHFNoseEM.product(), "EMn" );
        analyzeTICLTrackster ( selected_calotruths, *handle_TracksterHFNoseTrkEM.product(), "TrkEMn" );
        analyzeTICLTrackster ( selected_calotruths, *handle_TracksterHFNoseHAD.product(), "HADn" );
        analyzeTICLTrackster ( selected_calotruths, *handle_TracksterHFNoseMIP.product(), "MIPn" );
        comparePropagators ( generalTracks, b_field, *handle_RKPropagator, *handle_SteppingHelixPropagator );
        //analyzeTrackWithTrackster ( selected_calotruths, simhits_globalpoints, *handle_Tracks.product(), *handle_TracksterHFNoseTrkEM.product(), b_field, *handle_AnalyticalPropagator );
    }
    else std::cout << "Handle(s) invalid! Please investigate." << std::endl;
}


std::vector<math::XYZTLorentzVectorF> PropagatorAnalyzer::getTruthP4 ( const std::vector<CaloParticle> & caloTruthParticles )
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


std::vector<GlobalPoint> PropagatorAnalyzer::getSimHitGlobalPoint ( const std::vector<CaloParticle> & caloTruthParticles, const HGCalGeometry * geom )
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
            DetId    simhit_detid;
            Int_t    num_sim_clusters_first_layer = 0; // # sim clusters in layer 1
            
            for ( auto const& it_sc: ct_sim_clusters )
            {
                const SimCluster sc = *it_sc;
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
                
                if ( simcluster_seed_layer == 1 ) num_sim_clusters_first_layer++;                
            }
            
            //if ( num_sim_clusters_first_layer > 1 ) continue;
            
            try
            {
                const GlobalPoint & hit_globalPosition = geom->getPosition(simhit_detid);
                if ( abs(hit_globalPosition.eta()) > select_EtaLow_
                        && abs(hit_globalPosition.eta()) < select_EtaHigh_ )
                {
                    container.push_back ( hit_globalPosition );
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


void PropagatorAnalyzer::fillTruthHistograms ( const std::vector<math::XYZTLorentzVectorF> & caloTruthP4, const std::vector<GlobalPoint> & simHitsGlobalPoints )
{
    for ( auto const& truth: caloTruthP4 )
    {
        // Fill truth histograms
        TH1Container_["truthE"]->Fill( truth.energy() );
        TH1Container_["truthEta"]->Fill( truth.eta() );   
    }
    for ( auto const& hit_point: simHitsGlobalPoints )
    {
        TH1Container_["first_layer_SimHit_z"]->Fill( hit_point.z() );
        TH1Container_["first_layer_SimHit_eta"]->Fill( hit_point.eta() );
    }
}


void PropagatorAnalyzer::buildFirstLayersAndSteps ( const HGCalDDDConstants* hgcons )
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
    
    // TODO: buildsteps
    // for ()
}


void PropagatorAnalyzer::analyzeTrackPosition ( const std::vector<math::XYZTLorentzVectorF> & caloTruthP4, const std::vector<reco::Track> & generalTracks, const MagneticField * BField, const Propagator & prop )
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
            if ( !cutTk_(track) ) continue;
        
            FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(track, BField);
            int iSide = int(track.eta() > 0);
            TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
            
            if (!tsos.isValid())
                std::cout << "  Message in " << __func__ << ":  tsos invalid!" << std::endl;
            else
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


void PropagatorAnalyzer::analyzeTrackPosition ( const std::vector<GlobalPoint> & simHitsGlobalPoints, const std::vector<reco::Track> & generalTracks, const MagneticField * BField, const Propagator & prop )
{ // with simhits

    for ( auto const& hit_point: simHitsGlobalPoints )
    {
        Float_t truth_eta = hit_point.eta();
        Float_t truth_phi = hit_point.phi();
        
        Float_t min_R = 999.;
        Float_t truth_prop_dEta = 999.;
        Float_t truth_prop_dPhi = 999.;
        
        Float_t track_prop_dEta = 999.;
        Float_t track_prop_dPhi = 999.;
        
        for ( auto const& track: generalTracks )
        {
            if ( !cutTk_(track) ) continue;
        
            FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(track, BField);
            int iSide = int(track.eta() > 0);
            TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
            
            if (!tsos.isValid()) 
                std::cout << "  Message in " << __func__ << ":  tsos invalid!" << std::endl;
            else
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


void PropagatorAnalyzer::analyzeTICLTrackster ( const std::vector<math::XYZTLorentzVectorF> & caloTruthP4, const std::vector<ticl::Trackster> & tracksters, std::string tag )
{
    for ( auto const& trs: tracksters )
    {
        Float_t trackster_raw_energy = trs.raw_energy();
        Float_t trackster_eta = trs.barycenter().eta();
        Float_t trackster_phi = trs.barycenter().phi();
                
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


void PropagatorAnalyzer::analyzeTrackWithTrackster ( const std::vector<math::XYZTLorentzVectorF> & caloTruthP4, const std::vector<GlobalPoint> & simHitsGlobalPoints, const std::vector<reco::Track> & generalTracks, const std::vector<ticl::Trackster> & tracksters, const MagneticField * BField, const Propagator & prop )
{ // apply trackster raw energy cut

    std::vector<ticl::Trackster> selected_tracksters;
    
    // select tracksters
    for ( auto const& trs: tracksters )
    {
        //Float_t trackster_raw_energy = trs.raw_energy();
        Float_t trackster_eta = trs.barycenter().eta();
        Float_t trackster_phi = trs.barycenter().phi();
                
        for ( auto const& truth: caloTruthP4 )
        {
            Float_t dR = reco::deltaR ( trackster_eta, trackster_phi, truth.eta(), truth.phi() );

            if ( dR < (truth_matching_deltaR_) )
            {
                selected_tracksters.push_back(trs);
            }
        }
    }
    
    
    // select tracksters on same zside as simhit
    // + select trackster closest to simhit
    // TODO: apply cut --> abs(simhit_phi - prop_phi) < certain num
    for ( auto const& hit_point: simHitsGlobalPoints )
    {
        Float_t truth_eta = hit_point.eta();
        Float_t truth_phi = hit_point.phi();
        
        // First, get the distance from simhit to propagation
        Float_t truth_prop_dEta = 999.;
        Float_t truth_prop_dPhi = 999.;
        Float_t truth_prop_minR = 999.;
        
        for ( auto const& track: generalTracks )
        {
            if ( !cutTk_(track) ) continue;
        
            FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(track, BField);
            int iSide = int(track.eta() > 0);
            TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
            
            if (!tsos.isValid())
                std::cout << "  Message in " << __func__ << ":  tsos invalid!" << std::endl;
            else
            {            
                Float_t tsos_eta = tsos.globalPosition().eta();
                Float_t tsos_phi = tsos.globalPosition().phi();
                
                Float_t dR = reco::deltaR ( tsos_eta, tsos_phi, truth_eta, truth_phi );

                if ( dR < truth_matching_deltaR_ && dR < truth_prop_minR ) // within cone of truth
                {
                    truth_prop_minR = dR;
                    truth_prop_dEta = truth_eta - tsos_eta;
                    truth_prop_dPhi = truth_phi - tsos_phi;
                }
            }
            
        }
        
        if ( truth_prop_minR > truth_matching_deltaR_ ) continue;
//        // temporary! run check
        TH1Container_["dEta_Truth_Prop"]->Fill(truth_prop_dEta);
        TH1Container_["dPhi_Truth_Prop"]->Fill(truth_prop_dPhi);
        TH1Container_["dR_Truth_Prop"]->Fill(truth_prop_minR);
        

        // next, select trackster close/closest to simhit
        // and use the truth_prop_* variables above as plots
        Float_t minR_trackster_raw_energy = 0;
        for ( auto const& trs: selected_tracksters )
        {
            Float_t trackster_eta = trs.barycenter().eta();
            Float_t trackster_phi = trs.barycenter().phi();
            Float_t trackster_raw_energy = trs.raw_energy();
            
            Float_t dR = reco::deltaR ( trackster_eta, trackster_phi, hit_point.eta(), hit_point.phi() );
            Float_t truth_trackster_minR = 999.;
            
            // select same zside trackster
            if ( truth_eta * trackster_eta > 0 )
            {
                // add to TH2 histogram
                TH2Container_["dEta_Truth_Prop_EDist_zSide_Trackster"]->Fill(truth_prop_dEta, trackster_raw_energy);
                TH2Container_["dPhi_Truth_Prop_EDist_zSide_Trackster"]->Fill(truth_prop_dPhi, trackster_raw_energy);
            }
            // select closest trackster
            if ( dR < truth_trackster_minR )
            {
                truth_trackster_minR = dR;
                minR_trackster_raw_energy = trackster_raw_energy;
                // add to TH2 histogram (outside loop; see below)
            }
        }
        if ( minR_trackster_raw_energy > 0.1 )
        {
            TH2Container_["dEta_Truth_Prop_EDist_minR_Trackster"]->Fill(truth_prop_dEta, minR_trackster_raw_energy);
            TH2Container_["dPhi_Truth_Prop_EDist_minR_Trackster"]->Fill(truth_prop_dPhi, minR_trackster_raw_energy);
        }
    }
}


void PropagatorAnalyzer::comparePropagators ( const std::vector<reco::Track> & generalTracks, const MagneticField* BField, const Propagator & RungeKutta, const Propagator & SteppingHelix )
{
    for ( auto const& track: generalTracks )
    {
        if ( !cutTk_(track) ) continue;
    
        FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(track, BField);
        int iSide = int(track.eta() > 0);
        
        TrajectoryStateOnSurface tsos_RK = RungeKutta.propagate(fts, firstDisk_[iSide]->surface());
        TrajectoryStateOnSurface tsos_SH = SteppingHelix.propagate(fts, firstDisk_[iSide]->surface());
        
        if (!tsos_RK.isValid())
            std::cout << "  Message in " << __func__ << ": RK tsos invalid!" << std::endl;
        if (!tsos_SH.isValid())
            std::cout << "  Message in " << __func__ << ": SteppingHelix tsos invalid!" << std::endl;
        if (tsos_RK.isValid() && tsos_SH.isValid())
        {
            Float_t tsos_RK_eta = tsos_RK.globalPosition().eta();
            Float_t tsos_RK_phi = tsos_RK.globalPosition().phi();
            Float_t tsos_SH_eta = tsos_SH.globalPosition().eta();
            Float_t tsos_SH_phi = tsos_SH.globalPosition().phi();
            
            Float_t RK_SH_dEta = tsos_RK_eta - tsos_SH_eta;
            Float_t RK_SH_dPhi = tsos_RK_phi - tsos_SH_phi;
            
            TH1Container_["RK_SH_dEta"]->Fill(RK_SH_dEta);
            TH1Container_["RK_SH_dPhi"]->Fill(RK_SH_dPhi);
        }
    }
}


void PropagatorAnalyzer::beginJob ()
{
    edm::Service<TFileService> fs;
    
    // ---- CaloTruth ------
    TH1Container_["truthE"] = fs->make<TH1F>("truthE", "Truth Energy Distribution", 100, 0, 500);
    TH1Container_["truthEta"] = fs->make<TH1F>("truthEta", "Truth Eta Distribution", 40, (select_EtaLow_ + select_EtaHigh_)/2 - 0.5, (select_EtaLow_ + select_EtaHigh_)/2 + 0.5);
    
    // ---- SimHits ------
    TH1Container_["first_layer_SimHit_z"] = fs->make<TH1F>("first_layer_SimHit_z", "z of highest energy simhit in first layer", 150, 0, 1500);
    TH1Container_["first_layer_SimHit_eta"] = fs->make<TH1F>("first_layer_SimHit_eta", "eta of highest energy simhit in first layer", 150, 3.0, 4.5);
    
    // ---- Propagator & SimHit ------
    TH1Container_["dEta_Truth_Prop"] = fs->make<TH1F>("dEta_Truth_Prop", "#eta difference between reco and sim at HGCal interface", 600, -0.3, 0.3);
    TH1Container_["dPhi_Truth_Prop"] = fs->make<TH1F>("dPhi_Truth_Prop", "#phi difference between reco and sim at HGCal interface", 600, -0.3, 0.3);
    TH1Container_["dR_Truth_Prop"] = fs->make<TH1F>("dR_Truth_Prop", "#Delta R between reco and sim at HGCal interface", 300, 0.0, 0.3);
    
    TH1Container_["dEta_Truth_Prop"]->GetXaxis()->SetTitle("#eta_{simhit} - #eta_{prop}");
    TH1Container_["dPhi_Truth_Prop"]->GetXaxis()->SetTitle("phi_{simhit} - #phi_{prop}");
    TH1Container_["dR_Truth_Prop"]->GetXaxis()->SetTitle("#Delta{R}_{simhit, prop}");
    
    // ---- Propagator ------
    TH1Container_["RK_SH_dEta"] = fs->make<TH1F>("RK_SH_dEta", "#eta difference between RK and SteppingHelix tsos", 600, -0.3, 0.3);
    TH1Container_["RK_SH_dPhi"] = fs->make<TH1F>("RK_SH_dPhi", "#phi difference between RK and SteppingHelix tsos", 600, -0.3, 0.3);
    
    TH1Container_["RK_SH_dEta"]->GetXaxis()->SetTitle("#eta_{RK} - #eta_{SteppingHelix}");
    TH1Container_["RK_SH_dPhi"]->GetXaxis()->SetTitle("#phi_{RK} - #phi_{SteppingHelix}");
    
    // ---- Trackster, Propagator, & SimHit ------
    TH2Container_["dEta_Truth_Prop_EDist_zSide_Trackster"] = fs->make<TH2F>("dEta_Truth_Prop_EDist_zSide_Trackster", "simhit #Delta#eta & same zSide trackster raw energy", 600, -0.3, 0.3, 100, 0, 500);
    TH2Container_["dPhi_Truth_Prop_EDist_zSide_Trackster"] = fs->make<TH2F>("dPhi_Truth_Prop_EDist_zSide_Trackster", "simhit #Delta#phi & same zSide trackster raw energy", 600, -0.3, 0.3, 100, 0, 500);
    
    TH2Container_["dEta_Truth_Prop_EDist_minR_Trackster"] = fs->make<TH2F>("dEta_Truth_Prop_EDist_minR_Trackster", "simhit #Delta#eta & closest trackster raw energy", 600, -0.3, 0.3, 100, 0, 500);
    TH2Container_["dPhi_Truth_Prop_EDist_minR_Trackster"] = fs->make<TH2F>("dPhi_Truth_Prop_EDist_minR_Trackster", "simhit #Delta#phi & closest trackster raw energy", 600, -0.3, 0.3, 100, 0, 500);

    TH2Container_["dEta_Truth_Prop_EDist_zSide_Trackster"]->GetXaxis()->SetTitle("#eta_{simhit} - #eta_{prop}");
    TH2Container_["dEta_Truth_Prop_EDist_zSide_Trackster"]->GetYaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    TH2Container_["dPhi_Truth_Prop_EDist_zSide_Trackster"]->GetXaxis()->SetTitle("#phi_{simhit} - #phi_{prop}");
    TH2Container_["dPhi_Truth_Prop_EDist_zSide_Trackster"]->GetYaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    
    TH2Container_["dEta_Truth_Prop_EDist_zSide_Trackster"]->GetXaxis()->SetTitle("#eta_{simhit} - #eta_{prop}");
    TH2Container_["dEta_Truth_Prop_EDist_zSide_Trackster"]->GetYaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    TH2Container_["dPhi_Truth_Prop_EDist_zSide_Trackster"]->GetXaxis()->SetTitle("#phi_{simhit} - #phi_{prop}");
    TH2Container_["dPhi_Truth_Prop_EDist_zSide_Trackster"]->GetYaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");

    // ---- Tracks -------
    TH1Container_["dEta_Track_Prop"] = fs->make<TH1F>("dEta_Track_Prop", "#eta difference between #eta_{prop} and #eta_{track, outer}", 300, 0., 0.3);
    TH1Container_["dPhi_Track_Prop"] = fs->make<TH1F>("dPhi_Track_Prop", "#phi difference between #phi_{prop} and #phi_{track, outer}", 300, 0., 0.3);
    
    TH1Container_["dEta_Track_Prop"]->GetXaxis()->SetTitle("#eta_{track} - #eta_{prop}");
    TH1Container_["dPhi_Track_Prop"]->GetXaxis()->SetTitle("phi_{track} - #phi_{prop}");  
    
    // ---- Tracksters ------
    TH1Container_["RawEDist_tracksterHFNoseEM"] = fs->make<TH1F>("EDist_tracksterHFNoseEM", "TracksterHFNoseEM Raw Energy Distribution", 100, 0, 500);
    TH1Container_["RawEScale_tracksterHFNoseEM"] = fs->make<TH1F>("EScale_tracksterHFNoseEM", "TracksterHFNoseEM Raw Energy Scale", 100, 0, 500);
    TH1Container_["DeltaR_tracksterHFNoseEM"] = fs->make<TH1F>("DeltaR_tracksterHFNoseEM", "TracksterHFNoseEM #Delta R", 20, 0, 0.5);
    
    TH1Container_["RawEDist_tracksterHFNoseEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    TH1Container_["RawEScale_tracksterHFNoseEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    TH1Container_["DeltaR_tracksterHFNoseEM"]->GetXaxis()->SetTitle("|R_{trackster - caloParticle}|");
    
    //
    TH1Container_["RawEDist_tracksterHFNoseTrkEM"] = fs->make<TH1F>("EDist_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM Raw Energy Distribution", 100, 0, 500);
    TH1Container_["RawEScale_tracksterHFNoseTrkEM"] = fs->make<TH1F>("EScale_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM Raw Energy Scale", 100, 0, 500);
    TH1Container_["DeltaR_tracksterHFNoseTrkEM"] = fs->make<TH1F>("DeltaR_tracksterHFNoseTrkEM", "TracksterHFNoseTrkEM #Delta R", 20, 0, 0.5);
    
    TH1Container_["RawEDist_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    TH1Container_["RawEScale_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    TH1Container_["DeltaR_tracksterHFNoseTrkEM"]->GetXaxis()->SetTitle("|R_{trackster} - R_{caloParticle}|");
    
    //
    TH1Container_["RawEDist_tracksterEM"] = fs->make<TH1F>("EDist_tracksterEM", "TracksterEM Raw Energy Distribution", 500, 0, 500);
    TH1Container_["RawEScale_tracksterEM"] = fs->make<TH1F>("EScale_tracksterEM", "TracksterEM Raw Energy Scale", 100, 0, 500);
    TH1Container_["DeltaR_tracksterEM"] = fs->make<TH1F>("DeltaR_tracksterEM", "TracksterEM #Delta R", 20, 0, 0.5);
    
    TH1Container_["RawEDist_tracksterEM"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    TH1Container_["RawEScale_tracksterEM"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    TH1Container_["DeltaR_tracksterEM"]->GetXaxis()->SetTitle("|R_{trackster} - R_{caloParticle}|");
    
    //
    TH1Container_["RawEDist_tracksterHFNoseHAD"] = fs->make<TH1F>("EDist_tracksterHFNoseHAD", "TracksterHFNoseHAD Raw Energy Distribution", 500, 0, 500);
    TH1Container_["RawEScale_tracksterHFNoseHAD"] = fs->make<TH1F>("EScale_tracksterHFNoseHAD", "TracksterHFNoseHAD Raw Energy Scale", 100, 0, 500);
    TH1Container_["DeltaR_tracksterHFNoseHAD"] = fs->make<TH1F>("DeltaR_tracksterHFNoseHAD", "TracksterHFNoseHAD #Delta R", 20, 0, 0.5);
    
    TH1Container_["RawEDist_tracksterHFNoseHAD"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    TH1Container_["RawEScale_tracksterHFNoseHAD"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    TH1Container_["DeltaR_tracksterHFNoseHAD"]->GetXaxis()->SetTitle("|R_{trackster} - R_{caloParticle}|");
    
    //
    TH1Container_["RawEDist_tracksterHFNoseMIP"] = fs->make<TH1F>("EDist_tracksterHFNoseMIP", "TracksterHFNoseMIP Raw Energy Distribution", 500, 0, 500);
    TH1Container_["RawEScale_tracksterHFNoseMIP"] = fs->make<TH1F>("EScale_tracksterHFNoseMIP", "TracksterHFNoseMIP Raw Energy Scale", 100, 0, 500);
    TH1Container_["DeltaR_tracksterHFNoseMIP"] = fs->make<TH1F>("DeltaR_tracksterHFNoseMIP", "TracksterHFNoseMIP #Delta R", 20, 0, 0.5);
    
    TH1Container_["RawEDist_tracksterHFNoseMIP"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    TH1Container_["RawEScale_tracksterHFNoseMIP"]->GetXaxis()->SetTitle("E_{caloParticle} - E_{trackster} [GeV/c^{2}]");
    TH1Container_["DeltaR_tracksterHFNoseMIP"]->GetXaxis()->SetTitle("|R_{trackster} - R_{caloParticle}|");

}


void PropagatorAnalyzer::endJob ()
{
    // Job ends
    std::cout << "The job, PropagatorAnalyzer, has ended. Thank you for your patience." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (PropagatorAnalyzer);

