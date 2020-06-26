#include "EMShowerStudies.h"

#include <iostream>
#include <cmath> // Switch to TMath.h if you need more physics-related functions
#include "DataFormats/Math/interface/deltaR.h"

#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Physics Objects
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

// HGCNose
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

// CMS Coordinate System
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

// Detector Geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"


EMShowerStudies::EMShowerStudies ( const edm::ParameterSet& iConfig ) :

    TH1_Container_ (),
    TH2_Container_ (),
    
    // (tag name, default value (label, instance, process) -- CHECK SPELLING!!!!!!!    
    // Truth objects
    // tag_CaloParticle_MergedCaloTruth_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_MergedCaloTruth", edm::InputTag ("mix", "MergedCaloTruth") ) ), 
    tag_GenParticle_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_GenParticle", edm::InputTag ("genParticles") ) ),
    
    // Sim level objects
    tag_g4SimHits_HFNoseHits_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_g4SimHits_HFNoseHits", edm::InputTag ( "g4SimHits", "HFNoseHits", "SIM" ) ) ),
    
    // Reco level objects
    tag_HGCHFNoseRecHits_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_HGCHFNoseRecHits", edm::InputTag ("HGCalRecHit", "HGCHFNoseRecHits") ) ),
    tag_HGCalLayerClustersHFNose_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_HGCalLayerClusterHFNose", edm::InputTag ("hgcalLayerClustersHFNose") ) ),
    
    // Pre-selection parameters
    select_PID_ ( 22 ),
    select_EtaLow_ ( 3.49 ),
    select_EtaHigh_ ( 3.51 ),
    max_iter_R_ ( 200 ),
    steps_iter_R_ ( 40 )
    
{
    // consumes: frequent request of additional data | mayConsume: infrequent
    // token_CaloParticle_MergedCaloTruth_ = mayConsume<std::vector<CaloParticle>> ( tag_CaloParticle_MergedCaloTruth_ );
    token_GenParticle_ = consumes<reco::GenParticleCollection> ( tag_GenParticle_ );
    
    token_g4SimHits_HFNoseHits_ = consumes<std::vector<PCaloHit>> ( tag_g4SimHits_HFNoseHits_ );
    
    token_HGCRecHits_ = consumes<HGCRecHitCollection> ( tag_HGCHFNoseRecHits_ );
    token_HGCalLayerClustersHFNose_ = consumes<std::vector<reco::CaloCluster>> ( tag_HGCalLayerClustersHFNose_ );
}


EMShowerStudies::~EMShowerStudies ()
{
    // Deconstructor
}


void EMShowerStudies::analyze ( const edm::Event & iEvent, const edm::EventSetup & iSetup )
{
    // Get HGCRecHits
    edm::Handle<HGCRecHitCollection> handle_HGCRecHits;
    iEvent.getByToken ( token_HGCRecHits_, handle_HGCRecHits );

    // Get CaloClusters
    edm::Handle<std::vector<reco::CaloCluster>> handle_HGCalLayerClustersHFNose;
    iEvent.getByToken ( token_HGCalLayerClustersHFNose_, handle_HGCalLayerClustersHFNose );
    
    // Get g4SimHits:HFNoseHits
    edm::Handle<std::vector<PCaloHit>> handle_g4SimHits_HFNoseHits;
    iEvent.getByToken ( token_g4SimHits_HFNoseHits_, handle_g4SimHits_HFNoseHits );

    // Get CaloParticles:CaloTruth
    // edm::Handle<std::vector<CaloParticle>> handle_CaloParticle_MergedCaloTruth;
    // iEvent.getByToken ( token_CaloParticle_MergedCaloTruth_, handle_CaloParticle_MergedCaloTruth );
    
    // Get GenParticles
    edm::Handle<reco::GenParticleCollection> handle_GenParticle;
    iEvent.getByToken ( token_GenParticle_, handle_GenParticle );

    // Get HGCalGeometry (from EventSetup -- that is why no tag or token exist in member data)
    edm::ESHandle<HGCalGeometry> handle_HGCalGeometry;
    iSetup.get<IdealGeometryRecord>().get( "HGCalHFNoseSensitive", handle_HGCalGeometry );
    // const HGCalGeometry* geom_HGCal = handle_HGCalGeometry.product();  
    
    
    // Find position and energy of the maximum E deposit HGCRecHit
    // const FrontBackEtaPhiE_perLayer MaxE_RecHits = find_EtaPhiE_Reference_perLayer ( *handle_HGCRecHits.product(), handle_HGCalGeometry.product() ); // Why not pointer?
    const FrontBackEtaPhiE_perLayer MaxE_CaloClusters = find_EtaPhiE_Reference_perLayer ( *handle_GenParticle.product(), handle_HGCalGeometry.product() );
    
    // Plot position of above
    plot_maxEtaPhi ( MaxE_CaloClusters );
    
    // Calculate sum of E deposit of SimHits in each layer 
    const FrontBackEtaPhiE_perLayer SumE_Clusters = get_SumEDeposit_perLayer ( *handle_HGCalLayerClustersHFNose.product(), handle_HGCalGeometry.product() );
    
    // Validate sum
//    std::array<bool, 2> bool_check_sums = check_SumEDeposit_allLayers ( *handle_GenParticle.product(), SumE_Clusters );
//    std::cout << std::boolalpha
//              << "Sum validation results (5%): " << bool_check_sums[0] << ", "
//              << bool_check_sums[1] << std::endl;
    
    // Plot sum of SimHits
    plot_sum_TotalE_perLayer ( SumE_Clusters );

    // Construct radii with HGCRecHits
    iterative_R_search ( *handle_HGCalLayerClustersHFNose.product(), handle_HGCalGeometry.product(), MaxE_CaloClusters, SumE_Clusters );
    
}


EMShowerStudies::FrontBackEtaPhiE_perLayer EMShowerStudies::find_EtaPhiE_Reference_perLayer ( const reco::GenParticleCollection & GenParticles, const HGCalGeometry * geom )
{ // Find (eta, phi, energy) of CaloClusters with maximum energy deposit in the calorimeter for each layer. (3 of 3 overloaded methods)
    
    // Container to return
    FrontBackEtaPhiE_perLayer EtaPhiE_MaximumEDeposit;
    
    bool check_front = false;
    bool check_back = false;
    for ( const auto& gen: GenParticles )
    {
        if ( gen.pdgId() == select_PID_
                && abs(gen.eta()) > select_EtaLow_
                && abs(gen.eta()) < select_EtaHigh_  )
        {
            if ( gen.eta() > 0 ) // +z scenario
            {
                for ( int layer = 1; layer <= 8; layer++ )
                {
                    EtaPhiE_MaximumEDeposit[layer-1][0] = gen.p4().x();
                    EtaPhiE_MaximumEDeposit[layer-1][1] = gen.p4().y();
                    EtaPhiE_MaximumEDeposit[layer-1][2] = gen.energy();
                    check_front = true;
                }
            }
            else // -z scenario
            {
                for ( int layer = 1; layer <= 8; layer++ )
                {
                    EtaPhiE_MaximumEDeposit[layer-1][3] = gen.p4().x();
                    EtaPhiE_MaximumEDeposit[layer-1][4] = gen.p4().y();
                    EtaPhiE_MaximumEDeposit[layer-1][5] = gen.energy();
                    check_back = true;
                }        
            }
        }
        else std::cout << "GenParticle not found!" << std::endl;
    }
    if (!check_front && !check_back) std::cout << "GenParticle not found!" << std::endl;
    else if (!check_front || !check_back) std::cout << "GenParticle did not come in pairs!" << std::endl;
    
    return EtaPhiE_MaximumEDeposit;
}


void EMShowerStudies::plot_maxEtaPhi ( const FrontBackEtaPhiE_perLayer max_EtaPhiE_layer )
{ // Plot the (eta, phi) positions of maximum HGCRecHit on TH2F

    TH2_Container_["maxEtaPhi_layer1"]->Fill( max_EtaPhiE_layer[0][0], max_EtaPhiE_layer[0][1] );
    TH2_Container_["maxEtaPhi_layer1"]->Fill( max_EtaPhiE_layer[0][3], max_EtaPhiE_layer[0][4] );
    
    TH2_Container_["maxEtaPhi_layer2"]->Fill( max_EtaPhiE_layer[1][0], max_EtaPhiE_layer[1][1] );
    TH2_Container_["maxEtaPhi_layer2"]->Fill( max_EtaPhiE_layer[1][3], max_EtaPhiE_layer[1][4] );
    
    TH2_Container_["maxEtaPhi_layer3"]->Fill( max_EtaPhiE_layer[2][0], max_EtaPhiE_layer[2][1] );
    TH2_Container_["maxEtaPhi_layer3"]->Fill( max_EtaPhiE_layer[2][3], max_EtaPhiE_layer[2][4] );
    
    TH2_Container_["maxEtaPhi_layer4"]->Fill( max_EtaPhiE_layer[3][0], max_EtaPhiE_layer[3][1] );
    TH2_Container_["maxEtaPhi_layer4"]->Fill( max_EtaPhiE_layer[3][3], max_EtaPhiE_layer[3][4] );
    
    TH2_Container_["maxEtaPhi_layer5"]->Fill( max_EtaPhiE_layer[4][0], max_EtaPhiE_layer[4][1] );
    TH2_Container_["maxEtaPhi_layer5"]->Fill( max_EtaPhiE_layer[4][3], max_EtaPhiE_layer[4][4] );
    
    TH2_Container_["maxEtaPhi_layer6"]->Fill( max_EtaPhiE_layer[5][0], max_EtaPhiE_layer[5][1] );
    TH2_Container_["maxEtaPhi_layer6"]->Fill( max_EtaPhiE_layer[5][3], max_EtaPhiE_layer[5][4] );
    
    TH2_Container_["maxEtaPhi_layer7"]->Fill( max_EtaPhiE_layer[6][0], max_EtaPhiE_layer[6][1] );
    TH2_Container_["maxEtaPhi_layer7"]->Fill( max_EtaPhiE_layer[6][3], max_EtaPhiE_layer[6][4] );
    
    TH2_Container_["maxEtaPhi_layer8"]->Fill( max_EtaPhiE_layer[7][0], max_EtaPhiE_layer[7][1] );
    TH2_Container_["maxEtaPhi_layer8"]->Fill( max_EtaPhiE_layer[7][3], max_EtaPhiE_layer[7][4] );

}


EMShowerStudies::FrontBackEtaPhiE_perLayer EMShowerStudies::get_SumEDeposit_perLayer ( const std::vector<reco::CaloCluster> & Clusters, const HGCalGeometry * geom )
{ // Sum all the energies of the CaloClusters in each layer (3 of 3 overloaded methods)
    
    // Container to return
    FrontBackEtaPhiE_perLayer SumEDeposit_perLayer;
    SumEDeposit_perLayer.fill( std::array<Float_t, 6> { { 0., 0., 0., 0., 0., 0. } } );
    
    for ( const auto& cl: Clusters )
    {
        // Get layer number of the cluster
        HFNoseDetId cl_DetId = HFNoseDetId( cl.hitsAndFractions().at(0).first );
        Int_t layer = cl_DetId.layer();
        
        // Troubleshoot: if index out of range, doublecheck HGCNose_NLayers_ value in constructor.
        
        if ( cl.eta() > 0 ) SumEDeposit_perLayer[layer-1][2] += cl.energy();
        if ( cl.eta() < 0 ) SumEDeposit_perLayer[layer-1][5] += cl.energy();
    }
    
    Float_t sum_EDeposit_front = 0;
    for ( int i = 0; i < HGCNose_NLayers_; i++ )
    {
        sum_EDeposit_front += SumEDeposit_perLayer[i][2];
    }
    
    return SumEDeposit_perLayer;
}


void EMShowerStudies::plot_sum_TotalE_perLayer ( const FrontBackEtaPhiE_perLayer SumEDeposit_perLayer )
{ // Plot sum of jet energy deposits per layer on TH1F histograms

    TH1_Container_["Total_EDeposit_layer1"]->Fill( SumEDeposit_perLayer[0][2] );
    TH1_Container_["Total_EDeposit_layer1"]->Fill( SumEDeposit_perLayer[0][5] );
    
    TH1_Container_["Total_EDeposit_layer2"]->Fill( SumEDeposit_perLayer[1][2] );
    TH1_Container_["Total_EDeposit_layer2"]->Fill( SumEDeposit_perLayer[1][5] );
    
    TH1_Container_["Total_EDeposit_layer3"]->Fill( SumEDeposit_perLayer[2][2] );
    TH1_Container_["Total_EDeposit_layer3"]->Fill( SumEDeposit_perLayer[2][5] );
    
    TH1_Container_["Total_EDeposit_layer4"]->Fill( SumEDeposit_perLayer[3][2] );
    TH1_Container_["Total_EDeposit_layer4"]->Fill( SumEDeposit_perLayer[3][5] );
    
    TH1_Container_["Total_EDeposit_layer5"]->Fill( SumEDeposit_perLayer[4][2] );
    TH1_Container_["Total_EDeposit_layer5"]->Fill( SumEDeposit_perLayer[4][5] );
    
    TH1_Container_["Total_EDeposit_layer6"]->Fill( SumEDeposit_perLayer[5][2] );
    TH1_Container_["Total_EDeposit_layer6"]->Fill( SumEDeposit_perLayer[5][5] );
    
    TH1_Container_["Total_EDeposit_layer7"]->Fill( SumEDeposit_perLayer[6][2] );
    TH1_Container_["Total_EDeposit_layer7"]->Fill( SumEDeposit_perLayer[6][5] );
    
    TH1_Container_["Total_EDeposit_layer8"]->Fill( SumEDeposit_perLayer[7][2] );
    TH1_Container_["Total_EDeposit_layer8"]->Fill( SumEDeposit_perLayer[7][5] );
    
}


std::array<bool, 2> EMShowerStudies::check_SumEDeposit_allLayers ( const reco::GenParticleCollection & GenParticles, const FrontBackEtaPhiE_perLayer SumEDeposit_perLayer, const Float_t fractional_error )
{ // Check if the sum of sim hits in each layer really add up to the GenParticle energy

    std::array<bool, 2> return_bool;
    return_bool.fill(false);
    
    Float_t sum_EDeposit_front = 0;
    for ( int i = 0; i < HGCNose_NLayers_; i++ )
    {
        sum_EDeposit_front += SumEDeposit_perLayer[i][2];
    }
                                 
//    Float_t sum_EDeposit_back = lambda_sum_E(5);
    Float_t sum_EDeposit_back = 0;
    for ( int i = 0; i < HGCNose_NLayers_; i++ )
    {
        sum_EDeposit_back += SumEDeposit_perLayer[i][5];
    }

    for ( const auto& gen: GenParticles )
    {
        if ( gen.pdgId() == select_PID_ )
        { // Select by pdgId first
            if ( gen.eta() > select_EtaLow_ && gen.eta() < select_EtaHigh_ )
            { // Select by eta in the +z direction
                if ( abs(sum_EDeposit_front - gen.energy()) / gen.energy() < fractional_error ) return_bool[0] = true;
            }
            else if ( gen.eta() > -select_EtaHigh_ && gen.eta() < -select_EtaLow_ )
            { // Select by eta in the -z direction
                if ( abs(sum_EDeposit_back - gen.energy()) / gen.energy() < fractional_error ) return_bool[1] = true;
            }
        }
    }
    
    return return_bool;
}


Float_t EMShowerStudies::getContainmentR ( const std::vector<Float_t> iter_R, const std::vector<Float_t> R_frac_EDeposit, const Float_t FracContainment )
{ // Return R at which energy contained is a specific fraction of the total E deposit.

    // Check validity of vector
    bool nan_exists = false;
    bool all_zeros = true; // start with true and get false if any one is false
    bool all_ones = true;
    
    for ( auto const& it: R_frac_EDeposit )
    {
        if ( isnan(it) == true ) nan_exists = true;
        else if ( it != 0 ) all_zeros = false;
        else if ( it != 1. ) all_ones = false;
    }
    
    // Decide what to return
    if ( nan_exists || all_zeros || all_ones ) return 0;
    else
    { // Standard procedure when all things valid
        TGraph gr = TGraph ( steps_iter_R_, &iter_R[0], &R_frac_EDeposit[0] );
        TF1 fit = TF1 ( "R_containment_fit", "pol12", 0, max_iter_R_ );
        gr.Fit ( "R_containment_fit", "QS" );
        Float_t result = fit.GetX ( FracContainment, 0, 200, 1.E-8, 500 );
        return result;
    }
}


EMShowerStudies::FrontBackEtaPhiE_perLayer EMShowerStudies::getContainedEnergy ( const std::vector<reco::CaloCluster> & Clusters, const HGCalGeometry * geom, const FrontBackEtaPhiE_perLayer maxE_centers, const FrontBackEtaPhiE_perLayer frac_R )
{ // Get the total energy within the R given by the fraction of containment

    FrontBackEtaPhiE_perLayer EDeposit_layer;
    
    for ( const auto& cl: Clusters )
        {
            // Get layer number of the cluster
            HFNoseDetId cl_DetId = HFNoseDetId( cl.hitsAndFractions().at(0).first );
            Int_t layer = cl_DetId.layer();
            
            // Get GlobalPoint (position in global geometry)
            GlobalPoint cl_globalpoint = geom->getPosition( cl_DetId );
            
            if ( cl.eta() > 0 )
            { // +z position
                // GlobalPoint in mm, GenParticle.P4() in cm
                Float_t dx = ( cl_globalpoint.x() - 10*maxE_centers[layer-1][0] );
                Float_t dy = ( cl_globalpoint.y() - 10*maxE_centers[layer-1][1] );
                Float_t dR = sqrt ( dx*dx + dy*dy );
                // std::cout << "X: " << cl_globalpoint.x() << " Y: " << cl_globalpoint.y() << std::endl;
                // std::cout << dR << std::endl;
                if ( dR < frac_R[layer-1][0] ) EDeposit_layer[layer-1][2] += cl.energy();
            }
            else
            { // -z position
                Float_t dx = ( cl_globalpoint.x() - 10*maxE_centers[layer-1][3] );
                Float_t dy = ( cl_globalpoint.y() - 10*maxE_centers[layer-1][4] );
                Float_t dR = sqrt ( dx*dx + dy*dy );

                if ( dR < frac_R[layer-1][3] ) EDeposit_layer[layer-1][5] += cl.energy();
            }
        }
        
    return EDeposit_layer;

}


std::array<Int_t, 2> EMShowerStudies::getContainmentLayer ( const FrontBackEtaPhiE_perLayer EDeposit_layer, const FrontBackEtaPhiE_perLayer frac_R, Float_t FracContainment )
{ // Get layer which cumulatively contains x % of the total deposited energy

    // Get sum
    Float_t total_front = 0.;
    Float_t total_back = 0.;
    
    for ( int layer = 1; layer <= HGCNose_NLayers_; layer++ )
    {
        total_front += EDeposit_layer[layer-1][2];
        total_back += EDeposit_layer[layer-1][5];
    }
    
    // Get layer
    Int_t layer_front = 0;
    Int_t layer_back = 0;
    
    Float_t sum_front = 0.;
    Float_t sum_back = 0.;
    
    for ( int layer = 1; layer <= HGCNose_NLayers_; layer++ )
    {
        sum_front += EDeposit_layer[layer-1][2];
        layer_front++;
        if ( sum_front >= total_front * FracContainment ) break;
    }
    
    for ( int layer = 1; layer <= HGCNose_NLayers_; layer++ )
    {
        sum_back += EDeposit_layer[layer-1][5];
        layer_back++;
        if ( sum_back >= total_back * FracContainment ) break;
    }
    
    return std::array<Int_t, 2> { { layer_front, layer_back } };

}



void EMShowerStudies::iterative_R_search ( const std::vector<reco::CaloCluster> & Clusters, const HGCalGeometry * geom, const FrontBackEtaPhiE_perLayer maxE_centers, const FrontBackEtaPhiE_perLayer TotalE_perLayer )
{ // Iteratively expand dR from the reference and add all the CaloClusters' energies within it. Save to TH2F histograms. (2 of 2 overloaded methods.)

    // Containers for energy deposit within R
    Float_t increment = max_iter_R_/steps_iter_R_;
    std::array<std::vector<Float_t>, HGCNose_NLayers_> R_EDeposit_layer_front;
    std::array<std::vector<Float_t>, HGCNose_NLayers_> R_EDeposit_layer_back;
    
    // R values to iterate over
    std::vector<Float_t> iter_R ( steps_iter_R_ );
    
    // Fill container with zeros
    for ( int layer = 1; layer <= HGCNose_NLayers_; layer++ )
    {
        R_EDeposit_layer_front[layer-1].resize( steps_iter_R_, 0. );
        R_EDeposit_layer_back[layer-1].resize( steps_iter_R_, 0. );
        
        // std::cout << "X: " << maxE_centers[layer-1][0] << " Y: " << maxE_centers[layer-1][1] << std::endl;
    }

    for ( int step = 1; step <= steps_iter_R_; step++ )
    {    
        Float_t R = step * increment;
        iter_R[step-1] = R;
        
        for ( const auto& cl: Clusters )
        {
            // Get layer number of the cluster
            HFNoseDetId cl_DetId = HFNoseDetId( cl.hitsAndFractions().at(0).first );
            Int_t layer = cl_DetId.layer();
            
            // Get GlobalPoint (position in global geometry)
            GlobalPoint cl_globalpoint = geom->getPosition( cl_DetId );
            
            if ( cl.eta() > 0 )
            { // +z position
                // GlobalPoint in mm, GenParticle.P4() in cm
                Float_t dx = ( cl_globalpoint.x() - 10*maxE_centers[layer-1][0] );
                Float_t dy = ( cl_globalpoint.y() - 10*maxE_centers[layer-1][1] );
                Float_t dR = sqrt ( dx*dx + dy*dy );
                // std::cout << "X: " << cl_globalpoint.x() << " Y: " << cl_globalpoint.y() << std::endl;
                // std::cout << dR << std::endl;
                if ( dR < R ) R_EDeposit_layer_front[layer-1][step-1] += cl.energy();
            }
            else
            { // -z position
                Float_t dx = ( cl_globalpoint.x() - 10*maxE_centers[layer-1][3] );
                Float_t dy = ( cl_globalpoint.y() - 10*maxE_centers[layer-1][4] );
                Float_t dR = sqrt ( dx*dx + dy*dy );

                if ( dR < R ) R_EDeposit_layer_back[layer-1][step-1] += cl.energy();
            }
        }
    }
    
    // Container to hold R values (abuse of containers, really, but I'm lazy)
    FrontBackEtaPhiE_perLayer frac_R;
    
    for ( int layer = 1; layer <= HGCNose_NLayers_; layer++ )
    {
        std::vector<Float_t> R_frac_EDeposit_front;
        std::vector<Float_t> R_frac_EDeposit_back;
        
        for ( auto const& it: R_EDeposit_layer_front[layer-1] )
        {
            R_frac_EDeposit_front.emplace_back( it / TotalE_perLayer[layer-1][2] );
        }
        
        for ( auto const& it: R_EDeposit_layer_back[layer-1] )
        {
            R_frac_EDeposit_back.emplace_back( it / TotalE_perLayer[layer-1][5] );
        }
        
        frac_R[layer-1][0] = getContainmentR ( iter_R, R_frac_EDeposit_front, 0.9 );
        frac_R[layer-1][3] = getContainmentR ( iter_R, R_frac_EDeposit_back, 0.9 );

        // Save histograms
        std::string key = "Ninety_Percent_R_layer" + std::to_string(layer);
        
        if ( frac_R[layer-1][0] != 0 )
        {
            TH1_Container_[key.c_str()]->Fill( frac_R[layer-1][0] );
            TH1_Container_["Nevents_withJets_PerLayer"]->Fill( layer );
        }
        if ( frac_R[layer-1][3] != 0 )
        {
            TH1_Container_[key.c_str()]->Fill( frac_R[layer-1][3] );
            TH1_Container_["Nevents_withJets_PerLayer"]->Fill( layer );
        }
    }
    
    auto EDeposit_layer = getContainedEnergy ( Clusters, geom, maxE_centers, frac_R );
    
    for ( int layer = 1; layer <= HGCNose_NLayers_; layer++ )
    {
        TH2_Container_["Ninety_Percent_E"]->Fill( layer, EDeposit_layer[layer-1][2] );
        TH2_Container_["Ninety_Percent_E"]->Fill( layer, EDeposit_layer[layer-1][5] );
    }
    
    std::array<Int_t, 2> layer_front_back = getContainmentLayer ( EDeposit_layer, frac_R, 0.9);
    TH1_Container_["Ninety_Percent_Layer"]->Fill( layer_front_back[0] );
    TH1_Container_["Ninety_Percent_Layer"]->Fill( layer_front_back[1] );
}


void EMShowerStudies::beginJob ()
{
    edm::Service<TFileService> fs;
    
    // TH2 histograms about maximum E deposit location
    TH2_Container_["maxEtaPhi_layer1"] = fs->make<TH2F>("maxEtaPhi_layer1", "(#eta, #phi) of maximum energy deposit in layer1", 64, -4., 4., 52, -3.25, 3.25); // bin size is 0.125 on both axes
    TH2_Container_["maxEtaPhi_layer2"] = fs->make<TH2F>("maxEtaPhi_layer2", "(#eta, #phi) of maximum energy deposit in layer2", 64, -4., 4., 52, -3.25, 3.25);
    TH2_Container_["maxEtaPhi_layer3"] = fs->make<TH2F>("maxEtaPhi_layer3", "(#eta, #phi) of maximum energy deposit in layer3", 64, -4., 4., 52, -3.25, 3.25);
    TH2_Container_["maxEtaPhi_layer4"] = fs->make<TH2F>("maxEtaPhi_layer4", "(#eta, #phi) of maximum energy deposit in layer4", 64, -4., 4., 52, -3.25, 3.25);
    TH2_Container_["maxEtaPhi_layer5"] = fs->make<TH2F>("maxEtaPhi_layer5", "(#eta, #phi) of maximum energy deposit in layer5", 64, -4., 4., 52, -3.25, 3.25);
    TH2_Container_["maxEtaPhi_layer6"] = fs->make<TH2F>("maxEtaPhi_layer6", "(#eta, #phi) of maximum energy deposit in layer6", 64, -4., 4., 52, -3.25, 3.25);
    TH2_Container_["maxEtaPhi_layer7"] = fs->make<TH2F>("maxEtaPhi_layer7", "(#eta, #phi) of maximum energy deposit in layer7", 64, -4., 4., 52, -3.25, 3.25);
    TH2_Container_["maxEtaPhi_layer8"] = fs->make<TH2F>("maxEtaPhi_layer8", "(#eta, #phi) of maximum energy deposit in layer8", 64, -4., 4., 52, -3.25, 3.25);
    
    TH2_Container_["maxEtaPhi_layer1"]->GetXaxis()->SetTitle("#eta");
    TH2_Container_["maxEtaPhi_layer1"]->GetYaxis()->SetTitle("#phi");
    TH2_Container_["maxEtaPhi_layer2"]->GetXaxis()->SetTitle("#eta");
    TH2_Container_["maxEtaPhi_layer2"]->GetYaxis()->SetTitle("#phi");
    TH2_Container_["maxEtaPhi_layer3"]->GetXaxis()->SetTitle("#eta");
    TH2_Container_["maxEtaPhi_layer3"]->GetYaxis()->SetTitle("#phi");
    TH2_Container_["maxEtaPhi_layer4"]->GetXaxis()->SetTitle("#eta");
    TH2_Container_["maxEtaPhi_layer4"]->GetYaxis()->SetTitle("#phi");
    TH2_Container_["maxEtaPhi_layer5"]->GetXaxis()->SetTitle("#eta");
    TH2_Container_["maxEtaPhi_layer5"]->GetYaxis()->SetTitle("#phi");
    TH2_Container_["maxEtaPhi_layer6"]->GetXaxis()->SetTitle("#eta");
    TH2_Container_["maxEtaPhi_layer6"]->GetYaxis()->SetTitle("#phi");
    TH2_Container_["maxEtaPhi_layer7"]->GetXaxis()->SetTitle("#eta");
    TH2_Container_["maxEtaPhi_layer7"]->GetYaxis()->SetTitle("#phi");
    TH2_Container_["maxEtaPhi_layer8"]->GetXaxis()->SetTitle("#eta");
    TH2_Container_["maxEtaPhi_layer8"]->GetYaxis()->SetTitle("#phi");
    
    // TH2 histograms of FRACTION of energy included within cylinder radius
    TH2_Container_["R_frac_containment_layer1"] = fs->make<TH2F>("R_frac_containment_layer1", "E_{reco}/E_{gen} contained per #Delta R in layer1", steps_iter_R_, 0., max_iter_R_, 100, 0., 2.);
    TH2_Container_["R_frac_containment_layer2"] = fs->make<TH2F>("R_frac_containment_layer2", "E_{reco}/E_{gen} contained per #Delta R in layer2", steps_iter_R_, 0., max_iter_R_, 100, 0., 2.);
    TH2_Container_["R_frac_containment_layer3"] = fs->make<TH2F>("R_frac_containment_layer3", "E_{reco}/E_{gen} contained per #Delta Rin layer3", steps_iter_R_, 0., max_iter_R_, 100, 0., 2.);
    TH2_Container_["R_frac_containment_layer4"] = fs->make<TH2F>("R_frac_containment_layer4", "E_{reco}/E_{gen} contained per #Delta R in layer4", steps_iter_R_, 0., max_iter_R_, 100, 0., 2.);
    TH2_Container_["R_frac_containment_layer5"] = fs->make<TH2F>("R_frac_containment_layer5", "E_{reco}/E_{gen} contained per #Delta R in layer5", steps_iter_R_, 0., max_iter_R_, 100, 0., 2.);
    TH2_Container_["R_frac_containment_layer6"] = fs->make<TH2F>("R_frac_containment_layer6", "E_{reco}/E_{gen} contained per #Delta R in layer6", steps_iter_R_, 0., max_iter_R_, 100, 0., 2.);
    TH2_Container_["R_frac_containment_layer7"] = fs->make<TH2F>("R_frac_containment_layer7", "E_{reco}/E_{gen} contained per #Delta R in layer7", steps_iter_R_, 0., max_iter_R_, 100, 0., 2.);
    TH2_Container_["R_frac_containment_layer8"] = fs->make<TH2F>("R_frac_containment_layer8", "E_{reco}/E_{gen} contained per #Delta R in layer8", steps_iter_R_, 0., max_iter_R_, 100, 0., 2.);
    
    TH2_Container_["R_frac_containment_layer1"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_frac_containment_layer1"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_frac_containment_layer2"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_frac_containment_layer2"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_frac_containment_layer3"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_frac_containment_layer3"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_frac_containment_layer4"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_frac_containment_layer4"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_frac_containment_layer5"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_frac_containment_layer5"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_frac_containment_layer6"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_frac_containment_layer6"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_frac_containment_layer7"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_frac_containment_layer7"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_frac_containment_layer8"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_frac_containment_layer8"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    
    // TH2 histograms of energy included within cylinder radius
    TH2_Container_["R_containment_layer1"] = fs->make<TH2F>("R_containment_layer1", "E_{reco} contained per #Delta R in layer1", steps_iter_R_, 0., max_iter_R_, 100, 0., 100);
    TH2_Container_["R_containment_layer2"] = fs->make<TH2F>("R_containment_layer2", "E_{reco} contained per #Delta R in layer2", steps_iter_R_, 0., max_iter_R_, 100, 0., 100);
    TH2_Container_["R_containment_layer3"] = fs->make<TH2F>("R_containment_layer3", "E_{reco} contained per #Delta Rin layer3", steps_iter_R_, 0., max_iter_R_, 100, 0., 100);
    TH2_Container_["R_containment_layer4"] = fs->make<TH2F>("R_containment_layer4", "E_{reco} contained per #Delta R in layer4", steps_iter_R_, 0., max_iter_R_, 100, 0., 100);
    TH2_Container_["R_containment_layer5"] = fs->make<TH2F>("R_containment_layer5", "E_{reco} contained per #Delta R in layer5", steps_iter_R_, 0., max_iter_R_, 100, 0., 100);
    TH2_Container_["R_containment_layer6"] = fs->make<TH2F>("R_containment_layer6", "E_{reco} contained per #Delta R in layer6", steps_iter_R_, 0., max_iter_R_, 100, 0., 100);
    TH2_Container_["R_containment_layer7"] = fs->make<TH2F>("R_containment_layer7", "E_{reco} contained per #Delta R in layer7", steps_iter_R_, 0., max_iter_R_, 100, 0., 100);
    TH2_Container_["R_containment_layer8"] = fs->make<TH2F>("R_containment_layer8", "E_{reco} contained per #Delta R in layer8", steps_iter_R_, 0., max_iter_R_, 100, 0., 100);
    
    TH2_Container_["R_containment_layer1"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer1"]->GetYaxis()->SetTitle("E_{reco, #Delta R}");
    TH2_Container_["R_containment_layer2"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer2"]->GetYaxis()->SetTitle("E_{reco, #Delta R}");
    TH2_Container_["R_containment_layer3"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer3"]->GetYaxis()->SetTitle("E_{reco, #Delta R}");
    TH2_Container_["R_containment_layer4"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer4"]->GetYaxis()->SetTitle("E_{reco, #Delta R}");
    TH2_Container_["R_containment_layer5"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer5"]->GetYaxis()->SetTitle("E_{reco, #Delta R}");
    TH2_Container_["R_containment_layer6"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer6"]->GetYaxis()->SetTitle("E_{reco, #Delta R}");
    TH2_Container_["R_containment_layer7"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer7"]->GetYaxis()->SetTitle("E_{reco, #Delta R}");
    TH2_Container_["R_containment_layer8"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer8"]->GetYaxis()->SetTitle("E_{reco, #Delta R}");
    
    // TH1 histograms of total jet energy per layer
    TH1_Container_["Total_EDeposit_layer1"] = fs->make<TH1F>("Total_EDeposit_layer1", "Total E Deposit by Total in layer 1", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer2"] = fs->make<TH1F>("Total_EDeposit_layer2", "Total E Deposit by Total in layer 2", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer3"] = fs->make<TH1F>("Total_EDeposit_layer3", "Total E Deposit by Total in layer 3", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer4"] = fs->make<TH1F>("Total_EDeposit_layer4", "Total E Deposit by Total in layer 4", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer5"] = fs->make<TH1F>("Total_EDeposit_layer5", "Total E Deposit by Total in layer 5", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer6"] = fs->make<TH1F>("Total_EDeposit_layer6", "Total E Deposit by Total in layer 6", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer7"] = fs->make<TH1F>("Total_EDeposit_layer7", "Total E Deposit by Total in layer 7", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer8"] = fs->make<TH1F>("Total_EDeposit_layer8", "Total E Deposit by Total in layer 8", 40, 0, 200);
    
    // TH1 histograms of 90% energy containment radius per layer
    TH1_Container_["Ninety_Percent_R_layer1"] = fs->make<TH1F>("Ninety_Percent_R_layer1", "Ninenty Percent E Containment Radius in Layer 1", 100, 0, 100);
    TH1_Container_["Ninety_Percent_R_layer2"] = fs->make<TH1F>("Ninety_Percent_R_layer2", "Ninenty Percent E Containment Radius in Layer 2", 100, 0, 100);
    TH1_Container_["Ninety_Percent_R_layer3"] = fs->make<TH1F>("Ninety_Percent_R_layer3", "Ninenty Percent E Containment Radius in Layer 3", 100, 0, 100);
    TH1_Container_["Ninety_Percent_R_layer4"] = fs->make<TH1F>("Ninety_Percent_R_layer4", "Ninenty Percent E Containment Radius in Layer 4", 100, 0., 100);
    TH1_Container_["Ninety_Percent_R_layer5"] = fs->make<TH1F>("Ninety_Percent_R_layer5", "Ninenty Percent E Containment Radius in Layer 5", 100, 0., 100);
    TH1_Container_["Ninety_Percent_R_layer6"] = fs->make<TH1F>("Ninety_Percent_R_layer6", "Ninenty Percent E Containment Radius in Layer 6", 100, 0., 100);
    TH1_Container_["Ninety_Percent_R_layer7"] = fs->make<TH1F>("Ninety_Percent_R_layer7", "Ninenty Percent E Containment Radius in Layer 7", 100, 0., 100);
    TH1_Container_["Ninety_Percent_R_layer8"] = fs->make<TH1F>("Ninety_Percent_R_layer8", "Ninenty Percent E Containment Radius in Layer 8", 100, 0., 100);
    
    TH1_Container_["Ninety_Percent_R_layer1"]->GetXaxis()->SetTitle("R (mm)");
    TH1_Container_["Ninety_Percent_R_layer2"]->GetXaxis()->SetTitle("R (mm)");
    TH1_Container_["Ninety_Percent_R_layer3"]->GetXaxis()->SetTitle("R (mm)");
    TH1_Container_["Ninety_Percent_R_layer4"]->GetXaxis()->SetTitle("R (mm)");
    TH1_Container_["Ninety_Percent_R_layer5"]->GetXaxis()->SetTitle("R (mm)");
    TH1_Container_["Ninety_Percent_R_layer6"]->GetXaxis()->SetTitle("R (mm)");
    TH1_Container_["Ninety_Percent_R_layer7"]->GetXaxis()->SetTitle("R (mm)");
    TH1_Container_["Ninety_Percent_R_layer8"]->GetXaxis()->SetTitle("R (mm)");
    
    // TH2 histogram of energies within radius that contains 90% of energy per layer
    TH2_Container_["Ninety_Percent_E"] = fs->make<TH2F>("Ninety_Percent_E", "Energy Deposited within R_{90#%} vs Layer", 8, 1, 9, 72, 0, 180);
    TH2_Container_["Ninety_Percent_E"]->GetXaxis()->SetTitle("Layer");
    TH2_Container_["Ninety_Percent_E"]->GetYaxis()->SetTitle("E [GeV]");
    
    // Layer which cumulatively contains 90% of the total shower energy
    TH1_Container_["Ninety_Percent_Layer"] = fs->make<TH1F>("Ninety_Percent_Layer", "Layer which cumulatively contains 90#% of the deposited energy", 8, 1, 9);
    
    // Number of events which have jets developed in the layer
    TH1_Container_["Nevents_withJets_PerLayer"] = fs->make<TH1F>("Nevents_withJets_PerLayer", "Number of events with shower development, per layer", 8, 1, 9);
}


void EMShowerStudies::endJob ()
{   
    // Job ends
    std::cout << "The job, EMShowerStudies, has ended. Thank you for your patience." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (EMShowerStudies);
    
