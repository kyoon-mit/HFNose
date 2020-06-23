#include "EMShowerStudies.h"

#include <iostream>
#include <cmath> // Switch to TMath.h if you need more physics-related functions
#include "DataFormats/Math/interface/deltaR.h"

#include <TH1.h>
#include <TH2.h>

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
    select_coneR_ ( 0.5 )
    
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
    const FrontBackEtaPhiE_perLayer MaxE_RecHits = find_EtaPhiE_MaximumEDeposit_perLayer ( *handle_HGCRecHits.product(), handle_HGCalGeometry.product() ); // Why not pointer?
    
    // Plot position of above
    plot_maxEtaPhi ( MaxE_RecHits );
    
    // Calculate sum of E deposit of SimHits in each layer 
    // const FrontBackEtaPhiE_perLayer SumE_SimHits = get_SumEDeposit_perLayer ( *handle_g4SimHits_HFNoseHits.product(), handle_HGCalGeometry.product() );
    // const FrontBackEtaPhiE_perLayer SumE_SimHits = get_SumEDeposit_perLayer ( (const HGCRecHitCollection &) *handle_HGCRecHits.product(), handle_HGCalGeometry.product() );
    const FrontBackEtaPhiE_perLayer SumE_Clusters = get_SumEDeposit_perLayer ( *handle_HGCalLayerClustersHFNose.product(), handle_HGCalGeometry.product() );
    
    // Validate sum
    std::array<bool, 2> bool_check_sums = check_SumEDeposit_allLayers ( *handle_GenParticle.product(), SumE_Clusters );
    std::cout << std::boolalpha
              << "Sum validation results: " << bool_check_sums[0] << ", "
              << bool_check_sums[1] << std::endl;
    
    // Plot sum of SimHits
    plot_sum_TotalE_perLayer ( SumE_Clusters );

    // Construct radii with HGCRecHits
    iterative_R_search ( *handle_HGCalLayerClustersHFNose.product(), handle_HGCalGeometry.product(), MaxE_RecHits, SumE_Clusters );
}


EMShowerStudies::FrontBackEtaPhiE_perLayer EMShowerStudies::find_EtaPhiE_MaximumEDeposit_perLayer ( const HGCRecHitCollection & hits, const HGCalGeometry * geom )
{ // Find (eta, phi, energy) of HGCRecHit with maximum energy deposit in the calorimeter for each layer. (1 of 2 overloaded methods)

    // Containers to hold maximum energy in each layer, for comparison
    std::array<Float_t, HGCNose_NLayers_> max_E_deposit_front;
    std::array<Float_t, HGCNose_NLayers_> max_E_deposit_back;
    max_E_deposit_front.fill(0.);
    max_E_deposit_back.fill(0.);
    
    // Container to return
    FrontBackEtaPhiE_perLayer EtaPhiE_MaximumEDeposit;
    
    for ( const auto& hit : hits )
    {
        // Get layer number of the hit
        HFNoseDetId hit_DetId = HFNoseDetId ( hit.id() );
        Int_t layer = hit_DetId.layer();
        
        // Get position of the hit
        const GlobalPoint & hit_globalPosition = geom->getPosition(hit.id());
        
        // Troubleshoot: if index out of range, doublecheck HGCNose_NLayers_ value in constructor.
        
        if ( hit_globalPosition.eta() > 0 ) // +z scenario
        {   // Compare with maximum energy in same layer
            if ( max_E_deposit_front[layer-1] < hit.energy() )
            {
                max_E_deposit_front[layer-1] = hit.energy();
                EtaPhiE_MaximumEDeposit[layer-1][0] = hit_globalPosition.eta();
                EtaPhiE_MaximumEDeposit[layer-1][1] = hit_globalPosition.phi();
                EtaPhiE_MaximumEDeposit[layer-1][2] = hit.energy();
            } 
        }
        else // -z scenario
        {
            if ( max_E_deposit_back[layer-1] < hit.energy() )
            {
                max_E_deposit_back[layer-1] = hit.energy();
                EtaPhiE_MaximumEDeposit[layer-1][3] = hit_globalPosition.eta();
                EtaPhiE_MaximumEDeposit[layer-1][4] = hit_globalPosition.phi();
                EtaPhiE_MaximumEDeposit[layer-1][5] = hit.energy();
            }
        }
    }
    
    return EtaPhiE_MaximumEDeposit;
}


EMShowerStudies::FrontBackEtaPhiE_perLayer EMShowerStudies::find_EtaPhiE_MaximumEDeposit_perLayer ( const std::vector<reco::CaloCluster> & Clusters, const HGCalGeometry * geom )
{ // Find (eta, phi, energy) of CaloClusters with maximum energy deposit in the calorimeter for each layer. (2 of 2 overloaded methods)

    // Containers to hold maximum energy in each layer, for comparison
    std::array<Float_t, HGCNose_NLayers_> max_E_deposit_front;
    std::array<Float_t, HGCNose_NLayers_> max_E_deposit_back;
    max_E_deposit_front.fill(0.);
    max_E_deposit_back.fill(0.);
    
    // Container to return
    FrontBackEtaPhiE_perLayer EtaPhiE_MaximumEDeposit;
    
    for ( const auto& cl : Clusters )
    {
        // Get layer number of the cl
        HFNoseDetId cl_DetId = HFNoseDetId( cl.hitsAndFractions().at(0).first );
        Int_t layer = cl_DetId.layer();
                
        // Troubleshoot: if index out of range, doublecheck HGCNose_NLayers_ value in constructor.
        
        if ( cl.eta() > 0 ) // +z scenario
        {   // Compare with maximum energy in same layer
            if ( max_E_deposit_front[layer-1] < cl.energy() )
            {
                max_E_deposit_front[layer-1] = cl.energy();
                EtaPhiE_MaximumEDeposit[layer-1][0] = cl.eta();
                EtaPhiE_MaximumEDeposit[layer-1][1] = cl.phi();
                EtaPhiE_MaximumEDeposit[layer-1][2] = cl.energy();
            } 
        }
        else // -z scenario
        {
            if ( max_E_deposit_back[layer-1] < cl.energy() )
            {
                max_E_deposit_back[layer-1] = cl.energy();
                EtaPhiE_MaximumEDeposit[layer-1][3] = cl.eta();
                EtaPhiE_MaximumEDeposit[layer-1][4] = cl.phi();
                EtaPhiE_MaximumEDeposit[layer-1][5] = cl.energy();
            }
        }
    }
    
    return EtaPhiE_MaximumEDeposit;
}


/* Temporarily not using overloaded method with CaloParticles until it is figured out how to generate them properly

const std::array<std::array<Float_t, 3>, HGCNose_NLayers_> EMShowerStudies::find_PtEtaPhiE_MaximumEDeposit ( const std::vector<CaloParticle> & CaloParticles, const HGCalGeometry * geom )
{ // Find (eta, phi, energy) of CaloParticle with maximum energy deposit in the calorimeter for each layer. (1 of 2 overloaded methods)
    
    std::array<Float_t, HGCNose_NLayers_> max_E_deposit;
    max_E_deposit.fill(0.);

    std::array<CaloParticle, HGCNose_NLayers_> theseHits;
    std::array<std::array<Float_t, 3>, HGCNose_NLayers_> return_container;
    
    for ( const auto& cl : CaloParticles )
    {
        const SimClusterRefVector& simClusters = cl.simClusters();
        if ( cl.SimClusters().size() == 0) continue; // SimCluster must form around CaloTruth
        
        // Get layer number of the simClusters
        HFNoseDetId sc_DetId = HFNoseDetId( simClusters[0].hits_and_fractions()[0].first ); // ?????
        Int_t layer = detId.layer();
        
        // Troubleshoot: if index out of range, doublecheck HGCNose_NLayers_ value in constructor.
        
        // Compare with maximum energy in same layer
        if ( max_E_deposit[layer-1] < hit.energy() )
        {
            max_E_deposit[layer-1] = hit.energy();
            theseHits[layer-1] = hit;
        }
    }

    for ( int layer = 1; layer <= HGCNose_NLayers_; layer++ )
    { // Using index rather than iterator because need to access index in return_container
        thisHit = theseHits[layer-1];
        
        // Get coordinate of hit
        const GlobalPoint & thisHit_globalPosition = geom->getPosition(thisHit.id());
        
        // Fill return container with (eta, phi, energy)
        return_container[layer-1] = { { thisHit_globalPosition.eta(), thisHit_globalPosition.phi(), thisHit.energy() } };
    }
    
    return const_cast<std::array<std::array<Float_t, 3>, HGCNose_NLayers_>> return_container;
    
    if ( max_E_deposit < 0 ) std::cout << "findEtaPhi_CaloParticle_MaximumEDeposit: No hit found with positive energy deposit. Please revise!" << std::endl;
    
    return const_cast<CaloParticle> cl;
}
*/


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


EMShowerStudies::FrontBackEtaPhiE_perLayer EMShowerStudies::get_SumEDeposit_perLayer ( const std::vector<PCaloHit> & HFNose_SimHits, const HGCalGeometry * geom )
{ // Sum all the energies of the RecHits in each layer (1 of 3 overloaded methods)
    
    // Container to return
    FrontBackEtaPhiE_perLayer SumEDeposit_perLayer;
    SumEDeposit_perLayer.fill( std::array<Float_t, 6> { { 0., 0., 0., 0., 0., 0. } } );
    
    Float_t sum_EDeposit_front = 0;
    for ( int i = 0; i < HGCNose_NLayers_; i++ )
    {
        sum_EDeposit_front += SumEDeposit_perLayer[i][2];
    }
    std::cout << sum_EDeposit_front << std::endl;
    
    Float_t sum_EDeposit_back = 0;
    for ( int i = 0; i < HGCNose_NLayers_; i++ )
    {
        sum_EDeposit_back += SumEDeposit_perLayer[i][5];
    }
    std::cout << sum_EDeposit_back << std::endl;
    
    for ( const auto& hit: HFNose_SimHits )
    {
        // Get layer number of the hit
        HFNoseDetId hit_DetId = HFNoseDetId ( hit.id() );
        Int_t layer = hit_DetId.layer();
        
        // Get position of the hit
        const GlobalPoint & hit_globalPosition = geom->getPosition(hit.id());
        
        // Troubleshoot: if index out of range, doublecheck HGCNose_NLayers_ value in constructor.
        
        if ( hit_globalPosition.eta() > 0 ) SumEDeposit_perLayer[layer-1][2] += hit.energy();
        if ( hit_globalPosition.eta() < 0 ) SumEDeposit_perLayer[layer-1][5] += hit.energy();
    }
    
    sum_EDeposit_front = 0;
    for ( int i = 0; i < HGCNose_NLayers_; i++ )
    {
        sum_EDeposit_front += SumEDeposit_perLayer[i][2];
    }
    std::cout << sum_EDeposit_front << std::endl;
    
    sum_EDeposit_back = 0;
    for ( int i = 0; i < HGCNose_NLayers_; i++ )
    {
        sum_EDeposit_back += SumEDeposit_perLayer[i][5];
    }
    std::cout << sum_EDeposit_back << std::endl;
    
    return SumEDeposit_perLayer;
}


EMShowerStudies::FrontBackEtaPhiE_perLayer EMShowerStudies::get_SumEDeposit_perLayer ( const HGCRecHitCollection & HFNose_RecHits, const HGCalGeometry * geom )
{ // Sum all the energies of the RecHits in each layer (2 of 3 overloaded methods)
    
    // Container to return
    FrontBackEtaPhiE_perLayer SumEDeposit_perLayer;
    SumEDeposit_perLayer.fill( std::array<Float_t, 6> { { 0., 0., 0., 0., 0., 0. } } );
    
    for ( const auto& hit: HFNose_RecHits )
    {
        // Get layer number of the hit
        HFNoseDetId hit_DetId = HFNoseDetId ( hit.id() );
        Int_t layer = hit_DetId.layer();
        
        // Get position of the hit
        const GlobalPoint & hit_globalPosition = geom->getPosition(hit.id());
        
        // Troubleshoot: if index out of range, doublecheck HGCNose_NLayers_ value in constructor.
        
        if ( hit_globalPosition.eta() > 0 ) SumEDeposit_perLayer[layer-1][2] += hit.energy();
        if ( hit_globalPosition.eta() < 0 ) SumEDeposit_perLayer[layer-1][5] += hit.energy();
    }
    
    Float_t sum_EDeposit_front = 0;
    for ( int i = 0; i < HGCNose_NLayers_; i++ )
    {
        sum_EDeposit_front += SumEDeposit_perLayer[i][2];
    }
    std::cout << sum_EDeposit_front << std::endl;
    
    Float_t sum_EDeposit_back = 0;
    for ( int i = 0; i < HGCNose_NLayers_; i++ )
    {
        sum_EDeposit_back += SumEDeposit_perLayer[i][5];
    }
    std::cout << sum_EDeposit_back << std::endl;
    
    return SumEDeposit_perLayer;
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
    std::cout << sum_EDeposit_front << std::endl;
    
    Float_t sum_EDeposit_back = 0;
    for ( int i = 0; i < HGCNose_NLayers_; i++ )
    {
        sum_EDeposit_back += SumEDeposit_perLayer[i][5];
    }
    std::cout << sum_EDeposit_back << std::endl;
    
    return SumEDeposit_perLayer;
}


std::array<bool, 2> EMShowerStudies::check_SumEDeposit_allLayers ( const reco::GenParticleCollection & GenParticles, const FrontBackEtaPhiE_perLayer SumEDeposit_perLayer, const Float_t fractional_error )
{ // Check if the sum of sim hits in each layer really add up to the GenParticle energy

    std::array<bool, 2> return_bool;
    return_bool.fill(false);
    
//    auto lambda_sum_E = [SumEDeposit_perLayer] (int i)
//                 {
//                     Float_t s = 0;
//                     for ( const auto& iter: SumEDeposit_perLayer ) s += iter[i];
//                     return s;
//                 };
    
//     Sum of energies in all layers
//    Float_t sum_EDeposit_front = lambda_sum_E(2);
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
                std::cout << gen.energy() << std::endl;
                if ( abs(sum_EDeposit_front - gen.energy()) / gen.energy() < fractional_error ) return_bool[0] = true;
            }
            else if ( gen.eta() > -select_EtaHigh_ && gen.eta() < -select_EtaLow_ )
            { // Select by eta in the -z direction
                std::cout << gen.energy() << std::endl;
                if ( abs(sum_EDeposit_back - gen.energy()) / gen.energy() < fractional_error ) return_bool[1] = true;
            }
        }
    }
    
    return return_bool;
}


void EMShowerStudies::iterative_R_search ( const HGCRecHitCollection & HFNose_RecHits, const HGCalGeometry * geom, const FrontBackEtaPhiE_perLayer maxE_centers, const FrontBackEtaPhiE_perLayer TotalE_perLayer )
{ // Iteratively expand dR from the reference and add all the HGCRecHits' energies within it. Save to TH2F histograms. (1 of 2 overloaded methods.)

    for ( float R = 0.; R < 0.5; R += 0.025 )
    {
        // Containers for energy deposit within R
        std::array<Float_t, HGCNose_NLayers_> R_EDeposit_layer_front;
        std::array<Float_t, HGCNose_NLayers_> R_EDeposit_layer_back;
        R_EDeposit_layer_front.fill(0.);
        R_EDeposit_layer_back.fill(0.);
        
        if ( R == 0 )
        {
            for ( int layer = 1; layer <= HGCNose_NLayers_; layer++ )
            {   // Sum is equal to the energy of the "central" hit (i.e. hit with max E)
                R_EDeposit_layer_front[layer-1] = maxE_centers[layer-1][2];
                R_EDeposit_layer_back[layer-1]  = maxE_centers[layer-1][5];
            }
        }
        else
        {
            for ( const auto& hit: HFNose_RecHits )
            {
                // Get layer number of the hit
                HFNoseDetId hit_DetId = HFNoseDetId ( hit.id() );
                Int_t layer = hit_DetId.layer();
                
                // Get position of the hit
                const GlobalPoint & hit_globalPosition = geom->getPosition(hit.id());
                
                if ( hit_globalPosition.eta() > 0 )
                { // +z position
                    Float_t dR = reco::deltaR ( hit_globalPosition.eta(), hit_globalPosition.phi(),  maxE_centers[layer-1][0], maxE_centers[layer-1][1]);

                    if ( dR < R ) R_EDeposit_layer_front[layer-1] += hit.energy();

                }
                            
                else
                { // -z position
                    Float_t dR = reco::deltaR ( hit_globalPosition.eta(), hit_globalPosition.phi(), maxE_centers[layer-1][3], maxE_centers[layer-1][4]);

                    if ( dR < R ) R_EDeposit_layer_back[layer-1] += hit.energy();

                }
            }
        }
        
        std::cout << "R = " << R << " depo: " << R_EDeposit_layer_front[0] <<
                                        " " << R_EDeposit_layer_front[1] <<
                                        " " << R_EDeposit_layer_front[2] <<
                                        " " << R_EDeposit_layer_front[3] <<
                                        " " << R_EDeposit_layer_front[4] <<
                                        " " << R_EDeposit_layer_front[5] <<
                                        " " << R_EDeposit_layer_front[6] <<
                                        " " << R_EDeposit_layer_front[7] << std::endl;
                            
        std::cout << "R = " << R << " depo: " << R_EDeposit_layer_back[0] <<
                                        " " << R_EDeposit_layer_back[1] <<
                                        " " << R_EDeposit_layer_back[2] <<
                                        " " << R_EDeposit_layer_back[3] <<
                                        " " << R_EDeposit_layer_back[4] <<
                                        " " << R_EDeposit_layer_back[5] <<
                                        " " << R_EDeposit_layer_back[6] <<
                                        " " << R_EDeposit_layer_back[7] << std::endl;
        
        // Fill TH2 histograms
        TH2_Container_["R_containment_layer1"]->Fill( R, R_EDeposit_layer_front[0] / TotalE_perLayer[0][2] );
        TH2_Container_["R_containment_layer1"]->Fill( R, R_EDeposit_layer_back[0] / TotalE_perLayer[0][5] );
        TH2_Container_["R_containment_layer2"]->Fill( R, R_EDeposit_layer_front[1] / TotalE_perLayer[1][2] );
        TH2_Container_["R_containment_layer2"]->Fill( R, R_EDeposit_layer_back[1] / TotalE_perLayer[1][5] );
        TH2_Container_["R_containment_layer3"]->Fill( R, R_EDeposit_layer_front[2] / TotalE_perLayer[2][2] );
        TH2_Container_["R_containment_layer3"]->Fill( R, R_EDeposit_layer_back[2] / TotalE_perLayer[2][5] );
        TH2_Container_["R_containment_layer4"]->Fill( R, R_EDeposit_layer_front[3] / TotalE_perLayer[3][2] );
        TH2_Container_["R_containment_layer4"]->Fill( R, R_EDeposit_layer_back[3] / TotalE_perLayer[3][5] );
        TH2_Container_["R_containment_layer5"]->Fill( R, R_EDeposit_layer_front[4] / TotalE_perLayer[4][2] );
        TH2_Container_["R_containment_layer5"]->Fill( R, R_EDeposit_layer_back[4] / TotalE_perLayer[4][5] );
        TH2_Container_["R_containment_layer6"]->Fill( R, R_EDeposit_layer_front[5] / TotalE_perLayer[5][2] );
        TH2_Container_["R_containment_layer6"]->Fill( R, R_EDeposit_layer_back[5] / TotalE_perLayer[5][5] );
        TH2_Container_["R_containment_layer7"]->Fill( R, R_EDeposit_layer_front[6] / TotalE_perLayer[6][2] );
        TH2_Container_["R_containment_layer7"]->Fill( R, R_EDeposit_layer_back[6] / TotalE_perLayer[6][5] );
        TH2_Container_["R_containment_layer8"]->Fill( R, R_EDeposit_layer_front[7] / TotalE_perLayer[7][2] );
        TH2_Container_["R_containment_layer8"]->Fill( R, R_EDeposit_layer_back[7] / TotalE_perLayer[7][5] );
    }
}


void EMShowerStudies::iterative_R_search ( const std::vector<reco::CaloCluster> & Clusters, const HGCalGeometry * geom, const FrontBackEtaPhiE_perLayer maxE_centers, const FrontBackEtaPhiE_perLayer TotalE_perLayer )
{ // Iteratively expand dR from the reference and add all the CaloClusters' energies within it. Save to TH2F histograms. (2 of 2 overloaded methods.)

    for ( float R = 0.; R < 0.5; R += 0.025 )
    {
        // Containers for energy deposit within R
        std::array<Float_t, HGCNose_NLayers_> R_EDeposit_layer_front;
        std::array<Float_t, HGCNose_NLayers_> R_EDeposit_layer_back;
        R_EDeposit_layer_front.fill(0.);
        R_EDeposit_layer_back.fill(0.);
        
        if ( R == 0 )
        {
            for ( int layer = 1; layer <= HGCNose_NLayers_; layer++ )
            {   // Sum is equal to the energy of the "central" cluster (i.e. cluster with max E)
                R_EDeposit_layer_front[layer-1] = maxE_centers[layer-1][2];
                R_EDeposit_layer_back[layer-1]  = maxE_centers[layer-1][5];
            }
        }
        else
        {
            for ( const auto& cl: Clusters )
            {
                // Get layer number of the cluster
                HFNoseDetId cl_DetId = HFNoseDetId( cl.hitsAndFractions().at(0).first );
                Int_t layer = cl_DetId.layer();
                
                if ( cl.eta() > 0 )
                { // +z position
                    Float_t dR = reco::deltaR ( cl.eta(), cl.phi(),  maxE_centers[layer-1][0], maxE_centers[layer-1][1]);

                    if ( dR < R ) R_EDeposit_layer_front[layer-1] += cl.energy();

                }
                            
                else
                { // -z position
                    Float_t dR = reco::deltaR ( cl.eta(), cl.phi(), maxE_centers[layer-1][3], maxE_centers[layer-1][4]);

                    if ( dR < R ) R_EDeposit_layer_back[layer-1] += cl.energy();

                }
            }
        }
        
//        std::cout << "R = " << R << " depo: " << R_EDeposit_layer_front[0] <<
//                                        " " << R_EDeposit_layer_front[1] <<
//                                        " " << R_EDeposit_layer_front[2] <<
//                                        " " << R_EDeposit_layer_front[3] <<
//                                        " " << R_EDeposit_layer_front[4] <<
//                                        " " << R_EDeposit_layer_front[5] <<
//                                        " " << R_EDeposit_layer_front[6] <<
//                                        " " << R_EDeposit_layer_front[7] << std::endl;
//                            
//        std::cout << "R = " << R << " depo: " << R_EDeposit_layer_back[0] <<
//                                        " " << R_EDeposit_layer_back[1] <<
//                                        " " << R_EDeposit_layer_back[2] <<
//                                        " " << R_EDeposit_layer_back[3] <<
//                                        " " << R_EDeposit_layer_back[4] <<
//                                        " " << R_EDeposit_layer_back[5] <<
//                                        " " << R_EDeposit_layer_back[6] <<
//                                        " " << R_EDeposit_layer_back[7] << std::endl;
        
        // Fill TH2 histograms
        TH2_Container_["R_containment_layer1"]->Fill( R, R_EDeposit_layer_front[0] / TotalE_perLayer[0][2] );
        TH2_Container_["R_containment_layer1"]->Fill( R, R_EDeposit_layer_back[0] / TotalE_perLayer[0][5] );
        TH2_Container_["R_containment_layer2"]->Fill( R, R_EDeposit_layer_front[1] / TotalE_perLayer[1][2] );
        TH2_Container_["R_containment_layer2"]->Fill( R, R_EDeposit_layer_back[1] / TotalE_perLayer[1][5] );
        TH2_Container_["R_containment_layer3"]->Fill( R, R_EDeposit_layer_front[2] / TotalE_perLayer[2][2] );
        TH2_Container_["R_containment_layer3"]->Fill( R, R_EDeposit_layer_back[2] / TotalE_perLayer[2][5] );
        TH2_Container_["R_containment_layer4"]->Fill( R, R_EDeposit_layer_front[3] / TotalE_perLayer[3][2] );
        TH2_Container_["R_containment_layer4"]->Fill( R, R_EDeposit_layer_back[3] / TotalE_perLayer[3][5] );
        TH2_Container_["R_containment_layer5"]->Fill( R, R_EDeposit_layer_front[4] / TotalE_perLayer[4][2] );
        TH2_Container_["R_containment_layer5"]->Fill( R, R_EDeposit_layer_back[4] / TotalE_perLayer[4][5] );
        TH2_Container_["R_containment_layer6"]->Fill( R, R_EDeposit_layer_front[5] / TotalE_perLayer[5][2] );
        TH2_Container_["R_containment_layer6"]->Fill( R, R_EDeposit_layer_back[5] / TotalE_perLayer[5][5] );
        TH2_Container_["R_containment_layer7"]->Fill( R, R_EDeposit_layer_front[6] / TotalE_perLayer[6][2] );
        TH2_Container_["R_containment_layer7"]->Fill( R, R_EDeposit_layer_back[6] / TotalE_perLayer[6][5] );
        TH2_Container_["R_containment_layer8"]->Fill( R, R_EDeposit_layer_front[7] / TotalE_perLayer[7][2] );
        TH2_Container_["R_containment_layer8"]->Fill( R, R_EDeposit_layer_back[7] / TotalE_perLayer[7][5] );
    }
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
    
    // TH2 histograms about fraction of energy included within cylinder radius
    TH2_Container_["R_containment_layer1"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #Delta R in layer1", 20, 0., 0.5, 50, 0., 0.5);
    TH2_Container_["R_containment_layer2"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #Delta R in layer2", 20, 0., 0.5, 50, 0., 0.5);
    TH2_Container_["R_containment_layer3"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #Delta Rin layer3", 20, 0., 0.5, 50, 0., 0.5);
    TH2_Container_["R_containment_layer4"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #Delta R in layer4", 20, 0., 0.5, 50, 0., 0.5);
    TH2_Container_["R_containment_layer5"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #Delta R in layer5", 20, 0., 0.5, 50, 0., 0.5);
    TH2_Container_["R_containment_layer6"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #Delta R in layer6", 20, 0., 0.5, 50, 0., 0.5);
    TH2_Container_["R_containment_layer7"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #Delta R in layer7", 20, 0., 0.5, 50, 0., 0.5);
    TH2_Container_["R_containment_layer8"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #Delta R in layer8", 20, 0., 0.5, 50, 0., 0.5);
    
    TH2_Container_["R_containment_layer1"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer1"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_containment_layer2"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer2"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_containment_layer3"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer3"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_containment_layer4"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer4"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_containment_layer5"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer5"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_containment_layer6"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer6"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_containment_layer7"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer7"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    TH2_Container_["R_containment_layer8"]->GetXaxis()->SetTitle("#Delta R");
    TH2_Container_["R_containment_layer8"]->GetYaxis()->SetTitle("E_{reco, #Delta R}/E_{gen}");
    
    // TH1 histograms of total jet energy per layer
    TH1_Container_["Total_EDeposit_layer1"] = fs->make<TH1F>("Total_EDeposit_layer1", "Total E Deposit by Total in layer 1", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer2"] = fs->make<TH1F>("Total_EDeposit_layer2", "Total E Deposit by Total in layer 2", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer3"] = fs->make<TH1F>("Total_EDeposit_layer3", "Total E Deposit by Total in layer 3", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer4"] = fs->make<TH1F>("Total_EDeposit_layer4", "Total E Deposit by Total in layer 4", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer5"] = fs->make<TH1F>("Total_EDeposit_layer5", "Total E Deposit by Total in layer 5", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer6"] = fs->make<TH1F>("Total_EDeposit_layer6", "Total E Deposit by Total in layer 6", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer7"] = fs->make<TH1F>("Total_EDeposit_layer7", "Total E Deposit by Total in layer 7", 40, 0, 200);
    TH1_Container_["Total_EDeposit_layer8"] = fs->make<TH1F>("Total_EDeposit_layer8", "Total E Deposit by Total in layer 8", 40, 0, 200);
    
    // TH1 histograms of total RecHits energy per layer
    TH1_Container_["RecHits_EDeposit_layer1"] = fs->make<TH1F>("RecHits_EDeposit_layer1", "Total E Deposit by RecHits in layer 1", 40, 0, 200);
    TH1_Container_["RecHits_EDeposit_layer2"] = fs->make<TH1F>("RecHits_EDeposit_layer2", "Total E Deposit by RecHits in layer 2", 40, 0, 200);
    TH1_Container_["RecHits_EDeposit_layer3"] = fs->make<TH1F>("RecHits_EDeposit_layer3", "Total E Deposit by RecHits in layer 3", 40, 0, 200);
    TH1_Container_["RecHits_EDeposit_layer4"] = fs->make<TH1F>("RecHits_EDeposit_layer4", "Total E Deposit by RecHits in layer 4", 40, 0, 200);
    TH1_Container_["RecHits_EDeposit_layer5"] = fs->make<TH1F>("RecHits_EDeposit_layer5", "Total E Deposit by RecHits in layer 5", 40, 0, 200);
    TH1_Container_["RecHits_EDeposit_layer6"] = fs->make<TH1F>("RecHits_EDeposit_layer6", "Total E Deposit by RecHits in layer 6", 40, 0, 200);
    TH1_Container_["RecHits_EDeposit_layer7"] = fs->make<TH1F>("RecHits_EDeposit_layer7", "Total E Deposit by RecHits in layer 7", 40, 0, 200);
    TH1_Container_["RecHits_EDeposit_layer8"] = fs->make<TH1F>("RecHits_EDeposit_layer8", "Total E Deposit by RecHits in layer 8", 40, 0, 200);
}


void EMShowerStudies::endJob ()
{
    // Job ends
    std::cout << "The job, EMShowerStudies, has ended. Thank you for your patience." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (EMShowerStudies);
    
