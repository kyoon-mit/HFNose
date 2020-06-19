#include <iostream>


// CMS Coordinate System
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

EMShowerStudies::EMShowerStudies ( const edm::ParameterSet& iConfig ) :

    histContainer_ (),
    
    // (tag name, default value (label, instance, process) -- CHECK SPELLING!!!!!!!    
    // Truth objects
    tag_CaloParticle_MergedCaloTruth_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_MergedCaloTruth", edm::InputTag ("mix", "MergedCaloTruth") ) ), 
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
    select_coneR_ ( 0.5 ),
    
    // Detector parameter
    HGCNose_NLayers_ ( 8 )
    
{
    // consumes: frequent request of additional data | mayConsume: infrequent
    token_CaloParticle_MergedCaloTruth_ = mayConsume<std::vector<CaloParticle>> ( tag_CaloParticle_MergedCaloTruth_ );
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
    edm::Handle<std::vector<CaloParticle>> handle_CaloParticle_MergedCaloTruth;
    iEvent.getByToken ( token_CaloParticle_MergedCaloTruth_, handle_CaloParticle_MergedCaloTruth );
    
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
    plot_maxEtaPhi ( EtaPhiE_MaxE_RecHits );
    
    // Calculate sum of E deposit of SimHits in each layer 
    const FrontBackEtaPhiE_perLayer SumE_SimHits = get_SumEDeposit_perLayer ( *handle_g4SimHits_HFNoseHits.product() );
    
    // Validate sum
    std::array<bool, 2> bool_check_sums = check_SumEDeposit_allLayers ( *handle_GenParticle.product(), SumEDeposit_perLayer );
    std::cout << "Sum validation results: " << bool_check_sums[0] << " " << bool_check_sums[1] << std::endl;
    
    // Construct radii with HGCRecHits
    void EMShowerStudies::iterative_R_search ( *handle_HGCRecHits.product(), MaxE_RecHits, SumE_SimHits );
}



FrontBackEtaPhiE_perLayer EMShowerStudies::find_EtaPhiE_MaximumEDeposit_perLayer ( const HGCRecHitCollection & hits, const HGCalGeometry * geom )
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
        Int_t layer = hit_detId.layer();
        
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


FrontBackEtaPhiE_perLayer EMShowerStudies::get_SumEDeposit_perLayer ( const std::array<PCaloHit> & HFNose_SimHits )
{ // Sum all the energies of the sim hits in each layer (should not contain noise)
    
    // Container to return
    FrontBackEtaPhiE_perLayer SumEDeposit_perLayer;

    for ( const auto& hit: HFNose_SimHits )
    {
        // Get layer number of the hit
        HFNoseDetId hit_DetId = HFNoseDetId ( hit.id() );
        Int_t layer = hit_detId.layer();
        
        // Get position of the hit
        const GlobalPoint & hit_globalPosition = geom->getPosition(hit.id());
        
        // Troubleshoot: if index out of range, doublecheck HGCNose_NLayers_ value in constructor.
        
        if ( hit_globalPosition.eta() > 0 ) SumEDeposit_perLayer[layer-1][2] += hit.energy();
        else SumEDeposit_container[layer-1][5] += hit.energy();
    }
    
    return SumEDeposit_perLayer;
}


std::array<bool, 2> EMShowerStudies::check_SumEDeposit_allLayers ( const reco::GenParticleCollection & GenParticles, const FrontBackEtaPhiE_perLayer SumEDeposit_perLayer, const Float_t fractional_error )
{ // Check if the sum of sim hits in each layer really add up to the GenParticle energy

    std::array<bool, 2> return_bool;
    return_bool.fill(false);
    
    // Sum of energies in all layers
    Float_t sum_EDeposit_front = [SumEDeposit_perLayer] ()
                                 {
                                     Float_t s = 0;
                                     for ( const auto& iter: SumEDeposit_container ) s += iter[2];
                                     return s;
                                 }
                                 
    Float_t sum_EDeposit_back = [SumEDeposit_perLayer] ()
                                 {
                                     Float_t s = 0;
                                     for ( const auto& iter: SumEDeposit_container ) s += iter[5];
                                     return s;
                                 }
    
    for ( const auto& gen: GenParticles )
    {
        if ( gen.pdgId() == select_PID_ )
        { // Select by pdgId first
            if ( gen.eta() > select_EtaLow_ && gen.eta() < select )
            { // Select by eta in the +z direction
                if ( abs(sum_EDeposit_front - gen.energy()) / gen.energy() > fractional_error ) return_bool[0] = true;
            }
            else if ( gen.eta() > -select_EtaHigh_ && gen.eta() < -select_EtaLow_ )
            { // Select by eta in the -z direction
                if ( abs(sum_EDeposit_back - gen.energy()) / gen.energy() > fractional_error ) return_bool[1] = true;
            }
        }
    }
    
    return return_bool;
}


void EMShowerStudies::iterative_R_search ( const HGCRecHitCollection & HFNose_RecoHits, const FrontBackEtaPhiE_perLayer maxE_centers, const FrontBackEtaPhiE_perLayer SumE_SimHits )
{ // Iteratively expand dR from the reference HGCRecHit and add all the HGCRecHits' energies within it. Save to TH2F histograms. (1 of 2 overloaded methods.)

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
            for ( const auto& hit: HFNose_RecoHits )
            {
            // Get layer number of the hit
            HFNoseDetId hit_DetId = HFNoseDetId ( hit.id() );
            Int_t layer = hit_detId.layer();
            
            // Get position of the hit
            const GlobalPoint & hit_globalPosition = geom->getPosition(hit.id());
            
            if ( hit_globalPosition.eta() > 0 )
            { // +z position
                Float_t dR = reco::deltaR ( hit_globalPosition.eta(), hit_globalPosition.phi(), maxE_centers[layer-1][0], maxE_centers[layer-1][1]);

                if ( dR < R ) R_EDeposit_layer_front[layer-1] += hit.energy();

            }
                        
            if ( hit_globalPosition.eta() < 0 )
            { // -z position
                Float_t dR = reco::deltaR ( hit_globalPosition.eta(), hit_globalPosition.phi(), maxE_centers[layer-1][3], maxE_centers[layer-1][4]);

                if ( dR < R ) R_EDeposit_layer_back[layer-1] += hit.energy();

            }
        }            
        
        // Fill TH2 histograms
        TH2_Container_["R_containment_layer1"]->Fill( R, R_EDeposit_layer_front[0] / SumE_SimHits[0][2] );
        TH2_Container_["R_containment_layer1"]->Fill( R, R_EDeposit_layer_back[0] / SumE_SimHits[0][5] );
        TH2_Container_["R_containment_layer2"]->Fill( R, R_EDeposit_layer_front[0] / SumE_SimHits[1][2] );
        TH2_Container_["R_containment_layer2"]->Fill( R, R_EDeposit_layer_back[0] / SumE_SimHits[1][5] );
        TH2_Container_["R_containment_layer3"]->Fill( R, R_EDeposit_layer_front[0] / SumE_SimHits[2][2] );
        TH2_Container_["R_containment_layer3"]->Fill( R, R_EDeposit_layer_back[0] / SumE_SimHits[2][5] );
        TH2_Container_["R_containment_layer4"]->Fill( R, R_EDeposit_layer_front[0] / SumE_SimHits[3][2] );
        TH2_Container_["R_containment_layer4"]->Fill( R, R_EDeposit_layer_back[0] / SumE_SimHits[3][5] );
        TH2_Container_["R_containment_layer5"]->Fill( R, R_EDeposit_layer_front[0] / SumE_SimHits[4][2] );
        TH2_Container_["R_containment_layer5"]->Fill( R, R_EDeposit_layer_back[0] / SumE_SimHits[4][5] );
        TH2_Container_["R_containment_layer6"]->Fill( R, R_EDeposit_layer_front[0] / SumE_SimHits[5][2] );
        TH2_Container_["R_containment_layer6"]->Fill( R, R_EDeposit_layer_back[0] / SumE_SimHits[5][5] );
        TH2_Container_["R_containment_layer7"]->Fill( R, R_EDeposit_layer_front[0] / SumE_SimHits[6][2] );
        TH2_Container_["R_containment_layer7"]->Fill( R, R_EDeposit_layer_back[0] / SumE_SimHits[6][5] );
        TH2_Container_["R_containment_layer8"]->Fill( R, R_EDeposit_layer_front[0] / SumE_SimHits[7][2] );
        TH2_Container_["R_containment_layer8"]->Fill( R, R_EDeposit_layer_back[0] / SumE_SimHits[7][5] );
    }
}


/*
void EMShowerStudies::iterative_R_search ( const std::vector<reco::CaloCluster> & Clusters, const FrontBackEtaPhiE_perLayer max_EtaPhiE_layer, const FrontBackEtaPhiE_perLayer SumE_SimHits )
{ // Iteratively expand dR from the reference HGCRecHit and add all the CaloClusters' energies centered within it. Save to TH2F histograms. (2 of 2 overloaded methods.)


}
*/


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
    TH2_Container_["R_containment_layer1"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #delta{R} in layer1", 20, 0., 0.5, 20, 0., 1.);
    TH2_Container_["R_containment_layer2"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #delta{R} in layer2", 20, 0., 0.5, 20, 0., 1.);
    TH2_Container_["R_containment_layer3"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #delta{R}in layer3", 20, 0., 0.5, 20, 0., 1.);
    TH2_Container_["R_containment_layer4"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #delta{R} in layer4", 20, 0., 0.5, 20, 0., 1.);
    TH2_Container_["R_containment_layer5"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #delta{R} in layer5", 20, 0., 0.5, 20, 0., 1.);
    TH2_Container_["R_containment_layer6"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #delta{R} in layer6", 20, 0., 0.5, 20, 0., 1.);
    TH2_Container_["R_containment_layer7"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #delta{R} in layer7", 20, 0., 0.5, 20, 0., 1.);
    TH2_Container_["R_containment_layer8"] = fs->make<TH2F>("R_containment_layer1", "E_{reco}/E_{gen} contained per #delta{R} in layer8", 20, 0., 0.5, 20, 0., 1.);
    
    TH2_Container_["R_containment_layer1"]->GetXaxis()->SetTitle("#delta{R}");
    TH2_Container_["R_containment_layer1"]->GetYaxis()->SetTitle("E_{reco, #delta{R}}/E_{gen}");
    TH2_Container_["R_containment_layer2"]->GetXaxis()->SetTitle("#delta{R}");
    TH2_Container_["R_containment_layer2"]->GetYaxis()->SetTitle("E_{reco, #delta{R}}/E_{gen}");
    TH2_Container_["R_containment_layer3"]->GetXaxis()->SetTitle("#delta{R}");
    TH2_Container_["R_containment_layer3"]->GetYaxis()->SetTitle("E_{reco, #delta{R}}/E_{gen}");
    TH2_Container_["R_containment_layer4"]->GetXaxis()->SetTitle("#delta{R}");
    TH2_Container_["R_containment_layer4"]->GetYaxis()->SetTitle("E_{reco, #delta{R}}/E_{gen}");
    TH2_Container_["R_containment_layer5"]->GetXaxis()->SetTitle("#delta{R}");
    TH2_Container_["R_containment_layer5"]->GetYaxis()->SetTitle("E_{reco, #delta{R}}/E_{gen}");
    TH2_Container_["R_containment_layer6"]->GetXaxis()->SetTitle("#delta{R}");
    TH2_Container_["R_containment_layer6"]->GetYaxis()->SetTitle("E_{reco, #delta{R}}/E_{gen}");
    TH2_Container_["R_containment_layer7"]->GetXaxis()->SetTitle("#delta{R}");
    TH2_Container_["R_containment_layer7"]->GetYaxis()->SetTitle("E_{reco, #delta{R}}/E_{gen}");
    TH2_Container_["R_containment_layer8"]->GetXaxis()->SetTitle("#delta{R}");
    TH2_Container_["R_containment_layer8"]->GetYaxis()->SetTitle("E_{reco, #delta{R}}/E_{gen}");
}


void EMShowerStudies::endJob ()
{
    // Job ends
    std::cout << "The job, EMShowerStudies, has ended. Thank you for your patience." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (EMShowerStudies);
    
