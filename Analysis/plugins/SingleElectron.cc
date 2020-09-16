#include "SingleElectron.h"

#include <iostream>
#include <array>

#include <cmath> // Switch to TMath.h if you need more physics-related functions

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Physics objects
// #include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
// #include "DataFormats/HGCalReco/interface/TICLCandidate.h"

// ROOT headers

using namespace edm;

SingleElectron::SingleElectron ( const edm::ParameterSet& iConfig ) :
    // histContainer_ (),
    
    // (tag name, default value (label, instance, process) -- CHECK SPELLING!!!!!!!
    // tag_GenParticle_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_GenParticle", edm::InputTag ("genParticles") ) ),
    tag_CaloParticle_MergedCaloTruth_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_MergedCaloTruth", edm::InputTag ("mix", "MergedCaloTruth") ) ),
    tag_Trackster_HFNoseTrk_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_Trackster", edm::InputTag ("ticlTrackstersHFNoseTrk") ) ),
    // tag_TICLCandidate_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_TICLCandidate", edm::InputTag ("ticlCandidateFromTracksters") ) ),
    
    // Pre-selection parameters
    select_PID_ ( 11 ),
    select_EtaLow_ ( 3.49 ),
    select_EtaHigh_ ( 3.51 )
{
    // consumes: frequent request of additional data | mayConsume: infrequent
    // token_GenParticle_ = consumes<reco::GenParticleCollection> ( tag_GenParticle_ );
    token_CaloParticle_MergedCaloTruth_ = mayConsume<std::vector<CaloParticle>> ( tag_CaloParticle_MergedCaloTruth_ );
    token_Trackster_HFNoseTrk_ = consumes<std::vector<ticl::Trackster>> ( tag_Trackster_HFNoseTrk_ );
    // token_TICLCandidate_ = consumes<std::vector<TICLCandidate>> ( tag_TICLCandidate_ );
}


SingleElectron::~SingleElectron ()
{
    // Deconstructor
}


void SingleElectron::analyze ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{    
    // Get MC truth
    // edm::Handle<reco::GenParticleCollection> handle_GenParticle;
    // iEvent.getByToken ( token_GenParticle_, handle_GenParticle );
    edm::Handle<std::vector<CaloParticle>> handle_CaloParticle_MergedCaloTruth;
    iEvent.getByToken ( token_CaloParticle_MergedCaloTruth_, handle_CaloParticle_MergedCaloTruth );
    
    // Get TICLCandidate
    // edm::Handle<std::vector<TICLCandidate>> handle_TICLCandidate;
    // iEvent.getByToken ( token_TICLCandidate_, handle_TICLCandidate );
    
    // Get ticlTrackstersHFNoseTrk
    edm::Handle<std::vector<ticl::Trackster>> handle_Trackster_HFNose_Trk;
    iEvent.getByToken ( token_Trackster_HFNoseTrk_, handle_Trackster_HFNose_Trk );
    
    if ( handle_CaloParticle_MergedCaloTruth.isValid() && handle_Trackster_HFNose_Trk.isValid() )
    {

        // const std::vector<math::XYZTLorentzVectorF> truth_container = getTruthP4 ( *handle_CaloParticle_MergedCaloTruth.product() );
        analyzeTICLTrackster ( *handle_CaloParticle_MergedCaloTruth.product(), *handle_Trackster_HFNose_Trk.product() );
        
        /*for ( auto const& truth: truth_container )
        {            
            // fillHist_HGCalRecHitsEnergy_coneR ( truth, *handle_HGCRecHits.product(), handle_HGCalGeometry.product() );
            // fillHist_CaloClustersEnergy_coneR ( truth, *handle_HGCalLayerClustersHFNose.product() );
            // histContainer_["truthEDist"]->Fill ( truth.energy() );
            
        }*/
        
    }
    else std::cout << "Handle(s) invalid!" << std::endl;
}


std::vector<math::XYZTLorentzVectorF> SingleElectron::getTruthP4 ( const std::vector<CaloParticle> & caloTruthParticles )
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


void SingleElectron::analyzeTICLTrackster ( const std::vector<CaloParticle> & caloTruthParticles, const std::vector<ticl::Trackster> & tracksters )
{

    const std::vector<math::XYZTLorentzVectorF> selected_calotruths = getTruthP4 ( caloTruthParticles );

    for ( auto const& trs: tracksters )
    {
        Float_t trackster_raw_energy = trs.raw_energy();
        Float_t trackster_eta = trs.barycenter().eta();
        
        for ( auto const& clt: selected_calotruths )
        {
            histContainer_["truthE"]->Fill( clt.E() );
            histContainer_["truthEta"]->Fill( clt.Eta() );
            
            if ( clt.Eta() > 0 && trackster_eta > 0 )
            {
                histContainer_["tracksterRawEDist"]->Fill( trackster_raw_energy );
                histContainer_["tracksterAbsEtaDist"]->Fill( abs(trackster_eta) );
                histContainer_["tracksterRawEScale"]->Fill( clt.E() - trackster_raw_energy );
                histContainer_["tracksterAbsEtaScale"]->Fill( abs(clt.Eta() - trackster_eta) );
            }
            else if ( clt.Eta() < 0 && trackster_eta < 0 )
            {
                histContainer_["tracksterRawEDist"]->Fill( trackster_raw_energy );
                histContainer_["tracksterAbsEtaDist"]->Fill( abs(trackster_eta) );
                histContainer_["tracksterRawEScale"]->Fill( clt.E() - trackster_raw_energy );
                histContainer_["tracksterAbsEtaScale"]->Fill( abs(clt.Eta() - trackster_eta) );
            }
        }
    }
}


/*
std::vector<math::XYZTLorentzVectorF> SingleElectron::getTICLCandidateP4 ( const std::vector<TICLCandidate> & TICLCandidates )
{
    std::vector<math::XYZTLorentzVectorF> container;
    
    for ( auto const& candidate: TICLCandidates )
    {
        if ( abs(candidate.pdgId()) == select_PID_
                && abs(candidate.eta()) > select_EtaLow_
                && abs(candidate.eta()) < select_EtaHigh_ )
        {
            container.push_back ( (math::XYZTLorentzVectorF) candidate.p4() );
        }
    }
    
    return container;
}
*/


void SingleElectron::beginJob ()
{
    edm::Service<TFileService> fs;
    
    histContainer_["truthE"] = fs->make<TH1F>("truthE", "Truth Energy Distribution", 100, 0, 1000);
    histContainer_["truthEta"] = fs->make<TH1F>("truthEta", "Truth Eta Distribution", 40, 3.0, 4.0);
    histContainer_["tracksterRawEDist"] = fs->make<TH1F>("tracksterRawEDist", "Trackster Raw Energy Distribution", 100, 0, 1000);
    histContainer_["tracksterAbsEtaDist"] = fs->make<TH1F>("tracksterAbsEtaDist", "Trackster Absolute Eta Distribution", 40, 3.0, 4.0);
    histContainer_["tracksterRawEScale"] = fs->make<TH1F>("tracksterRawEScale", "Trackster Raw Energy Scale", 20, 0, 100);
    histContainer_["tracksterAbsEtaScale"] = fs->make<TH1F>("tracksterAbsEtaScale", "Trackster Absolute Eta Scale", 20, 0, 0.5);
    
    // Label axes
    histContainer_["tracksterRawEDist"]->GetXaxis()->SetTitle("E_{trackster} [GeV/c^{2}]");
    histContainer_["tracksterAbsEtaDist"]->GetXaxis()->SetTitle("|#eta_{trackster}|");
    histContainer_["tracksterRawEScale"]->GetXaxis()->SetTitle("E_{trackster}-E_{caloParticle} [GeV/c^{2}]");
    histContainer_["tracksterAbsEtaScale"]->GetXaxis()->SetTitle("|#eta_{trackster} - #eta_{caloParticle}| [GeV/c^{2}]");
}


void SingleElectron::endJob ()
{
    // Job ends
    std::cout << "The job, SingleElectron, has ended. Thank you for your patience." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (SingleElectron);
