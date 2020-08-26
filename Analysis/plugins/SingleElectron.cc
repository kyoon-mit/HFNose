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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"

// ROOT headers

using namespace edm;

SingleElectron::SingleElectron ( const edm::ParameterSet& iConfig ) :
    // histContainer_ (),
    
    // (tag name, default value (label, instance, process) -- CHECK SPELLING!!!!!!!
    tag_GenParticle_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_GenParticle", edm::InputTag ("genParticles") ) ),
    tag_TICLCandidate_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_TICLCandidate", edm::InputTag ("ticlCandidateFromTracksters") ) ),
    
    // Pre-selection parameters
    select_PID_ ( 11 ),
    select_EtaLow_ ( 3.49 ),
    select_EtaHigh_ ( 3.51 )
{
    // consumes: frequent request of additional data | mayConsume: infrequent
    token_GenParticle_ = consumes<reco::GenParticleCollection> ( tag_GenParticle_ );
    token_TICLCandidate_ = consumes<std::vector<TICLCandidate>> ( tag_TICLCandidate_ );
}


SingleElectron::~SingleElectron ()
{
    // Deconstructor
}


void SingleElectron::analyze ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{    
    // Get MC truth
    // edm::Handle<std::vector<CaloParticle>> handle_CaloParticle_MergedCaloTruth_;
    // iEvent.getByToken ( token_CaloParticle_MergedCaloTruth_, handle_CaloParticle_MergedCaloTruth_ );
    edm::Handle<reco::GenParticleCollection> handle_GenParticle;
    iEvent.getByToken ( token_GenParticle_, handle_GenParticle );
    
    // Get TICLCandidate
    edm::Handle<std::vector<TICLCandidate>> handle_TICLCandidate;
    iEvent.getByToken ( token_TICLCandidate_, handle_TICLCandidate );
    
    if ( handle_GenParticle.isValid() && handle_TICLCandidate.isValid() )
    {

        const std::vector<math::XYZTLorentzVectorF> truth_container = getTruthP4 ( *handle_GenParticle.product() );
        const std::vector<math::XYZTLorentzVectorF> reco_container = getTICLCandidateP4 ( *handle_TICLCandidate.product() );
        
        for ( auto const& truth: truth_container )
        {            
            // fillHist_HGCalRecHitsEnergy_coneR ( truth, *handle_HGCRecHits.product(), handle_HGCalGeometry.product() );
            // fillHist_CaloClustersEnergy_coneR ( truth, *handle_HGCalLayerClustersHFNose.product() );
            // histContainer_["truthEDist"]->Fill ( truth.energy() );
            std::cout << "GenParticle mass:" << truth.M() << std::endl;
        }
        
        for ( auto const& rc: reco_container )
        {
            std::cout << "TICL Candidate mass:" << rc.M() << std::endl;
        }
        
    }
    else std::cout << "Handle(s) invalid!" << std::endl;
}


std::vector<math::XYZTLorentzVectorF> SingleElectron::getTruthP4 ( const reco::GenParticleCollection & GenParticles )
{
    std::vector<math::XYZTLorentzVectorF> container;
    
    for ( auto const& gen: GenParticles )
    {
        if ( gen.pdgId() == select_PID_
                && abs(gen.eta()) > select_EtaLow_
                && abs(gen.eta()) < select_EtaHigh_ 
                && gen.isPromptFinalState() ) // Check if genparticle is final state
        {
            container.push_back ( (math::XYZTLorentzVectorF) gen.p4() );
        }
    }
    
    return container;
}


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


void SingleElectron::beginJob ()
{
    edm::Service<TFileService> fs;
    
    // histContainer_["truthEDist"] = fs->make<TH1F>("truthEDist", "Truth Energy Distribution", 200, 0, 200);
}


void SingleElectron::endJob ()
{
    // Job ends
    std::cout << "The job, SingleElectron, has ended. Thank you for your patience." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (SingleElectron);
