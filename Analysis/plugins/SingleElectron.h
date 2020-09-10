#ifndef SINGLE_ELECTRON_H
#define SINGLE_ELECTRON_H

/*
* Single electron reconstruction analysis for HGCal using TICL
* Author: K. Yoon (Phd Candidate, MIT)
* Year: 2020
*/


// Standard libraries
#include <map>
#include <string>
#include <vector>

// ROOT Rtypes
#include <RtypesCore.h>

// EDAnalyzer base class
#include "FWCore/Framework/interface/EDAnalyzer.h"

// Other CMSSW includes that are needed in this header
#include "DataFormats/Math/interface/LorentzVector.h"

// Forward declarations
// class TH1F;

class CaloParticle;
class TICLCandidate;

namespace reco
{
    // class GenParticle;
    // typedef std::vector<GenParticle> GenParticleCollection;
}


// SingleElectron class definition
class SingleElectron : public edm::EDAnalyzer
{

  public:
    explicit SingleElectron ( const edm::ParameterSet & );
    ~SingleElectron ();
    
  private:
    virtual void beginJob () override;
    virtual void analyze ( const edm::Event&, const edm::EventSetup & );
    virtual void endJob () override;
  
    std::vector<math::XYZTLorentzVectorF> getTruthP4 ( const std::vector<CaloParticle> & );
    std::vector<math::XYZTLorentzVectorF> getTICLTracksterP4 ( );
    std::vector<math::XYZTLorentzVectorF> getTICLCandidateP4 ( const std::vector<TICLCandidate> & );

    // Container
    // std::map <std::string, TH1F*> histContainer_; // map of histograms

    // ------ Data members -------
    // Tokens
    // edm::EDGetTokenT<reco::GenParticleCollection> token_GenParticle_;
    edm::EDGetTokenT<std::vector<CaloParticle>> token_CaloParticle_MergedCaloTruth_;
    edm::EDGetTokenT<std::vector<TICLCandidate>> token_TICLCandidate_;
    
    // Input Tags
    // edm::InputTag tag_GenParticle_;
    edm::InputTag tag_CaloParticle_MergedCaloTruth_;
    edm::InputTag tag_TICLCandidate_;
    
    // Others
    Int_t   select_PID_;
    Float_t select_EtaLow_;
    Float_t select_EtaHigh_;
};

#endif
