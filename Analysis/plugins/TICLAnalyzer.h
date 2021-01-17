#ifndef TICL_ANALYZER_H
#define TICL_ANALYZER_H

/*
* Author: K. Yoon
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
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

// Forward declarations
class TH1F;

class CaloParticle;
// class TICLCandidate;
class HGCalGeometry;

/*namespace edm*/
/*{*/
/*    class PCaloHit;*/
/*    typedef std::vector<PCaloHit> PCaloHitContainer;*/
/*}*/

namespace reco
{
    class GenParticle;
    typedef std::vector<GenParticle> GenParticleCollection;
    
    class CaloCluster;
}

namespace ticl
{
    class Trackster;
}


// TICLAnalyzer class definition
class TICLAnalyzer : public edm::EDAnalyzer
{

  public:
    explicit TICLAnalyzer ( const edm::ParameterSet & );
    ~TICLAnalyzer ();
    
  private:
    virtual void beginJob () override;
    virtual void analyze ( const edm::Event&, const edm::EventSetup & );
    virtual void endJob () override;
  
    std::vector<math::XYZTLorentzVectorF> getTruthP4 ( const std::vector<CaloParticle> & );
    void analyzeTICLTrackster ( const std::vector<CaloParticle> &, const std::vector<ticl::Trackster> &, std::string );
/*    void analyzeSimHits ( const std::vector<CaloParticle> &, const edm::PCaloHitContainer &, const HGCalGeometry* );*/
    void analyzeRecHits ( const std::vector<CaloParticle> &, const HGCRecHitCollection &, const HGCalGeometry* );
    void analyzeLayerClusters ( const std::vector<CaloParticle> &, const std::vector<reco::CaloCluster> & );

    // Container
    std::map <std::string, TH1F*> histContainer_; // map of histograms

    // ------ Data members -------
    // Tokens
    edm::EDGetTokenT<std::vector<CaloParticle>> token_CaloParticle_MergedCaloTruth_;
    edm::EDGetTokenT<reco::GenParticleCollection> token_GenParticle_;
    
/*    edm::EDGetTokenT<edm::PCaloHitContainer> token_SimHits_HFNose_;*/
    edm::EDGetTokenT<HGCRecHitCollection> token_RecHits_HFNose_;
    edm::EDGetTokenT<std::vector<reco::CaloCluster>> token_LayerClusters_HFNose_;
    edm::EDGetTokenT<std::vector<ticl::Trackster>> token_Trackster_;
    
    // Input Tags
    // edm::InputTag tag_GenParticle_;
    edm::InputTag tag_CaloParticle_MergedCaloTruth_;
    edm::InputTag tag_GenParticle_;
    
/*    edm::InputTag tag_SimHits_HFNose_;*/
    edm::InputTag tag_RecHits_HFNose_;
    edm::InputTag tag_LayerClusters_HFNose_;
    edm::InputTag tag_Trackster_;
    
    // Others
    std::string trackster_itername_;
    int   select_PID_;
    double select_EtaLow_;
    double select_EtaHigh_;
    double truth_matching_deltaR_;
};

#endif
