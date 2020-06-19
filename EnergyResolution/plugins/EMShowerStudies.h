#ifndef EMSHOWERSTUDIES_H
#define EMSHOWERSTUDIES_H
#ifndef ROOT_RtypesCore
#include <RtypesCore.h>
#endif

#include <map>
#include <string>
#include <array>
#include <TH1.h>
#include <TH2.h>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// #include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

// Forward declarations


class EMShowerStudies : public edm::EDAnalyzer
{

  public:
    explicit EMShowerStudies ( const edm::ParameterSet& );
    ~EMShowerStudies ();
    
  private:
    /* This typedef is the (eta, phi, E) coordinates of two unique objects per each 
     * layer of the HGCNose, one on the +z and the other on the -z side. std::array<Float_t, 6>    
     * holds (eta1, phi1, E1, eta2, phi2, E2), where 1 stands for +z and 2 stands for -z.
     */
    typedef std::array<std::array<Float_t, 6>, HGCNose_NLayers_> FrontBackEtaPhiE_perLayer;
  
    virtual void beginJob () override;
    virtual void analyze ( const edm::Event &, const edm::EventSetup & );
    virtual void endJob () override;

    FrontBackEtaPhiE_perLayer find_EtaPhiE_MaximumEDeposit_perLayer ( const HGCRecHitCollection &, const HGCalGeometry * );
    FrontBackEtaPhiE_perLayer get_SumEDeposit_perLayer ( const std::array<PCaloHit> & );
    
    void plot_maxEtaPhi ( const FrontBackEtaPhiE_perLayer );
    
    std::array<bool, 2> check_SumEDeposit_allLayers ( const reco::GenParticleCollection &, const FrontBackEtaPhiE_perLayer, const Float_t = 0.01 );
    
    void iterative_R_search ( const HGCRecHitCollection &, const FrontBackEtaPhiE_perLayer, const FrontBackEtaPhiE_perLayer );
    void iterative_R_search ( const std::vector<reco::CaloCluster> &, const FrontBackEtaPhiE_perLayer, const FrontBackEtaPhiE_perLayer );
    
    // Container
    std::map <std::string, TH1F*> TH1_Container_; // map of TH1 histograms
    std::map <std::string, TH2F*> TH2_Container_; // map of TH1 histograms

    // Input Tags
    edm::InputTag tag_CaloParticle_MergedCaloTruth_;
    edm::InputTag tag_GenParticle_;
    edm::InputTag tag_g4SimHits_HFNoseHits_;
    edm::InputTag tag_HGCHFNoseRecHits_;
    edm::InputTag tag_HGCalLayerClustersHFNose_;
    
    // Tokens
    edm::EDGetTokenT<std::vector<CaloParticle>> token_CaloParticle_MergedCaloTruth_;
    edm::EDGetTokenT<reco::GenParticleCollection> token_GenParticle_;
    edm::EDGetTokenT<std::vector<PCaloHit> token_g4SimHits_HFNoseHits_;
    edm::EDGetTokenT<HGCRecHitCollection> token_HGCRecHits_;
    edm::EDGetTokenT<std::vector<reco::CaloCluster>> token_HGCalLayerClustersHFNose_;
    
    // Others
    Int_t   select_PID_;
    Float_t select_EtaLow_;
    Float_t select_EtaHigh_;
    Float_t select_coneR_;
    
    const Int_t HGCNose_NLayers_;
};

#endif
