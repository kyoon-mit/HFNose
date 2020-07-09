#ifndef MOLIERERADIUS_H
#define MOLIERERADIUS_H

// Standard libraries
#include <map>
#include <string>
#include <array>
#include <vector>

// ROOT Rtypes
#include <RtypesCore.h>

// EDAnalyzer base class
#include "FWCore/Framework/interface/EDAnalyzer.h"

// Other CMSSW includes that are needed in this header
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

// Forward declarations
class TH1F;
class TH2F;

class PCaloHit;
class HGCalGeometry;

namespace reco
{
    class GenParticle;
    typedef std::vector<GenParticle> GenParticleCollection;
    class CaloCluster;
}


// MoliereRadius class definition
class MoliereRadius : public edm::EDAnalyzer
{

  public:
    explicit MoliereRadius ( const edm::ParameterSet& );
    ~MoliereRadius ();
    
  private:
    virtual void beginJob () override;
    virtual void analyze ( const edm::Event &, const edm::EventSetup & );
    virtual void endJob () override;
    
    /* This typedef is the (eta, phi, E) coordinates of two unique objects per each 
     * layer of the HGCNose, one on the +z and the other on the -z side. std::array<Float_t, 6>    
     * holds (eta1, phi1, E1, eta2, phi2, E2), where 1 stands for +z and 2 stands for -z.
     */
    const static Int_t HGCNose_NLayers_ = 8; // Number of layers in HGCNose
    typedef std::array<std::array<Float_t, 6>, HGCNose_NLayers_> FrontBackEtaPhiE_perLayer;

    // Functions
    // Rule for ordering parameters: 1. Event related, 2. EventSetup related, 3. others
    FrontBackEtaPhiE_perLayer find_EtaPhiE_Reference_perLayer ( const reco::GenParticleCollection &, const HGCalGeometry * );
    
    FrontBackEtaPhiE_perLayer get_SumEDeposit_perLayer ( const std::vector<reco::CaloCluster> &, const HGCalGeometry * );
    
    void plot_maxEtaPhi ( const FrontBackEtaPhiE_perLayer );
    
    void plot_sum_TotalE_perLayer ( const FrontBackEtaPhiE_perLayer );
    
    std::array<bool, 2> check_SumEDeposit_allLayers ( const reco::GenParticleCollection &, const FrontBackEtaPhiE_perLayer, const Float_t = 0.05 );
    
    Float_t getContainmentR ( const std::vector<Float_t>, const std::vector<Float_t>, const Float_t = 0.9 );
    
    std::array<Int_t, 2> getContainmentLayer ( const FrontBackEtaPhiE_perLayer, const FrontBackEtaPhiE_perLayer, Float_t FracContainment = 0.9 );
    
    FrontBackEtaPhiE_perLayer getContainedEnergy ( const std::vector<reco::CaloCluster> &, const HGCalGeometry *, const FrontBackEtaPhiE_perLayer, const FrontBackEtaPhiE_perLayer );
    
    void iterative_R_search ( const std::vector<reco::CaloCluster> &, const HGCalGeometry *, const FrontBackEtaPhiE_perLayer, const FrontBackEtaPhiE_perLayer);
    
    // Container
    std::map <std::string, TH1F*> TH1_Container_; // map of TH1 histograms
    std::map <std::string, TH2F*> TH2_Container_; // map of TH2 histograms

    // Input Tags
    // edm::InputTag tag_CaloParticle_MergedCaloTruth_;
    edm::InputTag tag_GenParticle_;
    edm::InputTag tag_g4SimHits_HFNoseHits_;
    edm::InputTag tag_HGCHFNoseRecHits_;
    edm::InputTag tag_HGCalLayerClustersHFNose_;
    
    // Tokens
    // edm::EDGetTokenT<std::vector<CaloParticle>> token_CaloParticle_MergedCaloTruth_;
    edm::EDGetTokenT<reco::GenParticleCollection> token_GenParticle_;
    edm::EDGetTokenT<std::vector<PCaloHit>> token_g4SimHits_HFNoseHits_;
    edm::EDGetTokenT<HGCRecHitCollection> token_HGCRecHits_;
    edm::EDGetTokenT<std::vector<reco::CaloCluster>> token_HGCalLayerClustersHFNose_;
    
    // Others
    Int_t   select_PID_;
    Float_t select_EtaLow_;
    Float_t select_EtaHigh_;
    Float_t max_iter_R_;
    Int_t steps_iter_R_;
};

#endif
