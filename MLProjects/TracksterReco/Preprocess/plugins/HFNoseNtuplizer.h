#ifndef HFNOSE_NTUPLIZER_H
#define HFNOSE_NTUPLIZER_H

/*
* Original Author: Bruno Alves
* Original Repo: https://gitlab.cern.ch/gouskos/hgcal_reco_analysis/HGCMLAnalyzer/plugins/Pid.cc
* Modified By: K. Yoon
* Modified In: 2021
* Plugin for creating ROOT Ntuple tree for preprocessing
*/

// Standard libraries
#include <map>
#include <vector>

#include "CommonDataFormats.h"

// Forward declarations
class CaloParticle;
class DetId;
class HGCRecHit;
typedef edm::SortedCollection<HGCRecHit> HGCRecHitCollection;
class layercluster;

namespace hgcal {
    class RecHitTools;
}

namespace reco {
    class CaloCluster;
    typedef std::vector<CaloCluster> CaloClusterCollection;
}

namespace ticl {
    class Trackster;
}


// HFNoseNtuplizer class definition
class HFNoseNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

  public:
    explicit HFNoseNtuplizer ( const edm::ParameterSet& );
    ~HFNoseNtuplizer() override;

    static void fillDescriptions( edm::ConfigurationDescriptions& );

  private:
    void beginJob() override;
    void analyze( const edm::Event &, const edm::EventSetup & ) override;
    void endJob() override;
    
    virtual void fillHitMap( std::map<DetId, const HGCRecHit*> &, 
                             const HGCRecHitCollection & ) const;
    std::vector<int> getClosestTrackstersToCPByIndex( const caloparticle &,
                                                      const std::vector<ticl::Trackster> &,
                                                      float );
    int isRecHitMatchedToCPRecHits( DetId, const std::vector<DetId>& );
    float getTracksterEnergyFromCP( const ticl::Trackster &,
                                    const reco::CaloClusterCollection &,
                                    const caloparticle & );
    std::vector<layercluster> getLCsFromTrackster( const ticl::Trackster &, const reco::CaloClusterCollection & );

    std::shared_ptr<hgcal::RecHitTools>;

    // ----------member data ---------------------------
    edm::EDGetTokenT<std::vector<CaloParticle>> token_CaloParticles_;
    edm::EDGetTokenT<HGCRecHitCollection> token_HFNoseRecHits_;
    edm::EDGetTokenT<std::vector<ticl::Trackster>> token_Tracksters_;
    edm::EDGetTokenT<reco::CaloClusterCollection> token_HFNoseLayerClusters_;

    TTree* tree;// = new TTree("tree", "tree");

    edm::RunNumber_t irun;
    edm::EventNumber_t ievent;
    edm::LuminosityBlockNumber_t ilumiblock;
    edm::Timestamp itime;

    //size_t run, event, lumi, time;
    size_t run, lumi, time;
    int event;
    bool withPU_;
    float ts_energy;
    float ts_sigma1;
    float ts_sigma2;
    float ts_sigma3;
    float cp_missingEnergyFraction;
    float cp_energy;
    float cp_eta;
    float cp_phi;
    int cp_pdgid;
    int trackster;

    std::vector<float> lc_energy;
    std::vector<float> lc_eta;
    std::vector<float> lc_phi;
    std::vector<int> lc_layer;
};

#endif
