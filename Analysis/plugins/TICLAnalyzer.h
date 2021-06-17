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
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

// Forward declarations

class CaloParticle;
class GeomDet;
class HGCalDDDConstants;
class HGCalGeometry;
class MagneticField;
class Propagator;
class TH1F;
class TH2F;

namespace reco
{
    class Track;    
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
    std::vector<std::pair<GlobalPoint, Int_t>> getSimHitGlobalPoint ( const std::vector<CaloParticle> &, const HGCalGeometry * );
    void fillTruthHistograms ( const std::vector<math::XYZTLorentzVectorF> &, const std::vector<std::pair<GlobalPoint, Int_t>> & );
    void buildFirstLayers ( const HGCalDDDConstants * );
    void analyzeTrackPosition ( const std::vector<math::XYZTLorentzVectorF> &, const std::vector<reco::Track> &, const MagneticField*, const Propagator & );
    void analyzeTrackPosition ( const std::vector<std::pair<GlobalPoint, Int_t>> &, const std::vector<reco::Track> &, const MagneticField*, const Propagator & ); // overloaded
    void analyzeRecHits ( const std::vector<math::XYZTLorentzVectorF> &, const HGCRecHitCollection &, const HGCalGeometry * );
    void analyzeLayerClusters ( const std::vector<math::XYZTLorentzVectorF> &, const std::vector<reco::CaloCluster> & );
    void analyzeTICLTrackster ( const std::vector<math::XYZTLorentzVectorF> &, const std::vector<ticl::Trackster> &, std::string );

    // Container
    std::map <std::string, TH1F*> TH1Container_; // map of TH1 histograms
    std::map <std::string, TH2F*> TH2Container_; // map of TH2 histograms

    // ------ Data members -------
    // Tokens
    edm::EDGetTokenT<std::vector<CaloParticle>> token_CaloParticle_MergedCaloTruth_;
    edm::EDGetTokenT<std::vector<reco::Track>> token_Tracks_;
    edm::EDGetTokenT<HGCRecHitCollection> token_RecHits_HFNose_;
    edm::EDGetTokenT<HGCRecHitCollection> token_RecHits_HF_;
    edm::EDGetTokenT<HGCRecHitCollection> token_RecHits_EE_;
    edm::EDGetTokenT<std::vector<reco::CaloCluster>> token_LayerClusters_HFNose_;
    edm::EDGetTokenT<std::vector<reco::CaloCluster>> token_LayerClusters_;
    edm::EDGetTokenT<std::vector<ticl::Trackster>> token_TracksterHFNoseEM_;
    edm::EDGetTokenT<std::vector<ticl::Trackster>> token_TracksterHFNoseTrkEM_;
    edm::EDGetTokenT<std::vector<ticl::Trackster>> token_TracksterHFNoseHAD_;
    edm::EDGetTokenT<std::vector<ticl::Trackster>> token_TracksterHFNoseMIP_;
    
    // Input Tags
    edm::InputTag tag_CaloParticle_MergedCaloTruth_;
    edm::InputTag tag_Tracks_;
    edm::InputTag tag_RecHits_HFNose_;
    edm::InputTag tag_RecHits_EE_;
    edm::InputTag tag_LayerClusters_HFNose_;
    edm::InputTag tag_LayerClusters_;
    edm::InputTag tag_Trackster_;
    
    // Others
    int   select_PID_;
    double select_EtaLow_;
    double select_EtaHigh_;
    double truth_matching_deltaR_;
    std::string trackster_itername_;
    hgcal::RecHitTools rhtools_;
    std::unique_ptr<GeomDet> firstDisk_[2];
    const StringCutObjectSelector<reco::Track> cutTk_;
};

#endif
