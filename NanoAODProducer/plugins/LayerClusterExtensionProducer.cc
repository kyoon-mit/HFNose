#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Common/interface/View.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

#include <vector>
#include <iostream>

class LayerClusterExtensionProducer : public edm::stream::EDProducer<> {
public:
  LayerClusterExtensionProducer(edm::ParameterSet const& params)
      : name_(params.getParameter<std::string>("name")),
        doc_(params.getParameter<std::string>("doc")),
        src_(consumes<std::vector<reco::CaloCluster>>(params.getParameter<edm::InputTag>("src"))),
        cut_(params.getParameter<std::string>("cut"), true) {
    produces<nanoaod::FlatTable>();
  }

  ~LayerClusterExtensionProducer() override {}

  void beginRun(const edm::Run&, const edm::EventSetup& iSetup) {
  }
  
  float layerFromCluster(const reco::CaloCluster& cluster) {
    hgcal::RecHitTools rhtools_;
    return rhtools_.getLayer(cluster.hitsAndFractions().at(0).first);
  }
  
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    edm::Handle<std::vector<reco::CaloCluster>> objs;
    iEvent.getByToken(src_, objs);
    
    std::vector<int> layernums;
    for (const auto& obj : *objs.product()) {
      if (cut_(obj)) {
        layernums.emplace_back(layerFromCluster(obj));
      }
    }

    auto tab = std::make_unique<nanoaod::FlatTable>(layernums.size(), name_, false, true);
    tab->addColumn<int>("layer", layernums, "layer number");

    iEvent.put(std::move(tab));
  }

protected:
  const std::string name_, doc_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> src_;
  const StringCutObjectSelector<reco::CaloCluster> cut_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LayerClusterExtensionProducer);
