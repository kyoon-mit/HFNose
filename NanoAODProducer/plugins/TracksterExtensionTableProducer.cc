#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Common/interface/View.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"

#include <vector>
#include <array>
#include <iostream>

class TracksterExtensionTableProducer : public edm::stream::EDProducer<> {
public:
  TracksterExtensionTableProducer(edm::ParameterSet const& params)
      : name_(params.getParameter<std::string>("name")),
        doc_(params.getParameter<std::string>("doc")),
        src_(consumes<std::vector<ticl::Trackster>>(params.getParameter<edm::InputTag>("src"))),
        cut_(params.getParameter<std::string>("cut"), true) {
    produces<nanoaod::FlatTable>();
  }

  ~TracksterExtensionTableProducer() override {}

  void beginRun(const edm::Run&, const edm::EventSetup& iSetup) {
    // TODO: check that the geometry exists
  }
  
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    edm::Handle<std::vector<ticl::Trackster>> objs;
    iEvent.getByToken(src_, objs);
    
    std::vector<float> sigma1;
    std::vector<float> sigma2;
    std::vector<float> sigma3;
    std::vector<float> sigmaPCA1;
    std::vector<float> sigmaPCA2;
    std::vector<float> sigmaPCA3;
    
    for (auto const& obj : *objs.product()) {
      if (cut_(obj)) {
        std::array<float, 3> sigmas = obj.sigmas();
        std::array<float, 3> sigmasPCA = obj.sigmasPCA();
        // std::array<float, 8> id_probabilities = obj.id_probabilities();
        sigma1.emplace_back(sigmas[0]);
        sigma2.emplace_back(sigmas[1]);
        sigma3.emplace_back(sigmas[2]);
        sigmaPCA1.emplace_back(sigmasPCA[0]);
        sigmaPCA2.emplace_back(sigmasPCA[1]);
        sigmaPCA3.emplace_back(sigmasPCA[2]);
      }
    }

    
        auto tab = std::make_unique<nanoaod::FlatTable>(sigma1.size(), name_, false, true);
        tab->addColumn<float>("sigma1", sigma1, "first component of sigmas");
        tab->addColumn<float>("sigma2", sigma2, "second component of sigmas");
        tab->addColumn<float>("sigma3", sigma3, "third component of sigmas");
        tab->addColumn<float>("sigmaPCA1", sigmaPCA1, "first component of PCA sigmas");
        tab->addColumn<float>("sigmaPCA2", sigmaPCA2, "second component of PCA sigmas");
        tab->addColumn<float>("sigmaPCA3", sigmaPCA3, "third component of PCA sigmas");

        iEvent.put(std::move(tab));
  }

protected:
  const std::string name_, doc_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> src_;
  const StringCutObjectSelector<ticl::Trackster> cut_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TracksterExtensionTableProducer);
