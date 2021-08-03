#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"f
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

//ROOT includes
//#include "TTree.h"
//#include <TFile.h>
//#include <TROOT.h>
//#include "TBranch.h"
//#include <string>
//#include "TSystem.h"
//#include "TVector3.h"
//#include "TH1.h"


bool sortLCsByEnergyAndLayer( const layercluster& a, const layercluster& b ) {
  return (a.energy_ > b.energy_) || ((a.energy_ == b.energy_) && (a.layer_ < b.layer_));
}


HFNoseNtuplizer::HFNoseNtuplizer(const edm::ParameterSet& iConfig) :
    token_CaloParticles_( consumes<std::vector<CaloParticle>>( iConfig.getParameter<edm::InputTag>("mix", "MergedCaloTruth", "HLT") ) ),
    token_HFNoseRecHits_( consumes<HGCRecHitCollection>( iConfig.getParameter<edm::InputTag>("HGCHFNoseRecHits") ) ),
    token_Tracksters_( consumes<std::vector<ticl::Trackster>>( iConfig.getParameter<edm::InputTag>("mergeTrackster") ) ),
    token_HFNoseLayerClusters_( consumes<reco::CaloClusterCollection>( iConfig.getParameter<edm::InputTag>("hgcalLayerClusters") ) ),
    withPU_( iConfig.getParameter<bool>("withPU") )
{
    recHitTools.reset( new hgcal::RecHitTools() );
    //now do what ever initialization is needed
    usesResource("TFileService");
    edm::Service<TFileService> file;

    tree = file->make<TTree>("pidtree", "pidtree");

    tree->Branch("run", &run, "run/I");
    tree->Branch("event", &event, "event/I");
    tree->Branch("lumi", &lumi, "lumi/I");
    tree->Branch("ts_energy", &ts_energy, "ts_energy/F");
    tree->Branch("ts_sigma1", &ts_sigma1, "ts_sigma1/F");
    tree->Branch("ts_sigma2", &ts_sigma2, "ts_sigma2/F");
    tree->Branch("ts_sigma3", &ts_sigma3, "ts_sigma3/F");
    tree->Branch("cp_missingEnergyFraction", &cp_missingEnergyFraction, "cp_missingEnergyFraction/F");
    tree->Branch("cp_energy", &cp_energy, "cp_energy/F");
    tree->Branch("cp_eta", &cp_eta, "cp_eta/F");
    tree->Branch("cp_phi", &cp_phi, "cp_phi/F");
    tree->Branch("cp_pdgid", &cp_pdgid, "cp_pdgid/I");
    tree->Branch("trackster", &trackster, "trackster/I");
    tree->Branch("lc_energy", &lc_energy);
    tree->Branch("lc_eta", &lc_eta);
    tree->Branch("lc_phi", &lc_phi);
    tree->Branch("lc_layer", &lc_layer);
}


HFNoseNtuplizer::~HFNoseNtuplizer() {}


void HFNoseNtuplizer::fillHitMap( std::map<DetId, const HGCRecHit*> & hitMap,
                                  const HGCRecHitCollection & recHitsHFNose ) const
{
    hitMap.clear();
    for ( const auto& hit : recHitsHFNose ) hitMap.emplace( hit.detid(), &hit );
}


std::vector<int> HFNoseNtuplizer::getClosestTrackstersToCPByIndex( const caloparticle & cp,
                                                       const std::vector<ticl::Trackster> & tracksters,
                                                       float maxDrTracksterCP )
{
    std::vector<int> closestTracksters_;
    int idx = 0;
    for ( auto const& t : tracksters )
    {
        float dr = reco::deltaR( cp.eta_, cp.phi_, t.barycenter().eta(), t.barycenter().phi() );
        if ( dr < maxDrTracksterCP ) closestTracksters_.push_back(idx);
        idx++;
    }
    return closestTracksters_;
}


int HFNoseNtuplizer::isRecHitMatchedToCPRecHits( DetId detid_, const std::vector<DetId> & rechitdetid_ )
{
    auto found = std::find( std::begin(rechitdetid_), std::end(rechitdetid_), detid_ );
    return ( found != rechitdetid_.end() ) ? std::distance( std::begin(rechitdetid_), found ) : -1;
}


float HFNoseNtuplizer::getTracksterEnergyFromCP( const ticl::Trackster & trackster,
                                                 const reco::CaloClusterCollection & lcs,
                                                 const caloparticle & cp )
{
    //  std::cout << " IN getTracksterEnFromCP \n";
    float enFromRecHits_ = 0.;

    // get the indices associated to the LCs of this trackster
    const std::vector<unsigned int>& lcIdxs_ = trackster.vertices();
    // loop over these idxs
    for ( unsigned int ilc = 0; ilc < lcIdxs_.size(); ++ilc )
    {
        // get the lc and the corresponding rechits to this lc
        const reco::BasicCluster& lc = lcs[lcIdxs_[ilc]];
        auto const& hf = lc.hitsAndFractions();
        // loop over the rechits of this specific layer cluster
        for (unsigned int j = 0; j < hf.size(); j++)
        {
            const DetId detid_ = hf[j].first;
            int detid_idx = isRecHitMatchedToCPRecHits(detid_, cp.rechitdetid_);
            if (detid_idx >= 0) enFromRecHits_ += cp.rechitenergy_[detid_idx];
        }  // end of looping over the rechits
    }    // end of looping over the idxs

    return enFromRecHits_;
}  // end of getTracksterEnFromCP


std::vector<layercluster> HFNoseNtuplizer::getLCsFromTrackster( const ticl::Trackster & trackster,
                                                                            const reco::CaloClusterCollection & lcs )
{
    std::vector<layercluster> layerclusters_;
    const std::vector<unsigned int> & lcIdxs_ = trackster.vertices();

    for ( unsigned int ilc = 0; ilc < lcIdxs_.size(); ++ilc )
    {
        const reco::BasicCluster & lc = lcs.at(lcIdxs_[ilc]);
        auto const & hf = lc.hitsAndFractions();

        int layer_ = recHitTools->getLayerWithOffset(hf[0].first);

        layerclusters_.push_back(layercluster());
        auto& layercluster_ = layerclusters_.back();
        layercluster_.energy_ = lc.energy() / (float)trackster.vertex_multiplicity(ilc);
        layercluster_.eta_ = lc.eta();
        layercluster_.phi_ = lc.phi();
        layercluster_.x_ = lc.position().x();
        layercluster_.y_ = lc.position().y();
        layercluster_.z_ = lc.position().z();
        layercluster_.nrechits_ = lc.hitsAndFractions().size();
        layercluster_.layer_ = abs(layer_);
    }

    return layerclusters_;
}  // end of getLCsFromTrackster


void HFNoseNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<HGCRecHitCollection> handle_HFNoseRecHits;
    iEvent.getByToken( token_HFNoseRecHits_, handle_HFNoseRecHits );

    edm::Handle<std::vector<ticl::Trackster>> handle_Tracksters;
    iEvent.getByToken( token_Tracksters_, handle_Tracksters );

    edm::Handle<std::vector<CaloParticle>> handle_CaloParticles;
    iEvent.getByToken( token_CaloParticles_, handle_CaloParticles );

    edm::Handle<reco::CaloClusterCollection> handle_HFNoseLayerClusters;
    iEvent.getByToken( token_HFNoseLayerClusters_, handle_HFNoseLayerClusters );

    std::map<DetId, const HGCRecHit*> hitMap;
    fillHitMap( hitMap, *handle_HFNoseRecHits );

    // init vars
    edm::ESHandle<HGCalGeometry> handle_HGCalGeometry; 
    iSetup.get<IdealGeometryRecord>().get( "HGCalHFNoseSensitive", handle_HGCalGeometry ); 
    recHitTools->setGeometry(*handle_HGCalGeometry);

    const std::vector<ticl::Trackster> & prod_tracksters = *handle_Tracksters;
    const CaloParticleCollection & cps = *handle_CaloParticles;
    const reco::CaloClusterCollection & lcs = *handle_HFNoseLayerClusters;

    // get edm CaloParticles & convert to caloparticles
    std::vector<caloparticle> caloparticles;
    for (const auto& it_cp : cps)
    {
        const CaloParticle& cp = ((it_cp));

        //only the main calo particle, and not the PU particles
        bool cp_cut = withPU_ ? cp.eventId().event() != 0 or cp.eventId().bunchCrossing() != 0 : cp.eventId().bunchCrossing() != 0;
        if (cp_cut) continue;

        // allow CP only within HFNose volume
        if ( !(3.0 <= abs(cp.eta()) <= 4.2) ) continue;

        caloparticle tmpcp_;
        tmpcp_.pdgid_ = cp.pdgId();
        tmpcp_.energy_ = cp.energy();
        tmpcp_.pt_ = cp.pt();
        tmpcp_.eta_ = cp.eta();
        tmpcp_.phi_ = cp.phi();

        std::vector<DetId> tmprechits_;
        std::vector<float> tmprechitenergy_;

        // get the simclusters
        const SimClusterRefVector& simclusters = cp.simClusters();
        for (const auto& it_simc : simclusters) {
            const SimCluster& simc = (*(it_simc));
            const auto& sc_haf = simc.hits_and_fractions();

            // get the rechits
            for (const auto& it_sc_haf : sc_haf)
            {
                DetId detid_ = (it_sc_haf.first);
                // need to map RecHits to the SimCluster
                // SimHits are not stored
                auto const itcheck = hitMap.find(detid_);
                // we need this check because some DetIDs assigned to CaloParticle do not always have
                // a RecHit -- due to thresholds or whatever
                if ( itcheck != hitMap.end() )
                {
                    const HGCRecHit* hit = itcheck->second;
                    tmprechits_.push_back(detid_);
                    tmprechitenergy_.push_back((it_sc_haf.second) * hit->energy());
                }
            }  // end of looping over the rechits
        }      // end of looping over the sim clusters

        tmpcp_.rechitdetid_ = tmprechits_;
        tmpcp_.rechitenergy_ = tmprechitenergy_;
        caloparticles.push_back(tmpcp_);

    }  // end of looping over the calo particles

    // get the relevant trackster collection
    // LG: For now always use the MergedTrackster
    //   int cp_pdgid_ = 0; if (caloparticles.size()>0) { cp_pdgid_ = caloparticles.at(0).pdgid_; }
    //   std::vector<ticl::Trackster> tracksters = getTracksterCollection(cp_pdgid_,emMCs, mipMCs, hadMCs, mergedMCs);

    // loop over the caloparticles and then find the closest trackster to it

    // keep tracksters [and the corresponding LC] that
    // have at least some ammount of energy from the CP
    //   std::vector<trackster> tracksterCollection; tracksterCollection.clear();

    float trackstersEnFromCP = 0.7; // expressed as fraction
    float maxDrTracksterCP = 0.3;
    auto const& tracksters = prod_tracksters;  // if we need a different trackster collection this should go in the loop

    for (unsigned int icp = 0; icp < caloparticles.size(); ++icp)
    {
        // get the relevant tracksterCollection based on the cp_id

        // find the tracksters within some DR from the CP
        std::vector<int> closestTracksters =
            getClosestTrackstersToCPByIndex(caloparticles[icp], tracksters, maxDrTracksterCP);

        // for those tracksters closest to the CP calculate the energy that is associated to the CP
        // and select the one with the closest Energy to the CP
        float tracksterCPEnDiffMin_ = std::numeric_limits<float>::max();
        int itracksterMin_ = -1;
        for ( int itrackster : closestTracksters )
        {
            float tracksterEnFromCP_ = getTracksterEnFromCP(tracksters[itrackster], lcs, caloparticles[icp]);
            float tracksterCPEnDiff = abs(caloparticles[icp].energy_ - tracksterEnFromCP_) / (caloparticles[icp].energy_);
            if (tracksterCPEnDiff < tracksterCPEnDiffMin_)
            {
                tracksterCPEnDiffMin_ = tracksterCPEnDiff;
                if ( tracksterCPEnDiff < trackstersEnFromCP )
                {
                    itracksterMin_ = itrackster;
                    cp_missingEnergyFraction = tracksterCPEnDiff;
                }
            }
        }

        if ( itracksterMin_ >= 0 )
        {
            irun = iEvent.id().run();
            ievent = iEvent.id().event();
            ilumiblock = iEvent.id().luminosityBlock();
            itime = iEvent.time();

            run = (size_t)irun;
            event = (size_t)ievent;
            lumi = (size_t)ilumiblock;

            lc_energy.clear();
            lc_eta.clear();
            lc_phi.clear();
            lc_layer.clear();

            std::vector<layercluster> lcsFromClosestTracksterToCP = getLCsFromTrackster(tracksters[itracksterMin_], lcs);
            std::sort(lcsFromClosestTracksterToCP.begin(), lcsFromClosestTracksterToCP.end(), sortLCsByEnergyAndLayer);

            ts_energy = tracksters[itracksterMin_].raw_energy();
            ts_sigma1 = tracksters[itracksterMin_].sigmasPCA()[0];
            ts_sigma2 = tracksters[itracksterMin_].sigmasPCA()[1];
            ts_sigma3 = tracksters[itracksterMin_].sigmasPCA()[2];
            cp_energy = caloparticles[icp].energy_;
            cp_eta = caloparticles[icp].eta_;
            cp_phi = caloparticles[icp].phi_;
            cp_pdgid = caloparticles[icp].pdgid_;
            trackster = itracksterMin_;

            for ( auto const& lc : lcsFromClosestTracksterToCP )
            {
                lc_energy.push_back(lc.energy_);
                lc_eta.push_back(abs(lc.eta_));
                lc_phi.push_back(lc.phi_);
                lc_layer.push_back(lc.layer_);
            }
            tree->Fill();     // Loop on trackster associated to the CaloParticle
        }
    }  // end of looping over the caloparticles

}

// ------------ method called once each job just before starting event loop  ------------
void HFNoseNtuplizer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void HFNoseNtuplizer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HFNoseNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFNoseNtuplizer);
