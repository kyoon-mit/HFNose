#include "EnergyResolution.h"

#include <typeinfo>

#include <cmath> // Switch to TMath.h if you need more physics-related functions
#include "DataFormats/Math/interface/LorentzVector.h"
// #include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// test #include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "FWCore/MessageLogger/interface/MessageLogger.h" ??

// HFNose (forward + HGCal)
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"

// HF (forward + HCal)
// #include "DataFormats/HCalRecHit/interface/HCalRecHitCollections.h"
// #include "DataFormats/HCalDetId/interface/HCalDetId.h"

// CMS Coordinate System
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

// Detector Geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

// ROOT headers

using namespace edm;

EnergyResolution::EnergyResolution ( const edm::ParameterSet& iConfig ) :

    histContainer_ (),
    
    // (tag name, default value (label, instance, process) -- CHECK SPELLING!!!!!!!
    tag_HGCHFNoseRecHits_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_HGCHFNoseRecHits", edm::InputTag ("HGCalRecHit", "HGCHFNoseRecHits") ) ),
    tag_genParticles_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_genParticles", edm::InputTag ("genParticles") ) )
    // tag_CaloParticle_MergedCaloTruth_ ( iConfig.getUntrackedParameter<edm::InputTag> ("TAG_MergedCaloTruth", edm::InputTag ("mix", "MergedCaloTruth") ) ), 
    
    // Pre-selection parameters
    select_PID_ ( 22 ),
    select_EtaLow_ ( 3.49 ),
    select_EtaHigh_ ( 3.51 ),
    select_coneR_ ( 0.5 )
    
{
    // consumes: frequent request of additional data | mayConsume: infrequent
    token_HGCRecHits_ = consumes<HGCRecHitCollection> ( tag_HGCHFNoseRecHits_ );
    token_genParticles_ = consumes<genParticles> ( tag_genParticles_ );
    // token_CaloParticle_MergedCaloTruth_ = mayConsume<std::vector<CaloParticle>> ( tag_CaloParticle_MergedCaloTruth_ );
}


EnergyResolution::~EnergyResolution ()
{
    // Deconstructor
}


void EnergyResolution::analyze ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
    // Get hits
    edm::Handle<HGCRecHitCollection> handle_HGCRecHits;
    iEvent.getByToken ( token_HGCRecHits_, handle_HGCRecHits );

    // Get Geometry
    edm::ESHandle<HGCalGeometry> handle_HGCalGeometry;
    iSetup.get<IdealGeometryRecord>().get( "HGCalHFNoseSensitive", handle_HGCalGeometry );
    const HGCalGeometry* geom_HGCal_ = handle_HGCalGeometry.product();

    // Get MC truth (maybe consider using a separate function if code becomes too long)
    // edm::Handle<std::vector<CaloParticle>> handle_CaloParticle_MergedCaloTruth_;
    // iEvent.getByToken ( token_CaloParticle_MergedCaloTruth_, handle_CaloParticle_MergedCaloTruth_ );
    edm::Handle<std::vector<genParticles>> handle_genParticles_;
    iEvent.getByToken ( token_genParticles_, handle_genParticles_ );
    
    // Int_t caloTruth_count_validation = 0; // Must be 1 in this event
    // math::XYZTLorentzVectorF this_truth_p4;
    
    if ( handle_CaloParticle_MergedCaloTruth_.isValid() )
    {
        for ( auto const& caloTruth : *handle_CaloParticle_MergedCaloTruth_.product() )
        {   
            if ( caloTruth.pdgId() == select_PID_ )
            { // Nested to save a bit of calculation
                if ( ( select_EtaLow_ < std::abs(caloTruth.eta()) ) &&
                     ( std::abs(caloTruth.eta()) < select_EtaHigh_ ) )
                {
                    caloTruth_count_validation++;
                    if ( caloTruth_count_validation > 1 )
                    {
                        std::cout << "CaloTruth count 2; skipping event" << std::endl;
                        return;
                    }
                    
                    this_truth_p4 = caloTruth.p4();
                    
                }
            }
        }
        
    }
    else std::cout << "Handle for CaloParticle:MergeCaloTruth invalid!" << std::endl;
    
    if ( caloTruth_count_validation == 0 )
    {
        std::cout << "CaloTruth count 0; skipping event" << std::endl;
        return;
    }
    
    if ( handle_HGCRecHits.isValid() )
    {
        Float_t sum_E = 0;
        Float_t truth_E = this_truth_p4.energy();
        
        for ( auto const& hit : *handle_HGCRecHits.product() )
        {
            Float_t hit_energy = hit.energy();
            uint32_t hit_id = hit.id();
 //           HFNoseDetId hit_DetId = HFNoseDetId ( hit_id );
//            Int_t hit_DetLayer = hit_DetId.layer();
            
            // std::shared_ptr<const CaloCellGeometry> thisCell = geom_HGCal_->getGeometry(hit_DetId);
            const GlobalPoint & hit_globalPosition = geom_HGCal_->getPosition(hit_id);
            Float_t dR = reco::deltaR ( hit_globalPosition.eta(), hit_globalPosition.phi(), this_truth_p4.eta(), this_truth_p4.phi() );
            
            if ( dR < select_coneR_ ) // hit within cone of truth
            {
                sum_E += hit_energy;
            }
        }
        
//        std::cout << sum_E << std::endl;
        histContainer_["EDist"]->Fill ( sum_E );
        histContainer_["truthEDist"]->Fill ( truth_E );
    }
    else std::cout << "Handle for HGCRecHits invalid!" << std::endl;
}


void EnergyResolution::beginJob ()
{
    edm::Service<TFileService> fs;
    
    histContainer_["EDist"] = fs->make<TH1F>("EDist", "Energy Distribution", 200, 0, 200);
    histContainer_["truthEDist"] = fs->make<TH1F>("truthEDist", "Truth Energy Distribution", 200, 0, 200);
}


void EnergyResolution::endJob ()
{
    // Per Job (or can combine multiple files in a separate module)
    // Bin size doesn't seem to affect
    // Alt, you can use RooFit (seems like the best option to me) <- but in separate script
    const Float_t det_E_Mean        = histContainer_["EDist"]->GetMean();
    const Float_t det_E_StdDev      = histContainer_["EDist"]->GetStdDev();
    const Float_t truth_E_Mean      = histContainer_["truthEDist"]->GetMean();
    const Float_t truth_E_StdDev    = histContainer_["truthEDist"]->GetStdDev();
    const Float_t E_Res             = det_E_StdDev / det_E_Mean;
    
    // Control
    std::cout << "Truth mean, std dev: " << truth_E_Mean << ", " << truth_E_StdDev << std::endl;
    std::cout << "Measured mean, std dev: " << det_E_Mean << ", " << det_E_StdDev << std::endl;
    std::cout << "ENERGY RESOLUTION: " << E_Res * 100. << " %" << std::endl;

    // Job ends
    std::cout << "The job, EnergyResolution, has ended. Thank you for your patience." << std::endl;
}

/*
void EnergyResolution::SelectTruth ()
{   // Select truth with data member select_*
    // Not sure how to go about this
}


Float_t EnergyResolution::SumHitEnergy_TruthConeR ()
{   // Sum the energy of all hits within Cone R of the truth
    // Do this if jet reco algo needed
}
*/

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (EnergyResolution);
