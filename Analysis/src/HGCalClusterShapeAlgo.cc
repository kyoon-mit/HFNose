// Modified version of RecoEcal/EgammaCoreTools/src/ClusterShapeAlgo.cc for HGCal

#include <iostream>

#include "HGCNose/Analysis/interface/HGCalClusterShapeAlgo.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
//#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

HGCalClusterShapeAlgo::HGCalClusterShapeAlgo ( const edm::ParameterSet& par ) : parameterSet_(par) {}


reco::ClusterShape HGCalClusterShapeAlgo::Calculate ( const reco::BasicCluster& passedCluster,
                                                      const HGCRecHitCollection* hits,
                                                      const CaloSubdetectorGeometry* geometry,
                                                      const CaloSubdetectorTopology* topology )
{
  Calculate_TopEnergy(passedCluster, hits);
  /*
  Calculate_2ndEnergy(passedCluster, hits);
  Create_Map(hits, topology);
  Calculate_e2x2();
  Calculate_e3x2();
  Calculate_e3x3();
  Calculate_e4x4();
  Calculate_e5x5();
  Calculate_e2x5Right();
  Calculate_e2x5Left();
  Calculate_e2x5Top();
  Calculate_e2x5Bottom();
  Calculate_Covariances(passedCluster, hits, geometry);
  Calculate_BarrelBasketEnergyFraction(passedCluster, hits, Eta, geometry);
  Calculate_BarrelBasketEnergyFraction(passedCluster, hits, Phi, geometry);
  Calculate_EnergyDepTopology(passedCluster, hits, geometry, true);
  Calculate_lat(passedCluster);
  Calculate_ComplexZernikeMoments(passedCluster);
  */

  return reco::ClusterShape ( covEtaEta_,
                              covEtaPhi_,
                              covPhiPhi_,
                              eMax_,
                              eMaxId_,
                              e2nd_,
                              e2ndId_,
                              e2x2_,
                              e3x2_,
                              e3x3_,
                              e4x4_,
                              e5x5_,
                              e2x5Right_,
                              e2x5Left_,
                              e2x5Top_,
                              e2x5Bottom_,
                              e3x2Ratio_,
                              lat_,
                              etaLat_,
                              phiLat_,
                              A20_,
                              A42_,
                              energyBasketFractionEta_,
                              energyBasketFractionPhi_ );
}


void HGCalClusterShapeAlgo::Calculate_TopEnergy ( const reco::BasicCluster& passedCluster, const HGCRecHitCollection* hits )
{ // Calculate the value of maximum E cluster and its detector Id
    double eMax = 0;
    DetId eMaxId(0);

    const std::vector<std::pair<DetId, float> >& clusterDetIds = passedCluster.hitsAndFractions();

    HGCRecHit testHGCRecHit;

    for (auto const& posCurrent : clusterDetIds)
    {
        if ( ( posCurrent.first != DetId(0) ) && ( hits->find(posCurrent.first) != hits->end() ) )
        {
            HGCRecHitCollection::const_iterator itt = hits->find(posCurrent.first);
            testHGCRecHit = *itt;
        }

        if ( testHGCRecHit.energy() * posCurrent.second > eMax )
        {
            eMax = testHGCRecHit.energy() * posCurrent.second;
            eMaxId = testHGCRecHit.id();
        }
    }

    eMax_ = eMax;
    eMaxId_ = eMaxId;
}


/*
void HGCalClusterShapeAlgo::Calculate_2ndEnergy(const reco::BasicCluster& passedCluster, const HGCRecHitCollection* hits) {
  double e2nd = 0;
  DetId e2ndId(0);

  const std::vector<std::pair<DetId, float> >& clusterDetIds = passedCluster.hitsAndFractions();

  EcalRecHit testHGCRecHit;

  for (auto const& posCurrent : clusterDetIds) {
    if ((posCurrent.first != DetId(0)) && (hits->find(posCurrent.first) != hits->end())) {
      HGCRecHitCollection::const_iterator itt = hits->find(posCurrent.first);
      testHGCRecHit = *itt;

      if (testHGCRecHit.energy() * posCurrent.second > e2nd && testHGCRecHit.id() != eMaxId_) {
        e2nd = testHGCRecHit.energy() * posCurrent.second;
        e2ndId = testHGCRecHit.id();
      }
    }
  }

  e2nd_ = e2nd;
  e2ndId_ = e2ndId;
}
*/
