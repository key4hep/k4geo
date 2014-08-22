//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  A.Sailer, CERN
//  $Id:$
//====================================================================
#include "SDCaloLCIO.h"
#include "G4Step.hh"
//LCIO
#include "IMPL/SimCalorimeterHitImpl.h"
#include "DD4hep/Segmentations.h"
#include "DD4hep/DD4hepUnits.h"

#include <G4TouchableHistory.hh>
#include <G4AffineTransform.hh>



namespace DDSim {

  inline G4ThreeVector localToGlobalCoordinates(G4Step *aStep, G4ThreeVector const& localPosition) {

    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    G4ThreeVector worldPosition = theTouchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(localPosition);
  
    return worldPosition;

  }


  bool SDCaloLCIO::process(G4Step* aStep, G4TouchableHistory* /*history*/) {

    const DD4hep::VolumeID myCellID = cellID(aStep);//inherited from Geant4Sensitive

    //Treat energy deposit
    G4double eDep = aStep->GetTotalEnergyDeposit();
    if (not (eDep > 0)) return true;

    //Check if the hit already exists
    MapCellHit::iterator existingHit = m_mapHits.find(myCellID);
    if ( existingHit != m_mapHits.end() ) {
      IMPL::SimCalorimeterHitImpl* hit = existingHit->second; 
      //add new contribution to the hit
      addContribution(aStep, hit);
    } else {
      createNewHit(aStep, myCellID);
    }

    return true;

  }

  inline void SDCaloLCIO::addContribution(G4Step *aStep, IMPL::SimCalorimeterHitImpl *hit) {
    const G4int pdg = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    const G4double time = aStep->GetTrack()->GetGlobalTime();
    const G4double eDep = aStep->GetTotalEnergyDeposit();
    //NULL for MCParticle, because I don't know which it will be
    hit->addMCParticleContribution(NULL, eDep, time, pdg);
  }


  inline bool SDCaloLCIO::createNewHit(G4Step *aStep, DD4hep::CellID myCellID){

    //Create the hit
    IMPL::SimCalorimeterHitImpl* hit = new IMPL::SimCalorimeterHitImpl ;

    //calculate cellID0 and 1
    const int cellID0 = ( myCellID >>  0 ) & 0x00000000ffffffff;
    const int cellID1 = ( myCellID >> 32 ) & 0x00000000ffffffff;
    hit->setCellID0( cellID0 );
    hit->setCellID1( cellID1 );

    //calculate global cell position
    const DD4hep::DDSegmentation::Vector3D& posVec = m_readout.segmentation()->position(myCellID);
    const G4ThreeVector localCellPosition( posVec.X / dd4hep::mm  , posVec.Y / dd4hep::mm , posVec.Z / dd4hep::mm );
    const G4ThreeVector& globalCellPosition = localToGlobalCoordinates(aStep, localCellPosition);
    float globalPosition[] = { (float)globalCellPosition[0], (float)globalCellPosition[1], (float)globalCellPosition[2]};
    hit->setPosition(globalPosition);

    //add energy deposit from this hit/particle
    addContribution(aStep, hit);

    //Add the hit to the collection
    collection(m_collectionID)->add(hit);

    //addToLocalCollection
    m_mapHits[myCellID] = hit;

    return true;

  }

} // namespace DDSim

//##############################################################################################

namespace DD4hep{
  namespace Simulation{
    typedef DDSim::SDCaloLCIO SDCaloLCIO;
  }
}
#include "DDG4/Factories.h"

DECLARE_GEANT4SENSITIVE(SDCaloLCIO)
