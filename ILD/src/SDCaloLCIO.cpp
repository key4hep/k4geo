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

#include <G4TouchableHistory.hh>
#include <G4AffineTransform.hh>



namespace DDSim {

  inline G4ThreeVector localToGlobalCoordinates(G4Step *aStep, G4ThreeVector& localPosition) {

    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    G4ThreeVector worldPosition = theTouchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(localPosition);
  
    return worldPosition;

  }


  bool SDCaloLCIO::process(G4Step* aStep, G4TouchableHistory* /*history*/) {

    DD4hep::VolumeID myCellID = cellID(aStep);//inherited from Geant4Sensitive

    //Treat energy deposit
    G4double eDep = aStep->GetTotalEnergyDeposit();
    if (not (eDep > 0)) return true;

    //FIXME Check if the cell already has a hit and then add to it, not always create a new hit
    //Create the hit
    IMPL::SimCalorimeterHitImpl* hit = new IMPL::SimCalorimeterHitImpl ;
    hit->setEnergy(eDep);

    const int cellID0a = ( myCellID >>  0 ) & 0x00000000ffffffff;
    const int cellID1a = ( myCellID >> 32 ) & 0x00000000ffffffff;
    hit->setCellID0( cellID0a );
    hit->setCellID1( cellID1a );

    const DD4hep::DDSegmentation::Vector3D& posVec = m_readout.segmentation()->position(myCellID);
    G4ThreeVector localCell(posVec.X, posVec.Y, posVec.Z);
    const G4ThreeVector& globalCell = localToGlobalCoordinates(aStep, localCell);

    float worldPosition[] = { (float)globalCell[0], (float)globalCell[1], (float)globalCell[2]};
    hit->setPosition(worldPosition);
    //Add at the hit to the collection
    collection(m_collectionID)->add(hit);
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
