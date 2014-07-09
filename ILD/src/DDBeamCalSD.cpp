//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  A.Sailer, CERN
//  $Id:$
//====================================================================
#include "DDBeamCalSD.h"
#include "G4Step.hh"
//LCIO
#include "IMPL/SimCalorimeterHitImpl.h"

namespace DDSim {

  bool DDBeamCalSD::process(G4Step* aStep,G4TouchableHistory* /*history*/) {

    DD4hep::VolumeID myCellID = cellID(aStep);//inherited from Geant4Sensitive
#pragma warning "FIXME: Use bitshifting intead of pointer casting"
    int *cellIDs = (int*)&myCellID;

    //Treat energy deposit
    G4double eDep = aStep->GetTotalEnergyDeposit();

    //Create the hit
    IMPL::SimCalorimeterHitImpl* hit = new IMPL::SimCalorimeterHitImpl ;
    hit->setEnergy(eDep);
    hit->setCellID0(cellIDs[0]);
    hit->setCellID1(cellIDs[1]);
    G4ThreeVector lp = (aStep->GetPreStepPoint()->GetPosition()+aStep->GetPostStepPoint()->GetPosition())*0.5;
    float localposition[] = { (float)lp[0], (float)lp[1], (float)lp[2]};
    hit->setPosition(localposition);
    collection(m_collectionID)->add(hit);
    //Add at the hit to the collection

    return true;

  }

} // namespace DDSim

//##############################################################################################

namespace DD4hep{
  namespace Simulation{
    typedef DDSim::DDBeamCalSD DDBeamCalSD;
  }
}
#include "DDG4/Factories.h"

DECLARE_GEANT4SENSITIVE(DDBeamCalSD)
