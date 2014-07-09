//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  A.Sailer, CERN
//  $Id:$
//====================================================================
#include "DDBeamCalSD.h"

#include "DDG4/Geant4StepHandler.h"
#include "DDG4/Geant4VolumeManager.h"
#include "DDG4/Geant4Mapping.h"

#include "DDSegmentation/BitField64.h"


//Geant4
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4TouchableHistory.hh"
#include "G4AffineTransform.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "globals.hh"

//LCIO
#include "IMPL/SimCalorimeterHitImpl.h"


#include <cassert>


namespace DDSim {

  bool DDBeamCalSD::process(G4Step* aStep,G4TouchableHistory* /*history*/) {

    DD4hep::VolumeID myCellID = cellID(aStep);//inherited from Geant4Sensitive
    int *cellIDs = (int*)&myCellID;
    //Treat energy deposit
    G4double eDep = aStep->GetTotalEnergyDeposit();
    //Add at the hit to the collection
    IMPL::SimCalorimeterHitImpl* hit = new IMPL::SimCalorimeterHitImpl ;
    hit->setEnergy(eDep);
    hit->setCellID0(cellIDs[0]);
    hit->setCellID1(cellIDs[1]);
    collection(m_collectionID)->add(hit);

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
