//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//  17.07.2014
//  $Id:$
//====================================================================
#include "SDTrackerLCIO.h"
#include "G4Step.hh"
//LCIO
#include "IMPL/SimTrackerHitImpl.h"

namespace DDSim {

  bool SDTrackerLCIO::process(G4Step* aStep,G4TouchableHistory* /*history*/) {

    DD4hep::VolumeID theCellID = cellID(aStep);//inherited from Geant4Sensitive

    //Treat energy deposit
    G4double eDep = aStep->GetTotalEnergyDeposit();

    //Create the hit
    IMPL::SimTrackerHitImpl* hit = new IMPL::SimTrackerHitImpl ;
    hit->setEDep(eDep);

    hit->setCellID0(   theCellID         & 0xffffffff );
    hit->setCellID1( ( theCellID >> 32 ) & 0xffffffff );

    G4ThreeVector lp = (aStep->GetPreStepPoint()->GetPosition()+aStep->GetPostStepPoint()->GetPosition())*0.5;

    double localposition[] = { lp[0], lp[1], lp[2] };
    hit->setPosition(localposition);
    collection(m_collectionID)->add(hit);
    //Add at the hit to the collection

    return true;

  }

} // namespace DDSim

//##############################################################################################

namespace DD4hep{
  namespace Simulation{
    typedef DDSim::SDTrackerLCIO SDTrackerLCIO;
  }
}
#include "DDG4/Factories.h"

DECLARE_GEANT4SENSITIVE(SDTrackerLCIO)
