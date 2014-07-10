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

namespace DDSim {

  bool SDCaloLCIO::process(G4Step* aStep,G4TouchableHistory* /*history*/) {

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
    std::cout << " CellID "<< std::hex << myCellID  << std::endl;
    std::cout << " CellID0 "<< std::hex << cellIDs[0]  << std::endl;
    std::cout << " CellID1 "<< std::hex << cellIDs[1]  << std::endl;
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
    typedef DDSim::SDCaloLCIO SDCaloLCIO;
  }
}
#include "DDG4/Factories.h"

DECLARE_GEANT4SENSITIVE(SDCaloLCIO)
