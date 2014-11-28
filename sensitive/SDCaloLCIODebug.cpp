//====================================================================
//  lcgeo - LC simulation based on DD4hep
//--------------------------------------------------------------------
//  A.Sailer, CERN
//  $Id$
//====================================================================
#include "SDCaloLCIODebug.h"
#include "G4Step.hh"
//LCIO
#include "IMPL/SimCalorimeterHitImpl.h"
#include "DD4hep/Segmentations.h"

#include <G4TouchableHistory.hh>
#include <G4AffineTransform.hh>



namespace lcgeo {

  inline G4ThreeVector localToGlobalCoordinates(G4Step *aStep, G4ThreeVector const& localPosition) {

    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    G4ThreeVector worldPosition = theTouchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(localPosition);

    return worldPosition;

  }


  bool SDCaloLCIODebug::process(G4Step* aStep, G4TouchableHistory* history) {

    SDCaloLCIO::process(aStep, history);

    DD4hep::VolumeID myCellID = cellID(aStep);//inherited from Geant4Sensitive
    G4double eDep = aStep->GetTotalEnergyDeposit();

    const DD4hep::DDSegmentation::Vector3D& posVec = m_readout.segmentation()->position(myCellID);
    const G4ThreeVector localCell(CM_2_MM*posVec.X, CM_2_MM*posVec.Y, CM_2_MM*posVec.Z);
    const G4ThreeVector& globalCell = localToGlobalCoordinates(aStep, localCell);
    G4ThreeVector hitPos = (aStep->GetPreStepPoint()->GetPosition()+aStep->GetPostStepPoint()->GetPosition())*0.5;

    hDeltaX->Fill(globalCell[0] - hitPos[0]);
    hDeltaY->Fill(globalCell[1] - hitPos[1]);
    hZ -> Fill (globalCell[2]);
    hEdep->Fill(localCell[0], localCell[1], eDep);

    G4cout << std::setw(13) << "stepPosition"
	   << std::setw(13) << hitPos[0]
	   << std::setw(13) << hitPos[1]
	   << std::setw(13) << hitPos[2]
	   << std::endl;
    G4cout << std::setw(13) << "cellGlobal"
	   << std::setw(13) << globalCell[0]
	   << std::setw(13) << globalCell[1]
	   << std::setw(13) << globalCell[2]
	   << std::endl;
    G4cout << std::setw(13) << "diff"
	   << std::setw(13) << globalCell[0] - hitPos[0]
	   << std::setw(13) << globalCell[1] - hitPos[1]
	   << std::setw(13) << globalCell[2] - hitPos[2]
	   << std::endl;
    G4cout << std::setw(13) << "cellLocal"
	   << std::setw(13) << localCell[0]
	   << std::setw(13) << localCell[1]
	   << std::setw(13) << localCell[2]
	   << std::endl;



    return true;

  }

} // namespace lcgeo

//##############################################################################################

namespace DD4hep{
  namespace Simulation{
    typedef lcgeo::SDCaloLCIODebug SDCaloLCIODebug;
  }
}
#include "DDG4/Factories.h"

DECLARE_GEANT4SENSITIVE(SDCaloLCIODebug)
