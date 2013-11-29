//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//====================================================================
#ifndef TRKSiSD00_h
#define TRKSiSD00_h 1
#include "TRKHit.hh"

#include "VSensitiveDetector.hh"
#include "G4Step.hh"

/** Sensitive detector for silicon trackers - ported from Mokka TRKSiSD00  
 *  (http://polzope.in2p3.fr:8081/MOKKA)
 *
 *  @author: F.Gaede, DESY, Nov 2013
 *  @version: $Id:$
 */
 class TRKSiSD00 : public VSensitiveDetector
{
  
public:
  TRKSiSD00(G4String TRKSiSD00name, 
	  G4double Threshold,
	  G4double thePrimaryTPCCut=0);
  virtual ~TRKSiSD00(){}
  
  void Initialize(G4HCofThisEvent*HCE);
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*HCE);
  void DrawAll();
  void PrintAll();

  // Basic reload hits from Ascii Files
  void LoadEvent(FILE* theSubDetectorEventHitsFileInput);

private:
  void StartNewHit(G4int aLayerNumber,
		   G4ThreeVector theEntryPoint,
		   G4ThreeVector theEntryMomentum,
		   G4int thePDG,
		   G4int theTrackPID);
  
  void UpdateHit(G4Step *aStep);

  void DumpHit(G4Step *aStep);

  void Clear();

  G4double Threshold;
  G4double PrimaryTPCCut;
  G4int HCID;
  G4int currentCylinder;
  G4int currentPID;
  G4int currentPDG;

  G4ThreeVector EntryPoint;
  G4ThreeVector ExitPoint;
  G4ThreeVector EntryMomentum;
  G4ThreeVector ExitMomentum;
  G4double DepositedEnergy;
  G4double HitTime ;
  G4double StepLength;

  G4bool detailedHitsStoring;
  G4bool defaultHitsStoring;

  TRKHitsCollection *CalCollection;

private:
};

#endif

