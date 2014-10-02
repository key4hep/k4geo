//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  A.Sailer, CERN
//====================================================================
#ifndef SDCaloLCIODebug_h
#define SDCaloLCIODebug_h 1


#include "SDCaloLCIO.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>


namespace DDSim {

  class SDCaloLCIODebug :  public SDCaloLCIO {
  
  protected:
    TFile *m_file;
    TH1D *hDeltaX, *hDeltaY;
    TH1D *hZ;
    TH2D *hEdep;
  public:
    //======================================================================================================
  
    SDCaloLCIODebug(dd4hep::Geant4Context* mContext, const std::string& mName, DetElement mDet, LCDD& mLcdd):
      SDCaloLCIO(mContext, mName, mDet, mLcdd),
      m_file(TFile::Open( TString(mName)+TString(".root") , "RECREATE" )),
      hDeltaX(new TH1D("hDeltaX", "hDeltaX;x[mm]", 6000, -30, 30)),
      hDeltaY(new TH1D("hDeltaY", "hDeltaY;y[mm]", 6000, -30, 30)),
      hZ(new TH1D("hZ", "hZ;z[mm]", 20000, -4000, 4000)),
      hEdep(new TH2D("hEdep", "Edep[MeV];x[mm];y[mm]", 400, -20, 20, 400, -20, 20))
    {
      defineCollections();
      dd4hep::InstanceCount::increment(this);
    }
  
    /// Default destructor
    virtual ~SDCaloLCIODebug(){
      dd4hep::InstanceCount::decrement(this);
      m_file->Write();
      m_file->Close();
    }
  
    /// Define collections created by this sensitive action object
    virtual void defineCollections() {
      
      m_collectionID = dd4hep::Geant4Sensitive::defineCollection<IMPL::SimCalorimeterHitImpl>( name() );
  }

    /// G4VSensitiveDetector interface: Method invoked at the begining of each event.
    virtual void begin(G4HCofThisEvent* hce) {
      dd4hep::Geant4Sensitive::begin(hce);
    }
  
    /// G4VSensitiveDetector interface: Method invoked at the end of each event.
    virtual void end(G4HCofThisEvent* hce) {
      dd4hep::Geant4Sensitive::end(hce);
      SDCaloLCIO::end(hce);
    }
  
    /// G4VSensitiveDetector interface: Method for generating hit(s) using the G4Step object.
    virtual bool process(G4Step* step,G4TouchableHistory* history) ;
  
    /// G4VSensitiveDetector interface: Method invoked if the event was aborted.
    virtual void clear(G4HCofThisEvent* hce) {
      dd4hep::Geant4Sensitive::clear(hce);
      SDCaloLCIO::clear(hce);
    }
  
    //======================================================================================================

  };

} // namespace
#endif

