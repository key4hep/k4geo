#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Geant4EventAction.h"
#include "DDG4/Geant4Mapping.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {
  
  /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
  namespace sim   {
    
    /**
     *  Geant4SensitiveAction<CalorimeterWithPreShowerLayer> sensitive detector for the special
     *  case of a calorimeter that has a pre-shower layer, i.e. one sensitive layer before
     *  the first absorber layer. This is for example used in the ILD Ecal. 
     *  Hits from the first layer are stored in a separate collection named READOUT_NAME_preShower.
     *
     *  \author  F.Gaede
     *  \version 1.0
     *  \ingroup DD4HEP_SIMULATION
     */
    struct CalorimeterWithPreShowerLayer: public Geant4Calorimeter{
      G4int _preShowerCollectionID ;
      G4int _firstLayerNumber ; 
      Geant4HitCollection *_preShowerCollection;
      CalorimeterWithPreShowerLayer() : Geant4Calorimeter(), 
					_preShowerCollectionID(0),
					_firstLayerNumber(1), //fixme: can we make this a parameter ?
					_preShowerCollection(0)
      {}
    };




   /// Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<CalorimeterWithPreShowerLayer>::defineCollections() {
      std::string name = m_sensitive.readout().name() ;
      m_collectionID = defineCollection<CalorimeterWithPreShowerLayer::Hit>(name);
      name += std::string("_preShower") ;
      m_userData._preShowerCollectionID = defineCollection<CalorimeterWithPreShowerLayer::Hit>( name );
    }


    /// template specialization for c'tor in order to define property: FirstLayerNumber
    template <> 
    Geant4SensitiveAction<CalorimeterWithPreShowerLayer>::Geant4SensitiveAction(Geant4Context* ctxt,
										const std::string& nam,
										DetElement det,
										Detector& lcdd_ref)
      : Geant4Sensitive(ctxt,nam,det,lcdd_ref), m_collectionID(0)
    {
      initialize();
      defineCollections();
      InstanceCount::increment(this);
      declareProperty("FirstLayerNumber", m_userData._firstLayerNumber = 1 );
    }

    /// Method for generating hit(s) using the information of G4Step object.
    template <> bool Geant4SensitiveAction<CalorimeterWithPreShowerLayer>::process(G4Step* step,G4TouchableHistory*) {
      typedef CalorimeterWithPreShowerLayer::Hit Hit;
      Geant4StepHandler h(step);
      HitContribution contrib = Hit::extractContribution(step);

      long long int cell;
      try {
        cell = cellID(step);
      } catch(std::runtime_error &e) {
        std::stringstream out;
        out << std::setprecision(20) << std::scientific;
        out << "ERROR: " << e.what()  << std::endl;
        out << "Position: "
            << "Pre (" << std::setw(24) << step->GetPreStepPoint()->GetPosition() << ") "
            << "Post (" << std::setw(24) << step->GetPostStepPoint()->GetPosition() << ") "
            << std::endl;
        out << "Momentum: "
            << " Pre (" <<std::setw(24) << step->GetPreStepPoint() ->GetMomentum()  << ") "
            << " Post (" <<std::setw(24) << step->GetPostStepPoint()->GetMomentum() << ") "
            << std::endl;

        std::cout << out.str();

        return true;
      }

      // get the layer number by decoding the cellID
      IDDescriptor idspec = m_sensitive.readout().idSpec() ;
      const DDSegmentation::BitFieldCoder& bc = *idspec.decoder() ;
      int layer = bc.get(cell, "layer") ;
      
      Geant4HitCollection*  coll = ( layer== m_userData._firstLayerNumber ?  collection( m_userData._preShowerCollectionID ) : collection(m_collectionID) ) ;
      
      Hit* hit = coll->find<Hit>(CellIDCompare<Hit>(cell));
      if ( h.totalEnergy() < std::numeric_limits<double>::epsilon() )  {
        return true;
      }
      else if ( !hit ) {
        Geant4TouchableHandler handler(step);
        DDSegmentation::Vector3D pos = m_segmentation.position(cell);
        Position global = h.localToGlobal(pos);
        hit = new Hit(global);
        hit->cellID = cell;
        coll->add(hit);
        printM2("%s> CREATE hit with deposit:%e MeV  Pos:%8.2f %8.2f %8.2f  %s",
                c_name(),contrib.deposit,pos.X,pos.Y,pos.Z,handler.path().c_str());
        if ( 0 == hit->cellID )  { // for debugging only!
          hit->cellID = cellID(step);
          except("+++ Invalid CELL ID for hit!");
        }
      }
      hit->truth.push_back(contrib);
      hit->energyDeposit += contrib.deposit;
      mark(step);
      return true;
    }



    typedef Geant4SensitiveAction<CalorimeterWithPreShowerLayer> CaloPreShowerSDAction;

  } // namespace
} // namespace



#include "DDG4/Factories.h"
DECLARE_GEANT4SENSITIVE( CaloPreShowerSDAction )
