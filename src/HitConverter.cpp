// Framework include files
#define DDG4_MAKE_INSTANTIATIONS
#include "DD4hep/LCDD.h"
#include "DDG4/Geant4HitCollection.h"
#include "DDG4/Geant4DataConversion.h"
#include "DDG4/Geant4Data.h"

// LCIO includes
#include "lcio.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/MCParticleImpl.h"
#include "IMPL/SimTrackerHitImpl.h"
#include "IMPL/SimCalorimeterHitImpl.h"
#include "UTIL/Operators.h"
#include "UTIL/ILDConf.h"


namespace DD4hep {
  namespace Simulation   {
    

    typedef Geometry::VolumeManager VolMgr;
    
    //==================================================================================
    // copied from DD4hep LCIOConversions.cpp

    template <> lcio::LCCollectionVec*
    Geant4DataConversion<lcio::LCCollectionVec,std::pair<VolMgr,G4VHitsCollection*>,Geant4HitCollection>::

    operator()(const arg_t& args)  const {
      G4VHitsCollection* c = args.second;
      Geant4HitCollection* coll = dynamic_cast<Geant4HitCollection*>(c);
      if ( coll )  {
	typedef std::pair<arg_t::first_type,Geant4HitCollection*> _A;
	typedef Geant4Conversion<output_t,_A> _C;
	const _C& cnv= _C::converter(coll->type().type);
	return cnv(_A(args.first,coll));
      }
      throw unrelated_type_error(typeid(Geant4HitCollection),typeid(*c),
				 "Cannot save the collection entries of:"+c->GetName());
    }

    //================== convert SimTrackerHits ========================================
    template <> lcio::LCCollectionVec*
    Geant4DataConversion<lcio::LCCollectionVec,std::pair<VolMgr,Geant4HitCollection*>,lcio::SimTrackerHitImpl>::operator()(const arg_t& args)  const
    {

      Geant4HitCollection* g4col = args.second;
      size_t     nhits = g4col->GetSize();
      
      lcio::LCCollectionVec* col = new lcio::LCCollectionVec(lcio::LCIO::SIMTRACKERHIT);
      col->reserve( nhits ) ;

      // need to write the description to the collection - look it up via the cellID of the first hit
      lcio::SimTrackerHit* hit   = g4col->hit(0) ;
      long cellID = (  ( ( hit->getCellID1() << 32 ) & 0xffffffff00000000 )  |  ( hit->getCellID0() & 0xffffffff ) ) ; 
      std::string     dsc   = encoding(  /*VolumeManager*/ args.first , cellID );
      UTIL::CellIDEncoder<lcio::SimTrackerHit> decoder(dsc,col);
      
      // simply copy all hits 
      for(size_t i=0; i<nhits; ++i)   {
	col->addElement(   g4col->hit(i)   );
      }
      
      return col;
    }

    //================== convert SimCalorimeterHits ========================================
    template <> lcio::LCCollectionVec*
    Geant4DataConversion<lcio::LCCollectionVec,std::pair<VolMgr,Geant4HitCollection*>,lcio::SimCalorimeterHitImpl>::operator()(const arg_t& args)  const
    {

      Geant4HitCollection* g4col = args.second;
      size_t     nhits = g4col->GetSize();
      
      lcio::LCCollectionVec* col = new lcio::LCCollectionVec(lcio::LCIO::SIMCALORIMETERHIT);
      col->reserve( nhits ) ;

      // need to write the description to the collection - look it up via the cellID of the first hit
      lcio::SimCalorimeterHit* hit   = g4col->hit(0) ;
      long cellID = (  ( ( hit->getCellID1() << 32 ) & 0xffffffff00000000 )  |  ( hit->getCellID0() & 0xffffffff ) ) ; 
      std::string     dsc   = encoding(  /*VolMgr*/ args.first , cellID );
      UTIL::CellIDEncoder<lcio::SimCalorimeterHit> decoder(dsc,col);
      
      // simply copy all hits 
      for(size_t i=0; i<nhits; ++i)   {
	col->addElement(   g4col->hit(i)   );
      }
      
      return col;
    }

    //================== convert MCParticles ========================================
    // can MCParticles be treated as G4VHits ?
    // ...probably not ...

    template <> lcio::LCCollectionVec*
    Geant4DataConversion<lcio::LCCollectionVec,std::pair<VolMgr,Geant4HitCollection*>,lcio::MCParticleImpl>::operator()(const arg_t& args)  const
    {

      Geant4HitCollection* g4col = args.second;
      size_t     nhits = g4col->GetSize();
      
      lcio::LCCollectionVec* col = new lcio::LCCollectionVec(lcio::LCIO::MCPARTICLE);
      col->reserve( nhits ) ;

      // simply copy all hits 
      for(size_t i=0; i<nhits; ++i)   {
	col->addElement(   g4col->hit(i)   );
      }
      
      return col;
    }
    //========================================================================================
  

    typedef std::pair<VolMgr,G4VHitsCollection*> _AA1;
    template class Geant4Conversion<lcio::LCCollectionVec,_AA1>;
    
    DECLARE_GEANT4_HITCONVERTER(lcio::LCCollectionVec,_AA1,Geant4HitCollection)
    
    typedef std::pair<VolMgr,Geant4HitCollection*> _AA2;
    template class Geant4Conversion<lcio::LCCollectionVec,_AA2>;
    
    DECLARE_GEANT4_HITCONVERTER(lcio::LCCollectionVec,_AA2, lcio::SimTrackerHitImpl)

    typedef std::pair<VolMgr,Geant4HitCollection*> _AA3;
    DECLARE_GEANT4_HITCONVERTER(lcio::LCCollectionVec,_AA3, lcio::SimCalorimeterHitImpl)

    typedef std::pair<VolMgr,Geant4HitCollection*> _AA4;
    DECLARE_GEANT4_HITCONVERTER(lcio::LCCollectionVec,_AA3, lcio::MCParticleImpl)



  }    // End namespace Simulation
}      // End namespace DD4hep
