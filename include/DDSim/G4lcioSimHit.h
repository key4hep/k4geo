//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//====================================================================
#ifndef G4lcioSimHit_h
#define G4lcioSimHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

#include "lcio.h"
#include "IMPL/SimTrackerHitImpl.h"
#include "IMPL/SimCalorimeterHitImpl.h"
#include "UTIL/Operators.h"

#include <iostream>

namespace DDSim {

  using namespace UTIL ;

  /** Template wrapper class for lcio SimTrackerHit/SimCalorimeterHits.
   *
   *  @author: F.Gaede, DESY, Nov 2013
   *  @version: $Id:$
   */
  template <class LSH>
  class G4lcioSimHit : public G4VHit{
  
    LSH* _lsh ;
    bool _owner ;

  public:
  
    G4lcioSimHit() : _lsh( new LSH ), _owner( true ) {}
  
    ~G4lcioSimHit() { if( _owner ) delete _lsh ; }
  
    /** The underlying lcio SimHit */
    LSH* lcio_hit() { return _lsh ; } 
  
    /** Return the lcio SimHit and take ownership away from this object
     *  as is needed when passing the hit to the lcio collection/event */
    LSH* take_lcio_hit() { 
      _owner = false ; 
      return _lsh ; 
    } 
  
    inline void *operator new(size_t){

      void *aHit;

      extern G4Allocator< G4lcioSimHit<LSH> > SimHitAllocator ;
    
      aHit = (void *) SimHitAllocator.MallocSingle();
    
      return aHit;
    }
  

    inline void operator delete(void *aHit){
    
      extern G4Allocator< G4lcioSimHit<LSH> > SimHitAllocator ;
    
      SimHitAllocator.FreeSingle(( G4lcioSimHit<LSH>*) aHit);
    }

    void Draw() {  // to be done ... 
    }

    void Print() {  std::cout << *_lsh << std::endl ;  }
  
  } ;

  typedef G4lcioSimHit<lcio::SimTrackerHitImpl> SimTrkHit ;
  typedef G4THitsCollection<SimTrkHit> SimTrkHitsCollection;

  typedef G4lcioSimHit<lcio::SimCalorimeterHitImpl> SimCaloHit ;
  typedef G4THitsCollection<SimCaloHit> SimCaloHitsCollection;



  // extern G4Allocator<SimTrkHit> SimTrkHitAllocator;
  // template<>
  // inline void* SimTrkHit::operator new(size_t)
  // {
  //   void *aHit;
  //   aHit = (void *) SimTrkHitAllocator.MallocSingle();
  //   return aHit;
  // }
  // template<>
  // inline void SimTrkHit::operator delete(void *aHit)
  // {
  //   SimTrkHitAllocator.FreeSingle((SimTrkHit*) aHit);
  // }


  //  extern G4Allocator<SimCaloHit> SimCaloHitAllocator;
  // template<>
  // inline void* SimCaloHit::operator new(size_t)
  // {
  //   void *aHit;
  //   aHit = (void *) SimCaloHitAllocator.MallocSingle();
  //   return aHit;
  // }
  // template<>
  // inline void SimCaloHit::operator delete(void *aHit)
  // {
  //   SimCaloHitAllocator.FreeSingle((SimCaloHit*) aHit);
  // }

} //namespace
#endif
