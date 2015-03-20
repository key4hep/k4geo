#ifndef envelope__h
#define envelope__h

#include "DD4hep/DetFactoryHelper.h"


/** Create an envelope volume that is placed into the world volume.
 *  
 *  
 *  @author S.Lu DESY, F. Gaede CERN/DESY 
 *  @version $Id:$
 */
DD4hep::Geometry::Ref_t create_placed_envelope(DD4hep::Geometry::LCDD& lcdd, DD4hep::XML::Handle_t e , 
					       DD4hep::Geometry::DetElement sdet  ) ;


#endif
