#include "detectorSegmentations/DRparamBarrel_k4geo.h"
#include "detectorSegmentations/DRparamEndcap_k4geo.h"
#include "detectorSegmentations/GridDRcaloHandle_k4geo.h"

#include "DRconstructor.h"

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Detector.h"

namespace ddDRcalo {
  static dd4hep::Ref_t create_detector( dd4hep::Detector &description, xml_h xmlElement, dd4hep::SensitiveDetector sensDet ) {    

    // Get the detector description from the xml-tree
    xml_det_t    x_det = xmlElement;
    std::string  name  = x_det.nameStr();
    // Create the detector element
    dd4hep::DetElement drDet( name, x_det.id() );
    // set the sensitive detector type to the DD4hep calorimeter
    dd4hep::xml::Dimension sensDetType = xmlElement.child(_Unicode(sensitive));
    sensDet.setType(sensDetType.typeStr());
    // Get the world volume
    dd4hep::Assembly experimentalHall("hall");
    // Get the dimensions defined in the xml-tree
    xml_comp_t x_barrel ( x_det.child( _Unicode(barrel) ) );
    xml_comp_t x_endcap ( x_det.child( _Unicode(endcap) ) );
    xml_comp_t x_structure ( x_det.child( _Unicode(structure) ) );
    xml_comp_t x_dim ( x_structure.child( _Unicode(dim) ) );
    xml_comp_t x_sipmDim ( x_det.child( _Unicode(sipmDim) ) );

    dd4hep::OpticalSurfaceManager surfMgr = description.surfaceManager();
    dd4hep::OpticalSurface sipmSurfProp = surfMgr.opticalSurface("/world/"+name+"#SiPMSurf");
    dd4hep::OpticalSurface mirrorSurfProp = surfMgr.opticalSurface("/world/"+name+"#MirrorSurf");
    surfMgr.opticalSurface("/world/"+name+"#FilterSurf"); // actual filtering applied in the stepping action

    auto segmentation = dynamic_cast<dd4hep::DDSegmentation::GridDRcalo_k4geo*>( sensDet.readout().segmentation().segmentation() );
    segmentation->setGridSize( x_dim.distance() );
    segmentation->setSipmSize( x_dim.dx() );

    auto paramBarrel = segmentation->paramBarrel();
    paramBarrel->SetInnerX(x_barrel.rmin());
    paramBarrel->SetTowerH(x_barrel.height());
    paramBarrel->SetNumZRot(x_barrel.nphi());
    paramBarrel->SetSipmHeight(x_sipmDim.height());

    auto paramEndcap = segmentation->paramEndcap();
    paramEndcap->SetInnerX(x_endcap.rmin());
    paramEndcap->SetTowerH(x_endcap.height());
    paramEndcap->SetNumZRot(x_endcap.nphi());
    paramEndcap->SetSipmHeight(x_sipmDim.height());

    auto constructor = DRconstructor(x_det);
    constructor.setExpHall(&experimentalHall);
    constructor.setDRparamBarrel(paramBarrel);
    constructor.setDRparamEndcap(paramEndcap);
    constructor.setDescription(&description);
    constructor.setDetElement(&drDet);
    constructor.setSipmSurf(&sipmSurfProp);
    constructor.setSensDet(&sensDet);
    constructor.construct(); // right

    dd4hep::Volume worldVol = description.pickMotherVolume(drDet);
    dd4hep::PlacedVolume hallPlace = worldVol.placeVolume(experimentalHall);
    hallPlace.addPhysVolID("system",x_det.id());
    // connect placed volume and physical volume
    drDet.setPlacement( hallPlace );

    paramBarrel->finalized();
    paramEndcap->finalized();

    return drDet;
  }
} // namespace detector
DECLARE_DETELEMENT(FiberDualReadoutCalo_o1_v01, ddDRcalo::create_detector) // factory method