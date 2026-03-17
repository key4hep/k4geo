#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"

#include <cassert>
#include <cmath>
#include <list>

#include "Math/AxisAngle.h"
#include "Math/Vector3D.h"

/*
 * Define a simple assert function which prints a message when the assert will fail
 * so that the user gets some printout
 */
void strawAssert(bool condition, std::string message) {
  if (!condition) {
    std::cout << message << std::endl;
    assert(condition);
  }
}

/*
 * Constructor for straw tube tracker detector
 *
 */
static dd4hep::Ref_t create_straw_tracker(dd4hep::Detector& theDetector, xml_h entities,
                                          dd4hep::SensitiveDetector sens) {

  // XML Detector Element (confusingly also XML::DetElement)
  xml_det_t x_det = entities;
  std::string detName = x_det.nameStr();
  dd4hep::DetElement sdet = dd4hep::DetElement(detName, x_det.id());

  sens.setType("tracker"); // use default "tracker" built-in

  // read detector-level attributes
  xml_dim_t dim = x_det.dimensions();
  double rmin = dim.rmin();
  double rmax = dim.rmax();
  double zmax = dim.zmax();
  double tube_gap = dim.gap();

  // Create top-level envelope for the straw tube tracker
  dd4hep::Tube envelope = dd4hep::Tube(rmin, rmax, zmax);
  dd4hep::Material gas = theDetector.material("Air");
  dd4hep::Volume envelopeVol = dd4hep::Volume(detName + "_envelope", envelope, gas);
  dd4hep::PlacedVolume physvol = theDetector.pickMotherVolume(sdet).placeVolume(envelopeVol);

  // add system ID and identify as barrel (as opposed to endcap +/-1)
  physvol.addPhysVolID("system", sdet.id()).addPhysVolID(_U(side), 0);
  sdet.setPlacement(physvol);

  // Initialize running variables which are updated per superlayer (SL*)
  double SLInnerRadius = rmin;
  int SLNum = 0;
  int tubeNum = 0;
  double tubeThickness = 0.0;
  double mloffset = 0.0;

  // convenience class to add all slice thicknesses together
  dd4hep::Layering layering(x_det);

  // begin the loop over all .xml elements called "layer"
  //
  // _U(layer) is a *macro* not a function.  Gives unicode from arg.
  //
  for (xml_coll_t c(x_det, _U(layer)); c; ++c, ++SLNum) {

    // read superlayer-level attributes
    //
    // Most attr have a top level function but no default value, so we check first if it is provided
    //
    xml_comp_t x_layer = c;
    tubeThickness = 2 * (layering.singleLayerThickness(x_layer));
    double layerRadius =
        SLInnerRadius + 0.5 * (tubeThickness + tube_gap); // first layer is half a tube outside of SL edge
    double delta_phi = asin((tube_gap + tubeThickness) * 0.5 / layerRadius);

    double SLThickness = x_layer.hasAttr(_U(thickness)) ? x_layer.thickness() : -1;
    double SLSectors = x_layer.hasAttr(_U(nsegments)) ? x_layer.nsegments() : 8;
    double SLLayers = x_layer.hasAttr(_U(count)) ? x_layer.count() : 10;
    double SLphiGap = x_layer.hasAttr(_U(gap)) ? x_layer.gap() : 1.5 * dd4hep::cm;
    double SLphiRepeat = x_layer.hasAttr(_U(repeat)) ? x_layer.repeat() : -1;

    double SLoffset = x_layer.hasAttr(_U(offset)) ? x_layer.offset() : 0.0;
    double Angle = x_layer.hasAttr(_U(angle)) ? x_layer.angle() : 0.0;

    // check that the thickness is set
    // and it is a sensible value
    strawAssert(SLThickness > 0 || SLLayers > 0, "ERROR: Each <layer> in straw tube tracker must have either\n"
                                                 "       thickness or count attribute defined.");
    // N.B. 0.866 = sqrt(3)/2 -> minimal radial distance between  staggered tubes centers 
    // in consecutive layers, in units of the diameter (i.e. tube thickness)
    double minThickness = (tubeThickness + tube_gap) * (1 + 0.866 * (SLLayers - 1));
    if (SLThickness < 0) {
      SLThickness = minThickness + SLphiGap; // add default gap to minimum possible thickness
    }

    strawAssert(SLThickness >= minThickness, "ERROR: The specified thickness for a superlayer is less than the\n"
                                             "       minimum allowable value, likely leading to overlaps!");

    // check that the phi repeat is set
    // and it is a sensible value
    double maxRepeat = std::floor((2 * 3.14159 * layerRadius / SLSectors - SLphiGap) / (tubeThickness + tube_gap));
    if (SLphiRepeat < 0) {
      SLphiRepeat = maxRepeat;
    }

    std::cout << "MAX REPEAT=" << maxRepeat << std::endl;

    strawAssert(SLphiRepeat <= maxRepeat, "ERROR: The specified number of tubes in the phi direction is greater"
                                          "       than the maximum allowable value, likely leading to overlaps!");

    // all asserts passed, can start building!

    // make a volume for the multi-layer
    std::string SLName = detName + dd4hep::_toString(SLNum, "_superlayer%d");
    dd4hep::Tube SLTube = dd4hep::Tube(SLInnerRadius, SLInnerRadius + SLThickness, zmax);
    dd4hep::Volume SLVol = dd4hep::Volume(SLName, SLTube, gas);

    // make envelope volume for the tube, which is the same on a per-superlayer basis
    std::string genericTubeName = SLName + "_genericTube";
    dd4hep::Tube singleTube = dd4hep::Tube(0, tubeThickness / 2, zmax);
    dd4hep::Volume singleVol = dd4hep::Volume(genericTubeName, singleTube, gas);

    // populate the envelope volume for the tube with slices: the annular rings of material
    double tubeInnerRadius = 0.0;
    int sliceNum = 0;

    for (xml_coll_t slice(x_layer, _U(slice)); slice; ++slice, ++sliceNum) {
      xml_comp_t x_slice = slice;
      double sliceThickness = x_slice.thickness();
      dd4hep::Material sliceMat = theDetector.material(x_slice.materialStr());
      std::string sliceName = SLName + x_slice.materialStr() + dd4hep::_toString(sliceNum, "slice%d");

      dd4hep::Tube sliceTube = dd4hep::Tube(tubeInnerRadius, tubeInnerRadius + sliceThickness, zmax);
      dd4hep::Volume sliceVol = dd4hep::Volume(sliceName, sliceTube, sliceMat);

      if (x_slice.isSensitive())
        sliceVol.setSensitiveDetector(sens);

      // place the slice in the layer
      singleVol.placeVolume(sliceVol);
      tubeInnerRadius += sliceThickness;

    } // slices (different materials within one tube)

    // place the superlayer in the world volume
    dd4hep::PlacedVolume layerVolPlaced = envelopeVol.placeVolume(SLVol);
    layerVolPlaced.addPhysVolID("superlayer", SLNum);

    // loop over layers, sectors, and tubes all within one superlayer
    // each tube gets a unique placement of the shared abstract "volume"
    // increment tube number on all loops so that it is unique for each placement
    for (int j = 0; j < SLLayers; ++j) {
      for (int l = 0; l < SLSectors; ++l) {
        for (int i = 0; i < SLphiRepeat; ++i) {
          // place the envelope volume in the superlayer envelope
          double phi = l * 2 * dd4hep::pi / SLSectors +
                       (j + 2 * i) * delta_phi * pow(-1, SLNum); // direction of diagonal gap changes per SL

          std::string placedTubeName = SLName + dd4hep::_toString(l, "sector%d") + dd4hep::_toString(j, "layer%d") +
                                       dd4hep::_toString(i, "tube%d");
          dd4hep::DetElement tubeElement = dd4hep::DetElement(sdet, placedTubeName, tubeNum++);

          // Position vector of the tube barycenter
          ROOT::Math::XYZVector axis(layerRadius * cos(phi + mloffset), layerRadius * sin(phi + mloffset), 0);

          // Rotation around the tube barycenter position vector --> stereo angle
          ROOT::Math::AxisAngle rot(axis, Angle);

          // Position (=displacement) vector of the tube barycenter, required for Transform3D
          ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> pos(layerRadius * cos(phi + mloffset),
                                                                                layerRadius * sin(phi + mloffset), 0);
          // Rotation + translation of the tube (ROOT::Math::Transform3D)
          dd4hep::Transform3D tubeTransform = dd4hep::Transform3D(rot, pos);
          dd4hep::PlacedVolume tubePlacedVolume = SLVol.placeVolume(singleVol, tubeTransform);

          // add physical volume IDs corresponding to <readout> <ids> in .xml
          // and set placement of sensitive detector
          tubePlacedVolume.addPhysVolID("sector", l).addPhysVolID("layer", j).addPhysVolID("tube", i);
          tubeElement.setPlacement(tubePlacedVolume);

        } // repeat (tubes in the phi direction)
      } // sectors within superlayer
      layerRadius += (tubeThickness + tube_gap) * 0.866; // equal tube gap
    } // layers within superlayer
    SLInnerRadius += SLThickness;
    mloffset += SLoffset;

  } // multi-layers aka _U(layer)
  return sdet;
}
DECLARE_DETELEMENT(StrawTubeTracker_o1_v01, create_straw_tracker)
