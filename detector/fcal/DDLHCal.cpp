#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "XML/Utilities.h"
#include <DD4hep/DetFactoryHelper.h>
#include <XML/Layering.h>

#include <string>

using dd4hep::_toString;
using dd4hep::Assembly;
using dd4hep::Box;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::Cone;
using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::DetType;
using dd4hep::IntersectionSolid;
using dd4hep::Layer;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::PolyhedraRegular;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::Rotation3D;
using dd4hep::RotationY;
using dd4hep::RotationZYX;
using dd4hep::Segmentation;
using dd4hep::SensitiveDetector;
using dd4hep::SubtractionSolid;
using dd4hep::Transform3D;
using dd4hep::Translation3D;
using dd4hep::Trapezoid;
using dd4hep::Tube;
using dd4hep::UnionSolid;
using dd4hep::Volume;

using dd4hep::rec::LayeredCalorimeterData;

static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens) {

  std::cout << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "Here is my LHCal" << std::endl;
  std::cout << " and this is the sensitive detector: " << &sens << std::endl;
  sens.setType("calorimeter");
  // Materials
  Material air = theDetector.air();

  // Access to the XML File
  //  dd4hep::xml::DetElement xmlLHCal = xmlHandle;
  //  const std::string detName = xmlLHCal.nameStr();

  xml_det_t xmlLHCal = element;
  std::string detName = xmlLHCal.nameStr();

  DetElement sdet(detName, xmlLHCal.id());

  // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, element, sdet);

  sdet.setTypeFlag(DetType::CALORIMETER | DetType::ENDCAP | DetType::HADRONIC | DetType::FORWARD);

  if (theDetector.buildType() == BUILD_ENVELOPE)
    return sdet;

  //-----------------------------------------------------------------------------------

  dd4hep::xml::Dimension dimensions = xmlLHCal.dimensions();

  // LHCal Dimensions
  const double lhcalInnerR = dimensions.inner_r();
  const double lhcalwidth = dimensions.width();
  const double lhcalheight = dimensions.height();
  const double lhcalInnerZ = dimensions.inner_z();
  const double lhcaloffset = dimensions.offset();
  const double lhcalThickness = Layering(xmlLHCal).totalThickness();
  const double lhcalCentreZ = lhcalInnerZ + lhcalThickness * 0.5;

  double LHcal_cell_size = theDetector.constant<double>("LHcal_cell_size");
  double LHCal_outer_radius = theDetector.constant<double>("LHCal_outer_radius");
  //========== fill data for reconstruction ============================
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData;
  caloData->layoutType = LayeredCalorimeterData::EndcapLayout;
  caloData->inner_symmetry = 0; // hard code cernter pipe hole
  caloData->outer_symmetry = 4; // outer box
  caloData->phi0 = 0;

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = lhcalInnerR;
  caloData->extent[1] = LHCal_outer_radius;
  caloData->extent[2] = lhcalInnerZ;
  caloData->extent[3] = lhcalInnerZ + lhcalThickness;

  // counter for the current layer to be placed
  int thisLayerId = 1;

  // Parameters we have to know about
  dd4hep::xml::Component xmlParameter = xmlLHCal.child(_Unicode(parameter));
  const double mradFullCrossingAngle = xmlParameter.attr<double>(_Unicode(crossingangle));
  std::cout << " The crossing angle is: " << mradFullCrossingAngle << " radian" << std::endl;

  // this needs the full crossing angle, because the beamCal is centred on the
  // outgoing beampipe, and the incoming beampipe is a full crossing agle away
  const Position incomingBeamPipePosition(-1.0 * lhcaloffset, 0.0, 0.0);
  // The extra parenthesis are paramount!
  const Rotation3D incomingBeamPipeRotation((RotationY(0)));
  const Transform3D incomingBPTransform(incomingBeamPipeRotation, incomingBeamPipePosition);

  // Envelope to place the layers in
  Box envelopeBox(lhcalwidth * 0.5, lhcalheight * 0.5, lhcalThickness * 0.5);
  Tube incomingBeamPipe(0.0, lhcalInnerR, lhcalThickness); // we want this to be longer than the LHCal

  SubtractionSolid LHCalModule(envelopeBox, incomingBeamPipe, incomingBPTransform);
  Volume envelopeVol(detName + "_module", LHCalModule, air);
  envelopeVol.setVisAttributes(theDetector, xmlLHCal.visStr());

  ////////////////////////////////////////////////////////////////////////////////
  // Create all the layers
  ////////////////////////////////////////////////////////////////////////////////

  // Loop over all the layer (repeat=NN) sections
  // This is the starting point to place all layers, we need this when we have more than one layer block
  double referencePosition = -lhcalThickness * 0.5;
  for (dd4hep::xml::Collection_t coll(xmlLHCal, _U(layer)); coll; ++coll) {
    dd4hep::xml::Component xmlLayer(coll); // we know this thing is a layer

    // This just calculates the total size of a single layer
    // Why no convenience function for this?
    double layerThickness = 0;

    for (dd4hep::xml::Collection_t l(xmlLayer, _U(slice)); l; ++l)
      layerThickness += xml_comp_t(l).thickness();

    std::cout << "Total Length " << lhcalThickness / dd4hep::cm << " cm" << std::endl;
    std::cout << "Layer Thickness " << layerThickness / dd4hep::cm << " cm" << std::endl;

    // Loop for repeat=NN
    for (int i = 0, repeat = xmlLayer.repeat(); i < repeat; ++i) {

      std::string layer_name = detName + dd4hep::xml::_toString(thisLayerId, "_layer%d");
      Box layer_base(lhcalwidth * 0.5, lhcalheight * 0.5, layerThickness * 0.5);
      SubtractionSolid layer_subtracted;
      layer_subtracted = SubtractionSolid(layer_base, incomingBeamPipe, incomingBPTransform);

      Volume layer_vol(layer_name, layer_subtracted, air);

      int sliceID = 1;
      double inThisLayerPosition = -layerThickness * 0.5;

      double radiator_thickness = 0.0;

      LayeredCalorimeterData::Layer caloLayer;
      caloLayer.cellSize0 = LHcal_cell_size;
      caloLayer.cellSize1 = LHcal_cell_size;

      double nRadiationLengths = 0.;
      double nInteractionLengths = 0.;
      double thickness_sum = 0;

      for (dd4hep::xml::Collection_t collSlice(xmlLayer, _U(slice)); collSlice; ++collSlice) {
        dd4hep::xml::Component compSlice = collSlice;
        const double sliceThickness = compSlice.thickness();
        const std::string sliceName = layer_name + dd4hep::xml::_toString(sliceID, "slice%d");
        Material slice_mat = theDetector.material(compSlice.materialStr());

        Box sliceBase(lhcalwidth * 0.5, lhcalheight * 0.5, sliceThickness * 0.5);
        SubtractionSolid slice_subtracted;

        slice_subtracted = SubtractionSolid(sliceBase, incomingBeamPipe, incomingBPTransform);

        Volume slice_vol(sliceName, slice_subtracted, slice_mat);

        if (compSlice.materialStr().compare("TungstenDens24") == 0)
          radiator_thickness = sliceThickness;

        nRadiationLengths += sliceThickness / (2. * slice_mat.radLength());
        nInteractionLengths += sliceThickness / (2. * slice_mat.intLength());
        thickness_sum += sliceThickness / 2.;

        if (compSlice.isSensitive()) {
          slice_vol.setSensitiveDetector(sens);

#if DD4HEP_VERSION_GE(0, 15)
          // Store "inner" quantities
          caloLayer.inner_nRadiationLengths = nRadiationLengths;
          caloLayer.inner_nInteractionLengths = nInteractionLengths;
          caloLayer.inner_thickness = thickness_sum;
          // Store scintillator thickness
          caloLayer.sensitive_thickness = sliceThickness;
#endif
          // Reset counters to measure "outside" quantitites
          nRadiationLengths = 0.;
          nInteractionLengths = 0.;
          thickness_sum = 0.;
        }

        nRadiationLengths += sliceThickness / (2. * slice_mat.radLength());
        nInteractionLengths += sliceThickness / (2. * slice_mat.intLength());
        thickness_sum += sliceThickness / 2.;

        slice_vol.setAttributes(theDetector, compSlice.regionStr(), compSlice.limitsStr(), compSlice.visStr());
        PlacedVolume pv = layer_vol.placeVolume(slice_vol, Position(0, 0, inThisLayerPosition + sliceThickness * 0.5));

        pv.addPhysVolID("slice", sliceID);
        inThisLayerPosition += sliceThickness;

        ++sliceID;
      } // For all slices in this layer

#if DD4HEP_VERSION_GE(0, 15)
      // Store "outer" quantities
      caloLayer.outer_nRadiationLengths = nRadiationLengths;
      caloLayer.outer_nInteractionLengths = nInteractionLengths;
      caloLayer.outer_thickness = thickness_sum;
#endif
      //-----------------------------------------------------------------------------------------
      caloLayer.distance = lhcalCentreZ + referencePosition; //+0.5*layerThickness ;
      caloLayer.absorberThickness = radiator_thickness;

      caloData->layers.push_back(caloLayer);
      //-----------------------------------------------------------------------------------------

      // Why are we doing this for each layer, this just needs to be done once and then placed multiple times
      // Do we need unique IDs for each piece?
      layer_vol.setVisAttributes(theDetector, xmlLayer.visStr());

      Position layer_pos(0, 0, referencePosition + 0.5 * layerThickness);
      referencePosition += layerThickness;
      PlacedVolume pv = envelopeVol.placeVolume(layer_vol, layer_pos);
      pv.addPhysVolID("layer", thisLayerId);

      ++thisLayerId;

    } // for all layers

  } // for all layer collections

  const Position bcForwardPos(lhcaloffset, 0.0, lhcalCentreZ);
  const Position bcBackwardPos(-1.0 * lhcaloffset, 0.0, -lhcalCentreZ);
  const Rotation3D bcForwardRot(RotationY(0));
  const Rotation3D bcBackwardRot(RotationY(M_PI));

  PlacedVolume pv = envelope.placeVolume(envelopeVol, Transform3D(bcForwardRot, bcForwardPos));
  pv.addPhysVolID("barrel", 1);

  PlacedVolume pv2 = envelope.placeVolume(envelopeVol, Transform3D(bcBackwardRot, bcBackwardPos));
  pv2.addPhysVolID("barrel", 2);

  sdet.addExtension<LayeredCalorimeterData>(caloData);

  return sdet;
}

DECLARE_DETELEMENT(LHCal, create_detector)
