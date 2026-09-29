/*
 @author: Mahmoud Althakeel
mahmoud.althakeel@cern.ch

Factory for IDEA muon system — o1_v02 (sequence-based layer configuration)

New in this version vs o1_v01:
  <Barrel> and <Endcap> now accept an explicit ordered sequence of <layer>
  child elements instead of fixed numDetectorLayers/numYokes counts.
  This allows arbitrary detector/radiator arrangements, e.g.:
    two detector layers → one radiator → one detector → one radiator

Expected XML structure:
<detector type="muonSystem_o1_v02" ...>
  <dimensions x="chamberHalfThickness" y="chamberHalfY" z="chamberHalfZ" .../>
  <sensitive type="tracker"/>
  <generalParameters numSides="..." overlapY="..." overlapZ="..." clearance="..."/>

  <Barrel rmin="..." length="...">
    <layer type="detector"/>
    <layer type="detector"/>
    <layer type="radiator" thickness="300*mm" material="G4_Fe"/>
    <layer type="detector"/>
    <layer type="radiator" thickness="300*mm" material="G4_Fe"/>
  </Barrel>

  <Endcap rmin="..." z_offset="...">
    <layer type="detector"/>
    <layer type="detector"/>
    <layer type="radiator" thickness="300*mm" material="G4_Fe"/>
    <layer type="detector"/>
    <layer type="radiator" thickness="300*mm" material="G4_Fe"/>
  </Endcap>

  <!-- mRWELL chamber slices (same as o1) -->
  <slice x="..." material="..." [sensitive="true"] [vis="..."]/>
  ...
</detector>

Sequence rules:
  - Layers are placed in the order listed (innermost = first).
  - Each <layer type="detector"> has the same chamber thickness derived from
    the <dimensions> element (same as o1).
  - Each <layer type="radiator"> requires thickness and material attributes.
  - Endcap sequence is mirrored (reversed) for the negative endcap so that
    detector layer 0 is always closest to the IP in both endcaps.
*/

#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/DocumentHandler.h"
#include "XML/Utilities.h"
#include "XML/XMLElements.h"
#include <cmath>
#include <sstream>
#include <vector>
#include <string>

using namespace std;
using namespace dd4hep;
using dd4hep::rec::SurfaceType;
using dd4hep::rec::Vector3D;
using dd4hep::rec::VolCylinder;
using dd4hep::rec::VolPlane;
using dd4hep::rec::volSurfaceList;

// ---------------------------------------------------------------------------
// Helper: one entry in the ordered layer sequence
// ---------------------------------------------------------------------------
struct LayerSpec {
  std::string      type;      // "detector" or "radiator(yoke)"
  double           thickness; // full thickness (not half-length)
  dd4hep::Material matObj;    // set only for type=="radiator"
};

// ---------------------------------------------------------------------------
static dd4hep::Ref_t createmuonSystem_o1_v02(dd4hep::Detector&       lcdd,
                                                     dd4hep::xml::Handle_t   xmlElement,
                                                     dd4hep::SensitiveDetector sensDet) {

  dd4hep::xml::DetElement xmlDet = static_cast<dd4hep::xml::DetElement>(xmlElement);
  std::string  name        = xmlDet.nameStr();
  dd4hep::DetElement detElement(name, xmlDet.id());
  dd4hep::Material   mat  = lcdd.material("Air");
  dd4hep::Volume     experimentalHall = lcdd.pickMotherVolume(detElement);

  xml_comp_t dimensions(xmlDet.dimensions());

  // -------------------------------------------------------------------------
  //  General parameters
  // -------------------------------------------------------------------------
  auto generalParameters = xmlElement.child(_Unicode(generalParameters));
  int    numSides   = generalParameters.attr<int>("numSides");
  double overlapY   = generalParameters.attr<double>("overlapY");
  double overlapZ   = generalParameters.attr<double>("overlapZ");
  double clearance  = generalParameters.attr<double>("clearance");

  double shapeAngle         = 360.0 / numSides;
  double shapeAngle_radians = M_PI / numSides;
  double angle_clearance    = 0.0;

  // -------------------------------------------------------------------------
  //  Parse Barrel layer sequence 
  // -------------------------------------------------------------------------
  auto BarrelXML = xmlElement.child(_Unicode(Barrel));
  double radius      = BarrelXML.attr<double>("rmin");
  double barrelLength= BarrelXML.attr<double>("length");

  //count detector layers so we can compute detectorVolumeThickness
  int numBarrelDetectorLayers = 0;
  for (xml_coll_t layerIt(BarrelXML, _Unicode(layer)); layerIt; ++layerIt) {
      xml_comp_t layer(layerIt);
      if (layer.attr<std::string>("type") == "detector") {
          ++numBarrelDetectorLayers;
      }
  }

  // -------------------------------------------------------------------------
  //  Chamber / detector-layer geometry 
  // -------------------------------------------------------------------------
  double sideLengthEstimate = 2.0 * radius * std::tan(shapeAngle_radians);
  double chamberAngle_rad   = std::atan(dimensions.x() / dimensions.y());
  double rectangleThickness = (2.0 * dimensions.y()) * std::sin(chamberAngle_rad)
                            + (2.0 * dimensions.x()) * std::cos(chamberAngle_rad)
                            + 1.2 * clearance;
  double rectangleAngle_rad = std::atan(rectangleThickness / dimensions.z());
  double detectorVolumeThickness, endcapDetectorEnvZ;
  if (sideLengthEstimate <= (2.0 * dimensions.y()) && numBarrelDetectorLayers == 1) {
    detectorVolumeThickness = rectangleThickness * 1.33;
    endcapDetectorEnvZ      = 2.0 * detectorVolumeThickness;
  } else {
    detectorVolumeThickness = (2.0 * dimensions.z()) * std::sin(rectangleAngle_rad)
                            + rectangleThickness       * std::cos(rectangleAngle_rad);
    endcapDetectorEnvZ      = detectorVolumeThickness;
  }

  // -------------------------------------------------------------------------
  //  Parse Barrel layer sequence 
  // -------------------------------------------------------------------------
  std::vector<LayerSpec> barrelSeq;
  for (xml_coll_t il(BarrelXML, _Unicode(layer)); il; ++il) {
    xml_comp_t le(il);
    LayerSpec s;
    s.type = le.attr<std::string>("type");
    if (s.type == "radiator") {
      s.thickness = le.attr<double>("thickness");
      s.matObj    = lcdd.material(le.attr<std::string>("material"));
    } else {
      s.thickness = 2.0 * detectorVolumeThickness;
    }
    barrelSeq.push_back(s);
  }

  // Pre-compute per-entry detector index (used for physVolID "layer")
  std::vector<int> barrelDetIdx(barrelSeq.size(), -1);
  {
    int dIdx = 0;
    for (int i = 0; i < (int)barrelSeq.size(); ++i)
      if (barrelSeq[i].type == "detector") barrelDetIdx[i] = dIdx++;
  }

  // Cumulative radial positions and locate last detector layer
  double totalBarrelThickness = 0.0;
  double lastDetLayerRMin     = radius;
  double lastDetLayerRMax     = radius;
  {
    double cur = radius;
    for (int i = 0; i < (int)barrelSeq.size(); ++i) {
      if (barrelSeq[i].type == "detector") {
        lastDetLayerRMin = cur;
        lastDetLayerRMax = cur + barrelSeq[i].thickness;
      }
      cur += barrelSeq[i].thickness;
    }
    totalBarrelThickness = cur - radius;
  }
  double barrelRMax = radius + totalBarrelThickness;

  // -------------------------------------------------------------------------
  //  Parse Endcap layer sequence
  // -------------------------------------------------------------------------
  auto EndcapXML = xmlElement.child(_Unicode(Endcap));
  double endcapDetectorLayerInnerRadius = EndcapXML.attr<double>("rmin");
  double endcapZOffset                  = EndcapXML.attr<double>("z_offset");

  std::vector<LayerSpec> endcapSeq;
  int numEndcapDetectorLayers = 0;
  for (xml_coll_t il(EndcapXML, _Unicode(layer)); il; ++il) {
    xml_comp_t le(il);
    LayerSpec s;
    s.type = le.attr<std::string>("type");
    if (s.type == "radiator") {
      s.thickness = le.attr<double>("thickness");
      s.matObj    = lcdd.material(le.attr<std::string>("material"));
    } else {
      s.thickness = 2.0 * endcapDetectorEnvZ;
      ++numEndcapDetectorLayers;
    }
    endcapSeq.push_back(s);
  }

  // Pre-compute per-entry detector index for endcap
  std::vector<int> endcapDetIdx(endcapSeq.size(), -1);
  {
    int dIdx = 0;
    for (int i = 0; i < (int)endcapSeq.size(); ++i)
      if (endcapSeq[i].type == "detector") endcapDetIdx[i] = dIdx++;
  }

  double EndcaptotalLength = 0.0;
  for (const auto& s : endcapSeq) EndcaptotalLength += s.thickness;

  double barreltotalLength = barrelLength + 2.0 * EndcaptotalLength;
  double endcapOffset      = endcapZOffset + EndcaptotalLength / 2.0;

  // Endcap outer radius = inner edge of last barrel detector layer (enclosed by it)
  double endcapDetectorLayerOuterRadius  = lastDetLayerRMin;
  double endcapRadiatorLayerInnerRadius  = endcapDetectorLayerInnerRadius;
  double endcapRadiatorLayerOuterRadius  = endcapDetectorLayerOuterRadius;

  // Endcap shape helpers (same as o1)
  double endcapDetectorSideLength =
      (2.0 * (endcapDetectorLayerInnerRadius + 2.0 * dimensions.y()) * std::tan(shapeAngle_radians)) + 2.0 * clearance;
  double endcapDetectorSideTrapLength =
      (2.0 * endcapDetectorLayerOuterRadius * std::tan(shapeAngle_radians)) + 2.0 * clearance;
  double endcapDetectorSideTrapYLength =
      endcapDetectorLayerOuterRadius - 2.0 * dimensions.z() - endcapDetectorLayerInnerRadius;
  double endcapDetectorSideBoxLength  = 2.0 * endcapDetectorLayerOuterRadius * std::tan(shapeAngle_radians);
  double endcapDetectorSideBoxYLength = 2.0 * dimensions.z();

  double endcapRadiatorSideLength  = 2.0 * endcapRadiatorLayerInnerRadius * std::tan(shapeAngle_radians);
  double endcapRadiatorSideLength2 = 2.0 * endcapRadiatorLayerOuterRadius * std::tan(shapeAngle_radians);

  double endcapDetectorYLength = endcapDetectorLayerOuterRadius - endcapDetectorLayerInnerRadius;
  double endcapYLength         = endcapRadiatorLayerOuterRadius - endcapRadiatorLayerInnerRadius;
  double endcapRemainderZ =
      std::fmod((endcapDetectorYLength - 2.0 * clearance), (2.0 * dimensions.z() - overlapZ)) / (2.0 * dimensions.z())
      - (2.0 * clearance / dimensions.z());

  int barrelIdCounter = 1;
  int endcapIdCounter = 1;

  // -------------------------------------------------------------------------
  //  System envelope
  // -------------------------------------------------------------------------
  dd4hep::PolyhedraRegular BarrelEnv(numSides, radius, barrelRMax, barreltotalLength);
  dd4hep::PolyhedraRegular EndcapEnv(numSides, endcapDetectorLayerInnerRadius,
                                      endcapDetectorLayerOuterRadius, EndcaptotalLength);

  dd4hep::Position    unionPos (0.0, 0.0,  endcapOffset);
  dd4hep::Position    unionPos2(0.0, 0.0, -endcapOffset);
  dd4hep::Rotation3D  unionRot(dd4hep::RotationY(0.0 * dd4hep::degree));
  dd4hep::Transform3D unionTransform (unionRot, unionPos);
  dd4hep::Transform3D unionTransform2(unionRot, unionPos2);

  dd4hep::UnionSolid barrelAndPositiveEndcap(BarrelEnv, EndcapEnv, unionTransform);
  dd4hep::UnionSolid systemEnvelope(barrelAndPositiveEndcap, EndcapEnv, unionTransform2);
  dd4hep::Volume     detectorVolume(name, systemEnvelope, mat);

  dd4hep::Position      detectorTrans(0., 0., 0.);
  dd4hep::PlacedVolume  detectorPhys = experimentalHall.placeVolume(
      detectorVolume, dd4hep::Transform3D(dd4hep::RotationZ(shapeAngle_radians), detectorTrans));
  detectorPhys.addPhysVolID("system", xmlDet.id());
  detElement.setPlacement(detectorPhys);
  detectorVolume.setVisAttributes(lcdd.visAttributes("no_vis"));

  // -------------------------------------------------------------------------
  //  BARREL volume
  // -------------------------------------------------------------------------
  dd4hep::PolyhedraRegular BarrelEnvBase(numSides, radius, barrelRMax, barrelLength);
  dd4hep::PolyhedraRegular BarrelLastDetExt(numSides, lastDetLayerRMin, lastDetLayerRMax, barreltotalLength);

  std::string barrelName = name + "-Barrel" + std::to_string(0);
  dd4hep::Position    barrelUnionPos(0., 0., 0.);
  dd4hep::Transform3D barrelUnionTransform(dd4hep::Rotation3D(dd4hep::RotationY(0.)), barrelUnionPos);
  dd4hep::UnionSolid  barrelUnion(BarrelEnvBase, BarrelLastDetExt, barrelUnionTransform);
  dd4hep::Volume      BarrelVolume(barrelName, barrelUnion, mat);

  dd4hep::PlacedVolume barrelPhys = detectorVolume.placeVolume(
      BarrelVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), dd4hep::Position(0., 0., 0.)));
  barrelPhys.addPhysVolID("type", 0);
  dd4hep::DetElement BarrelDE(detElement, name + "-BarrelDE", 0);
  BarrelDE.setPlacement(barrelPhys);
  BarrelVolume.setVisAttributes(lcdd.visAttributes("no_vis"));

  // -------------------------------------------------------------------------
  //  BARREL — unified sequence loop (detector layers + radiators)
  // -------------------------------------------------------------------------
  double currentRMin = radius;
  for (int seqIdx = 0; seqIdx < (int)barrelSeq.size(); ++seqIdx) {
    const auto& spec     = barrelSeq[seqIdx];
    double      layerRMin = currentRMin;
    double      layerRMax = layerRMin + spec.thickness;
    currentRMin = layerRMax;

    // ---- DETECTOR LAYER ------------------------------------------------
    if (spec.type == "detector") {
      int    numBarrelLayer    = barrelDetIdx[seqIdx];
      bool   isLastDetLayer    = (std::abs(layerRMax - lastDetLayerRMax) < 1e-6);
      double barrelLayerLength = isLastDetLayer ? barreltotalLength : barrelLength;
      double barrelLayerRMid   = (layerRMin + layerRMax) / 2.0;

      dd4hep::PolyhedraRegular BarrelDetectorLayer(numSides, layerRMin, layerRMax, barrelLayerLength);
      std::string barrelDetectorName = name + "-BarrelDetectorLayer" + std::to_string(numBarrelLayer + 1);
      dd4hep::Volume BarrelDetectorLayerVolume(barrelDetectorName, BarrelDetectorLayer, mat);

      double sideLength  = 2.0 * (barrelLayerRMid - detectorVolumeThickness / 2.0) * std::tan(shapeAngle_radians);
      double sideLength2 = 2.0 * (barrelLayerRMid + detectorVolumeThickness / 2.0) * std::tan(shapeAngle_radians);
      double sideEnvX    = detectorVolumeThickness / 2.0;
      double sideEnvY    = sideLength  / 2.0 + clearance;
      double sideEnvY2   = sideLength2 / 2.0 + clearance;
      double sideEnvZ    = barrelLayerLength / 2.0;
      double remainderZ  =
          std::fmod((barrelLayerLength - 2.0 * clearance), (2.0 * dimensions.z() - overlapZ)) / (2.0 * dimensions.z())
          - (2.0 * clearance / dimensions.z());

      dd4hep::PlacedVolume detectorLayerPhys = BarrelVolume.placeVolume(
          BarrelDetectorLayerVolume,
          dd4hep::Transform3D(dd4hep::RotationZ(0.), dd4hep::Position(0., 0., 0.)));
      detectorLayerPhys.addPhysVolID("layer", numBarrelLayer);
      dd4hep::DetElement BarrelDetectorLayerDE(
          BarrelDE, name + "-Barrel_DetectorLayerDE_" + std::to_string(numBarrelLayer + 1), numBarrelLayer + 1);
      BarrelDetectorLayerDE.setPlacement(detectorLayerPhys);
      BarrelDetectorLayerVolume.setVisAttributes(lcdd.visAttributes("no_vis"));

      // ---- Sides --------------------------------------------------------
      for (int side = 0; side < numSides; ++side) {
        int sideID = (seqIdx + 1) * 100 + (side + 1);

        dd4hep::Box sideEnvelope (sideEnvX, sideEnvY,  sideEnvZ);
        dd4hep::Box sideEnvelope2(sideEnvX, sideEnvY2, sideEnvZ);
        std::string sideName = dd4hep::xml::_toString(sideID, "side%d");
        dd4hep::Volume sideVol (sideName, sideEnvelope,  mat);
        dd4hep::Volume sideVol2(sideName, sideEnvelope2, mat);

        double angle_degrees = shapeAngle * side;
        double angle_radians = angle_degrees * M_PI / 180.0;
        double sideXOffset  = (barrelLayerRMid - detectorVolumeThickness / 2.0) * std::cos(angle_radians + shapeAngle_radians);
        double sideYOffset  = (barrelLayerRMid - detectorVolumeThickness / 2.0) * std::sin(angle_radians + shapeAngle_radians);
        double sideXOffset2 = (barrelLayerRMid + detectorVolumeThickness / 2.0) * std::cos(angle_radians + shapeAngle_radians);
        double sideYOffset2 = (barrelLayerRMid + detectorVolumeThickness / 2.0) * std::sin(angle_radians + shapeAngle_radians);

        dd4hep::RotationZ  sideRotationZ(angle_radians + shapeAngle_radians + angle_clearance);
        dd4hep::Rotation3D sideRotation = dd4hep::Rotation3D(sideRotationZ);

        dd4hep::Position     sideTrans;
        dd4hep::PlacedVolume sidePhys;
        if (side % 2 == 0) {
          sideTrans = dd4hep::Position(sideXOffset,  sideYOffset,  0.0);
          sidePhys  = BarrelDetectorLayerVolume.placeVolume(sideVol,  dd4hep::Transform3D(sideRotation, sideTrans));
        } else {
          sideTrans = dd4hep::Position(sideXOffset2, sideYOffset2, 0.0);
          sidePhys  = BarrelDetectorLayerVolume.placeVolume(sideVol2, dd4hep::Transform3D(sideRotation, sideTrans));
        }
        dd4hep::DetElement sideDE(BarrelDetectorLayerDE, sideName + "DE", sideID);
        sideDE.setPlacement(sidePhys);
        sideVol.setVisAttributes(lcdd,  xmlDet.visStr());
        sideVol2.setVisAttributes(lcdd, xmlDet.visStr());

        // ---- Rectangles -------------------------------------------------
        int numRectangles = (int)(barrelLayerLength / (2.0 * dimensions.z() - overlapZ));
        for (int rectangle = 0; rectangle <= numRectangles; ++rectangle) {
          double rectangleEnvY = (side % 2 == 0) ? sideEnvY : sideEnvY2;
          double rectangleEnvZ = (rectangle == numRectangles)
                                     ? (remainderZ * dimensions.z() + clearance / 4.0)
                                     : (dimensions.z() + clearance / 2.0);
          double rectangleEnvX = (rectangleEnvY <= dimensions.y())
                                     ? dimensions.x()
                                     : detectorVolumeThickness / 4.5;

          dd4hep::Box rectangleEnvelope(rectangleEnvX, rectangleEnvY, rectangleEnvZ);
          std::string rectangleEnvelopeName = (rectangle == numRectangles)
              ? dd4hep::xml::_toString(rectangle, "rectangleRemainderEnvelope%d")
              : dd4hep::xml::_toString(rectangle, "rectangleEnvelope%d");
          dd4hep::Volume rectangleEnvVol(rectangleEnvelopeName, rectangleEnvelope, mat);

          double rectangleEnvZPos, rectangleRemainderEnvXPos;
          if (rectangle == numRectangles) {
            rectangleRemainderEnvXPos = sideEnvX / 2.0;
            rectangleEnvZPos = barrelLayerLength / 2.0 - (remainderZ * dimensions.z())
                               - (rectangle * (2.0 * dimensions.z() - overlapZ)) - clearance;
          } else {
            rectangleRemainderEnvXPos = 0.0;
            rectangleEnvZPos = barrelLayerLength / 2.0 - dimensions.z()
                               - (rectangle * (2.0 * dimensions.z() - overlapZ)) - clearance;
          }

          double yRotation = std::atan(rectangleEnvX / (rectangleEnvZ - 2.0 * overlapZ));
          double yRemainderRotation;
          dd4hep::RotationY rotationY(yRotation);
          if (rectangleEnvZ <= (dimensions.z() + clearance / 2.0) / 10.0)
            yRemainderRotation = yRotation;
          else
            yRemainderRotation = std::atan(rectangleEnvX / (rectangleEnvZ - 2.0 * overlapZ));
          dd4hep::RotationY remainderRotationY(yRemainderRotation);

          dd4hep::Position rectangeEnvelopeTrans;
          if (rectangle == numRectangles) {
            rectangeEnvelopeTrans = (rectangleEnvZ <= (dimensions.z() + clearance / 2.0) / 10.0)
                ? dd4hep::Position(rectangleRemainderEnvXPos, 0.0, rectangleEnvZPos)
                : dd4hep::Position(0.0, 0.0, rectangleEnvZPos);
          } else {
            rectangeEnvelopeTrans = dd4hep::Position(0.0, 0.0, rectangleEnvZPos);
          }

          dd4hep::PlacedVolume rectangleEnvelopePhys;
          if (rectangle == numRectangles) {
            if (side % 2 == 0)
              rectangleEnvelopePhys = sideVol.placeVolume(rectangleEnvVol, dd4hep::Transform3D(remainderRotationY, rectangeEnvelopeTrans));
            else
              rectangleEnvelopePhys = sideVol2.placeVolume(rectangleEnvVol, dd4hep::Transform3D(remainderRotationY, rectangeEnvelopeTrans));
          } else {
            if (side % 2 == 0)
              rectangleEnvelopePhys = sideVol.placeVolume(rectangleEnvVol, dd4hep::Transform3D(rotationY, rectangeEnvelopeTrans));
            else
              rectangleEnvelopePhys = sideVol2.placeVolume(rectangleEnvVol, dd4hep::Transform3D(rotationY, rectangeEnvelopeTrans));
          }
          dd4hep::DetElement rectangleEnvelopeDE(sideDE, rectangleEnvelopeName + "DE", rectangle);
          rectangleEnvelopeDE.setPlacement(rectangleEnvelopePhys);
          rectangleEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

          // ---- Chambers -------------------------------------------------
          int numChambersInRectangle = (rectangleEnvY <= dimensions.y())
              ? 0
              : (int)(2.0 * rectangleEnvY / (2.0 * dimensions.y() - overlapY));

          for (int chamberIndex = 0; chamberIndex <= numChambersInRectangle; ++chamberIndex) {
            std::stringstream bns;
            bns << "-MuRWELL_Barrel_" << barrelIdCounter++;
            std::string BarrelChamberName = name + bns.str();

            double rectangleRemainderY, rectangleRemainderREnvYPos;
            if (numChambersInRectangle == 0) {
              rectangleRemainderY = std::abs(std::fmod((2.0 * rectangleEnvY - clearance),
                                                        (2.0 * dimensions.y() - clearance))) / (2.0 * dimensions.y());
              rectangleRemainderREnvYPos = (chamberIndex * 2.0 * dimensions.y()) - (overlapY * chamberIndex)
                                         + dimensions.y_offset() - rectangleEnvY
                                         + rectangleRemainderY * dimensions.y() + clearance / 20.0;
            } else {
              rectangleRemainderY = std::fmod(2.0 * (rectangleEnvY - clearance),
                                               (2.0 * dimensions.y() - overlapY)) / (2.0 * dimensions.y());
              rectangleRemainderREnvYPos = (chamberIndex * 2.0 * dimensions.y()) - (overlapY * chamberIndex)
                                         + dimensions.y_offset() - rectangleEnvY
                                         + rectangleRemainderY * dimensions.y() + 1.5 * clearance;
            }

            double chamberHalfZ = (rectangle == numRectangles) ? remainderZ * dimensions.z() : dimensions.z();
            dd4hep::Box envelope(dimensions.x(), dimensions.y(), chamberHalfZ);
            dd4hep::Volume envVolume(BarrelChamberName, envelope, lcdd.material(dimensions.materialStr()));

            dd4hep::Box rectangleRemainderYEnvelope(dimensions.x(), rectangleRemainderY * dimensions.y(), chamberHalfZ);
            dd4hep::Volume rectangleRemainderYEnvVolume(
                BarrelChamberName + "rectangleRemainderY", rectangleRemainderYEnvelope,
                lcdd.material(dimensions.materialStr()));

            double envYPos = (chamberIndex * 2.0 * dimensions.y()) - (overlapY * chamberIndex)
                           + dimensions.y_offset() - rectangleEnvY + dimensions.y() + clearance / 20.0;

            double zRotation = 0.0, rectangleRemainderZRotation = 0.0;
            if (numChambersInRectangle != 0) {
              zRotation = std::atan(dimensions.x() / (dimensions.y() - 2.0 * overlapY));
              rectangleRemainderZRotation = std::atan(dimensions.x() /
                  (rectangleRemainderY * dimensions.y() - overlapY * 1.5));
            }
            dd4hep::RotationZ chamberRotation(zRotation);
            dd4hep::RotationZ rectangleRemainderRotationZ(rectangleRemainderZRotation);

            auto Slices   = xmlElement.children(_Unicode(slice));
            auto numSlices= xmlElement.numChildren(_Unicode(slice), true);
            dd4hep::xml::Handle_t slice(Slices.reset());
            int sensitiveSliceIndex = 0;

            if (chamberIndex == numChambersInRectangle) {
              // --- Remainder chamber ---
              dd4hep::PlacedVolume rrmPhys = rectangleEnvVol.placeVolume(
                  rectangleRemainderYEnvVolume,
                  dd4hep::Transform3D(rectangleRemainderRotationZ,
                                      dd4hep::Position(dimensions.x_offset(), rectangleRemainderREnvYPos, 0.0)));
              rrmPhys.addPhysVolID("chamber", barrelIdCounter);
              dd4hep::DetElement rrmDE(rectangleEnvelopeDE, BarrelChamberName, barrelIdCounter);
              rrmDE.setPlacement(rrmPhys);
              rectangleRemainderYEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());

              double sliceXOffset = -dimensions.x();
              for (unsigned sliceIdx = 0; sliceIdx < numSlices; ++sliceIdx) {
                dd4hep::xml::DetElement sliceDet = static_cast<dd4hep::xml::DetElement>(slice);
                dd4hep::Box sliceShape(sliceDet.x(), rectangleRemainderY * dimensions.y(), chamberHalfZ);
                std::string sliceName = dd4hep::xml::_toString(sliceIdx, "slice%d");
                dd4hep::Volume sliceVolume(sliceName, sliceShape, lcdd.material(slice.attr<std::string>("material")));
                dd4hep::PlacedVolume slicePV = rectangleRemainderYEnvVolume.placeVolume(
                    sliceVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.),
                                                     dd4hep::Position(sliceXOffset + sliceDet.x(), 0.0, 0.0)));
                if (slice.hasAttr("vis"))       sliceVolume.setVisAttributes(lcdd, sliceDet.visStr());
                if (slice.hasAttr("sensitive") && sliceDet.isSensitive()) {
                  dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                  sensDet.setType(sdType.typeStr());
                  sliceVolume.setSensitiveDetector(sensDet);
                  slicePV.addPhysVolID("slice", sensitiveSliceIndex);
                  dd4hep::DetElement sliceDE(rrmDE, "slice_" + std::to_string(sensitiveSliceIndex), sensitiveSliceIndex);
                  sliceDE.setPlacement(slicePV);
                  ++sensitiveSliceIndex;
                }
                sliceXOffset += 2.0 * sliceDet.x();
                slice.m_node = Slices.next();
              }

            } else {
              // --- Full chamber ---
              dd4hep::PlacedVolume envPhys = rectangleEnvVol.placeVolume(
                  envVolume, dd4hep::Transform3D(chamberRotation, dd4hep::Position(dimensions.x_offset(), envYPos, 0.0)));
              envPhys.addPhysVolID("chamber", barrelIdCounter);
              dd4hep::DetElement envDE(rectangleEnvelopeDE, BarrelChamberName, barrelIdCounter);
              envDE.setPlacement(envPhys);
              envVolume.setVisAttributes(lcdd, xmlDet.visStr());

              double sliceXOffset = -dimensions.x();
              for (unsigned sliceIdx = 0; sliceIdx < numSlices; ++sliceIdx) {
                dd4hep::xml::DetElement sliceDet = static_cast<dd4hep::xml::DetElement>(slice);
                dd4hep::Box sliceShape(sliceDet.x(), dimensions.y(), chamberHalfZ);
                std::string sliceName = dd4hep::xml::_toString(sliceIdx, "slice%d");
                dd4hep::Volume sliceVolume(sliceName, sliceShape, lcdd.material(slice.attr<std::string>("material")));
                dd4hep::PlacedVolume slicePV = envVolume.placeVolume(
                    sliceVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.),
                                                     dd4hep::Position(sliceXOffset + sliceDet.x(), 0.0, 0.0)));
                if (slice.hasAttr("vis"))       sliceVolume.setVisAttributes(lcdd, sliceDet.visStr());
                if (slice.hasAttr("sensitive") && sliceDet.isSensitive()) {
                  dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                  sensDet.setType(sdType.typeStr());
                  sliceVolume.setSensitiveDetector(sensDet);
                  slicePV.addPhysVolID("slice", sensitiveSliceIndex);
                  dd4hep::DetElement sliceDE(envDE, "slice_" + std::to_string(sensitiveSliceIndex), sensitiveSliceIndex);
                  sliceDE.setPlacement(slicePV);
                  ++sensitiveSliceIndex;
                }
                sliceXOffset += 2.0 * sliceDet.x();
                slice.m_node = Slices.next();
              }
            }
          } // chamberIndex
        }   // rectangle
      }     // side

    // ---- RADIATOR LAYER ------------------------------------------------
    } else {
      double layerRMid = (layerRMin + layerRMax) / 2.0;
      double barrelRadiatorSideLength  = 2.0 * (layerRMid - spec.thickness / 2.0) * std::tan(shapeAngle_radians);
      double barrelRadiatorSideLength2 = 2.0 * (layerRMid + spec.thickness / 2.0) * std::tan(shapeAngle_radians);

      dd4hep::PolyhedraRegular BarrelRadiatorLayer(numSides, layerRMin, layerRMax, barrelLength);
      std::string barrelRadiatorEnvName = name + "-BarrelRadiatorLayer" + std::to_string(seqIdx + 1);
      dd4hep::Volume BarrelRadiatorLayerVolume(barrelRadiatorEnvName, BarrelRadiatorLayer, mat);

      for (int side = 0; side < numSides; ++side) {
        int sideID = (seqIdx + 1) * 100 + (side + 1);
        dd4hep::Trapezoid barrelRadiatorEnvelope(
            barrelLength / 2.0, barrelLength / 2.0,
            barrelRadiatorSideLength / 2.0, barrelRadiatorSideLength2 / 2.0,
            spec.thickness / 2.0);
        std::string barrelRadiatorEnvelopeName = name + "-BarrelRadiatorSide" + std::to_string(sideID);
        dd4hep::Volume barrelRadiatorEnvVol(barrelRadiatorEnvelopeName, barrelRadiatorEnvelope, mat);

        double angle_degrees = shapeAngle * side;
        double angle_radians = angle_degrees * M_PI / 180.0;
        double bRXOffset = layerRMid * std::cos(angle_radians + shapeAngle_radians);
        double bRYOffset = layerRMid * std::sin(angle_radians + shapeAngle_radians);

        dd4hep::RotationZ  bRRotZ(angle_radians + shapeAngle_radians);
        dd4hep::RotationY  bRRotY(90.0 * dd4hep::degree);
        dd4hep::Rotation3D bRRot = bRRotZ * bRRotY;

        dd4hep::PlacedVolume barrelRadiatorEnvelopePhys = BarrelRadiatorLayerVolume.placeVolume(
            barrelRadiatorEnvVol,
            dd4hep::Transform3D(bRRot, dd4hep::Position(bRXOffset, bRYOffset, 0.0)));
        dd4hep::DetElement barrelRadiatorEnvelopeDE(detElement, barrelRadiatorEnvelopeName + "DE", sideID);
        barrelRadiatorEnvelopeDE.setPlacement(barrelRadiatorEnvelopePhys);
        barrelRadiatorEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

        dd4hep::Trapezoid barrelRadiator(
            barrelLength / 2.0, barrelLength / 2.0,
            barrelRadiatorSideLength / 2.0, barrelRadiatorSideLength2 / 2.0,
            spec.thickness / 2.0);
        std::string barrelRadiatorName = dd4hep::xml::_toString(side, "barrelRadiator%d");
        dd4hep::Volume barrelRadiatorVol(barrelRadiatorName, barrelRadiator, spec.matObj);
        dd4hep::PlacedVolume barrelRadiatorPhys = barrelRadiatorEnvVol.placeVolume(
            barrelRadiatorVol, dd4hep::Transform3D(dd4hep::RotationY(0.), dd4hep::Position(0, 0, 0)));
        dd4hep::DetElement barrelRadiatorDE(barrelRadiatorEnvelopeDE, barrelRadiatorName + "DE", sideID);
        barrelRadiatorDE.setPlacement(barrelRadiatorPhys);
        barrelRadiatorVol.setVisAttributes(lcdd.visAttributes("yoke_vis"));
      }

      dd4hep::PlacedVolume radiatorLayerPhys = BarrelVolume.placeVolume(
          BarrelRadiatorLayerVolume,
          dd4hep::Transform3D(dd4hep::RotationZ(0.), dd4hep::Position(0., 0., 0.)));
      int radLayerDE_id = seqIdx + 1;
      dd4hep::DetElement BarrelRadiatorLayerDE(BarrelDE, barrelRadiatorEnvName + "DE", radLayerDE_id);
      BarrelRadiatorLayerDE.setPlacement(radiatorLayerPhys);
      BarrelRadiatorLayerVolume.setVisAttributes(lcdd.visAttributes("no_vis"));
    }
  }

  // =========================================================================
  //  ENDCAP
  // =========================================================================
  int numSeq = (int)endcapSeq.size();

  for (int endcapType = -1; endcapType < 2; ++endcapType) {
    if (endcapType == 0) continue;

    std::string          EndcapName  = name + "-Endcap" + std::to_string(endcapType);
    dd4hep::Volume       endcapVolume(EndcapName, EndcapEnv, mat);
    dd4hep::Position     endcapTrans(0., 0., endcapType * endcapOffset);
    dd4hep::PlacedVolume endcapPhys = detectorVolume.placeVolume(
        endcapVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), endcapTrans));
    endcapPhys.addPhysVolID("type", endcapType);
    dd4hep::DetElement EndcapDE(detElement, name + "EndcapDE_" + std::to_string(endcapType), endcapType);
    EndcapDE.setPlacement(endcapPhys);
    endcapVolume.setVisAttributes(lcdd.visAttributes("no_vis"));

    // For negative endcap iterate sequence in reverse so layer-0 stays
    // closest to the IP in both endcaps (mirror symmetry around z=0).
    double currentZ = -EndcaptotalLength / 2.0;
    int endcapDetLayerCounter = 0;
    int endcapRadLayerCounter = 0;

    for (int si = 0; si < numSeq; ++si) {
      int seqIdx = (endcapType == 1) ? si : (numSeq - 1 - si);
      const auto& spec        = endcapSeq[seqIdx];
      double      halfT       = spec.thickness / 2.0;
      double      layerZCenter= currentZ + halfT;
      currentZ               += spec.thickness;

      // ---- ENDCAP DETECTOR LAYER ----------------------------------------
      if (spec.type == "detector") {
        int numEndcapLayer = endcapDetLayerCounter++;

        dd4hep::PolyhedraRegular endcapDetectorEnvelope(numSides,
            endcapDetectorLayerInnerRadius, endcapDetectorLayerOuterRadius, 2.0 * endcapDetectorEnvZ);
        std::string endcapDetectorEnvelopeName =
            name + "-EndcapDetectorLayer" + std::to_string(numEndcapLayer + 1);
        dd4hep::Volume endcapDetectorEnvVol(endcapDetectorEnvelopeName, endcapDetectorEnvelope, mat);

        int detElementID = numEndcapLayer + numEndcapDetectorLayers + 1;
        dd4hep::PlacedVolume endcapDetectorEnvelopePhys = endcapVolume.placeVolume(
            endcapDetectorEnvVol,
            dd4hep::Transform3D(dd4hep::RotationZ(0.), dd4hep::Position(0., 0., layerZCenter)));
        dd4hep::DetElement endcapDetectorEnvelopeDE(EndcapDE, endcapDetectorEnvelopeName + "DE", detElementID);
        endcapDetectorEnvelopeDE.setPlacement(endcapDetectorEnvelopePhys);
        endcapDetectorEnvelopePhys.addPhysVolID("layer", numEndcapLayer);
        endcapDetectorEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

        // ---- Endcap sides -----------------------------------------------
        for (int side = 0; side < numSides; ++side) {
          int sideID = (numEndcapLayer + 1) * 100 + (side + 1);

          dd4hep::Trapezoid endcapDetectorSideTrap(
              detectorVolumeThickness / 2.0, detectorVolumeThickness / 2.0,
              endcapDetectorSideLength / 2.0, endcapDetectorSideTrapLength / 2.0,
              endcapDetectorSideTrapYLength / 2.0);
          dd4hep::Box endcapDetectorSideBox(
              detectorVolumeThickness / 2.0, endcapDetectorSideBoxLength / 2.0, endcapDetectorSideBoxYLength / 2.0);
          dd4hep::UnionSolid endcapDetectorSideEnvelope(
              endcapDetectorSideTrap, endcapDetectorSideBox,
              dd4hep::Transform3D(dd4hep::Rotation3D(dd4hep::RotationY(0.)),
                                  dd4hep::Position(0., 0., endcapDetectorYLength / 2.0)));

          std::string endcapDetectorSideEnvName = dd4hep::xml::_toString(sideID, "endcapDetectorSideEnv%d");
          dd4hep::Volume endcapDetectorSideEnvVol(endcapDetectorSideEnvName, endcapDetectorSideEnvelope, mat);

          double angle_degrees = shapeAngle * side;
          double angle_radians = angle_degrees * M_PI / 180.0;
          double endcapDetectorTrapCenterRadius =
              endcapDetectorLayerInnerRadius + (endcapDetectorSideTrapYLength / 2.0);
          double ecXOffset = endcapDetectorTrapCenterRadius * std::cos(angle_radians + shapeAngle_radians);
          double ecYOffset = endcapDetectorTrapCenterRadius * std::sin(angle_radians + shapeAngle_radians);

          dd4hep::RotationZ  ecRotZ(angle_radians + shapeAngle_radians);
          dd4hep::Rotation3D ecRot(ecRotZ);

          dd4hep::Position endcapDetectorSideEnvTrans;
          if (sideLengthEstimate <= (2.0 * dimensions.y())) {
            double z1 = -endcapDetectorEnvZ / 4.0;
            double z2 = -endcapDetectorEnvZ * 3.0 / 4.0;
            double z3 =  endcapDetectorEnvZ / 4.0;
            double z4 =  endcapDetectorEnvZ * 3.0 / 4.0;
            if      (side % 4 == 0) endcapDetectorSideEnvTrans = dd4hep::Position(ecXOffset, ecYOffset, z1);
            else if (side % 4 == 1) endcapDetectorSideEnvTrans = dd4hep::Position(ecXOffset, ecYOffset, z2);
            else if (side % 4 == 2) endcapDetectorSideEnvTrans = dd4hep::Position(ecXOffset, ecYOffset, z3);
            else                    endcapDetectorSideEnvTrans = dd4hep::Position(ecXOffset, ecYOffset, z4);
          } else {
            double zPos = (side % 2 == 0) ? -endcapDetectorEnvZ / 2.0 : endcapDetectorEnvZ / 2.0;
            endcapDetectorSideEnvTrans = dd4hep::Position(ecXOffset, ecYOffset, zPos);
          }

          dd4hep::PlacedVolume endcapDetectorSideEnvPhys = endcapDetectorEnvVol.placeVolume(
              endcapDetectorSideEnvVol,
              dd4hep::Transform3D(ecRot * dd4hep::RotationY(90.0 * dd4hep::degree), endcapDetectorSideEnvTrans));
          dd4hep::DetElement endcapDetectorSideEnvDE(
              endcapDetectorEnvelopeDE, endcapDetectorSideEnvName + "DE", sideID);
          endcapDetectorSideEnvDE.setPlacement(endcapDetectorSideEnvPhys);
          endcapDetectorSideEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

          // ---- Endcap rectangles ----------------------------------------
          int numRectangles = (int)(endcapDetectorYLength / (2.0 * dimensions.z() - overlapZ));
          for (int rectangle = 0; rectangle <= numRectangles; ++rectangle) {
            double rectangleEnvY =
                (endcapDetectorLayerOuterRadius - rectangle * (2.0 * dimensions.y() - overlapY))
                * std::tan(shapeAngle_radians);
            double rectangleEnvZ;
            if (rectangle == numRectangles)
              rectangleEnvZ = endcapRemainderZ * dimensions.z() + clearance / 4.0;
            else
              rectangleEnvZ = dimensions.z() + clearance / 2.0;
            double rectangleEnvX = (rectangleEnvY <= dimensions.y()) ? dimensions.x() : detectorVolumeThickness / 4.5;

            std::string rectangleEnvelopeName = (rectangle == numRectangles)
                ? dd4hep::xml::_toString(rectangle, "rectangleRemainderEnvelope%d")
                : dd4hep::xml::_toString(rectangle, "rectangleEnvelope%d");
            dd4hep::Box    rectangleEnvelope(rectangleEnvX, rectangleEnvY, rectangleEnvZ);
            dd4hep::Volume rectangleEnvVol(rectangleEnvelopeName, rectangleEnvelope, mat);

            double rectangleEnvZPos = (rectangle == numRectangles)
                ? endcapDetectorYLength / 2.0 + ((1.0 - endcapRemainderZ) * dimensions.z())
                  - (rectangle * (2.0 * dimensions.z() - overlapZ)) - clearance
                : endcapDetectorYLength / 2.0 - (rectangle * (2.0 * dimensions.y() - overlapY)) - clearance;

            double yRotation = std::atan(rectangleEnvX / (rectangleEnvZ - 2.0 * overlapY));
            dd4hep::PlacedVolume rectangleEnvelopePhys = endcapDetectorSideEnvVol.placeVolume(
                rectangleEnvVol,
                dd4hep::Transform3D(dd4hep::RotationY(yRotation),
                                    dd4hep::Position(0., 0., rectangleEnvZPos)));
            dd4hep::DetElement rectangleEnvelopeDE(
                endcapDetectorSideEnvDE, rectangleEnvelopeName + "DE", rectangle);
            rectangleEnvelopeDE.setPlacement(rectangleEnvelopePhys);
            rectangleEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

            // ---- Endcap chambers ----------------------------------------
            int numChambersInRectangle = (rectangleEnvY <= dimensions.y())
                ? 0 : (int)(2.0 * rectangleEnvY / (2.0 * dimensions.y() - overlapY));

            for (int chamberIndex = 0; chamberIndex <= numChambersInRectangle; ++chamberIndex) {
              std::stringstream ens;
              ens << "-MuRWELL_Endcap_" << endcapIdCounter++;
              std::string EndcapChamberName = name + ens.str();

              double rectangleRemainderY;
              if (numChambersInRectangle == 0)
                rectangleRemainderY = std::abs(std::fmod((2.0 * rectangleEnvY - clearance),
                    (2.0 * dimensions.y() - clearance))) / (2.0 * dimensions.y());
              else
                rectangleRemainderY = std::fmod(2.0 * (rectangleEnvY - clearance),
                    (2.0 * dimensions.y() - overlapY)) / (2.0 * dimensions.y());

              double chamberHalfZ = (rectangle == numRectangles) ? endcapRemainderZ * dimensions.z() : dimensions.z();
              dd4hep::Box    envelope(dimensions.x(), dimensions.y(), chamberHalfZ);
              dd4hep::Volume envVolume(EndcapChamberName, envelope, lcdd.material(dimensions.materialStr()));

              dd4hep::Box    rrmYEnv(dimensions.x(), rectangleRemainderY * dimensions.y(), chamberHalfZ);
              dd4hep::Volume rectangleRemainderYEnvVolume(
                  EndcapChamberName + "rectangleRemainderY", rrmYEnv, lcdd.material(dimensions.materialStr()));

              double envYPos = (chamberIndex * 2.0 * dimensions.y()) - (overlapY * chamberIndex)
                             + dimensions.y_offset() - rectangleEnvY + dimensions.y() + 0.005;
              double rectangleRemainderREnvYPos = (chamberIndex * 2.0 * dimensions.y()) - (overlapY * chamberIndex)
                             + dimensions.y_offset() - rectangleEnvY
                             + rectangleRemainderY * dimensions.y() + 0.005;

              double zRotation = 0.0, rrmZRotation = 0.0;
              if (numChambersInRectangle != 0) {
                zRotation    = std::atan(dimensions.x() / (dimensions.y() - 2.0 * overlapY));
                rrmZRotation = std::atan(dimensions.x() / (rectangleRemainderY * dimensions.y() - 2.0 * overlapY));
              }
              dd4hep::RotationZ chamberRotation(zRotation);
              dd4hep::RotationZ rrmRotationZ(rrmZRotation);

              auto Slices    = xmlElement.children(_Unicode(slice));
              auto numSlices = xmlElement.numChildren(_Unicode(slice), true);
              dd4hep::xml::Handle_t slice(Slices.reset());
              int sensitiveSliceIndex = 0;

              if (chamberIndex == numChambersInRectangle) {
                // remainder chamber
                dd4hep::PlacedVolume rrmPhys = rectangleEnvVol.placeVolume(
                    rectangleRemainderYEnvVolume,
                    dd4hep::Transform3D(rrmRotationZ,
                                        dd4hep::Position(dimensions.x_offset(), rectangleRemainderREnvYPos, 0.0)));
                rrmPhys.addPhysVolID("chamber", endcapIdCounter);
                dd4hep::DetElement rrmDE(rectangleEnvelopeDE, EndcapChamberName, endcapIdCounter);
                rrmDE.setPlacement(rrmPhys);
                rectangleRemainderYEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());

                double sliceXOffset = -dimensions.x();
                for (unsigned sliceIdx = 0; sliceIdx < numSlices; ++sliceIdx) {
                  dd4hep::xml::DetElement sliceDet = static_cast<dd4hep::xml::DetElement>(slice);
                  dd4hep::Box    sliceShape(sliceDet.x(), rectangleRemainderY * dimensions.y(), chamberHalfZ);
                  std::string    sliceName = dd4hep::xml::_toString(sliceIdx, "slice%d");
                  dd4hep::Volume sliceVolume(sliceName, sliceShape, lcdd.material(slice.attr<std::string>("material")));
                  dd4hep::PlacedVolume slicePV = rectangleRemainderYEnvVolume.placeVolume(
                      sliceVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.),
                                                       dd4hep::Position(sliceXOffset + sliceDet.x(), 0.0, 0.0)));
                  if (slice.hasAttr("vis"))       sliceVolume.setVisAttributes(lcdd, sliceDet.visStr());
                  if (slice.hasAttr("sensitive") && sliceDet.isSensitive()) {
                    dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                    sensDet.setType(sdType.typeStr());
                    sliceVolume.setSensitiveDetector(sensDet);
                    slicePV.addPhysVolID("slice", sensitiveSliceIndex);
                    dd4hep::DetElement sliceDE(rrmDE, "slice_" + std::to_string(sensitiveSliceIndex), sensitiveSliceIndex);
                    sliceDE.setPlacement(slicePV);
                    ++sensitiveSliceIndex;
                  }
                  sliceXOffset += 2.0 * sliceDet.x();
                  slice.m_node = Slices.next();
                }

              } else {
                // full chamber
                dd4hep::PlacedVolume envPhys = rectangleEnvVol.placeVolume(
                    envVolume, dd4hep::Transform3D(chamberRotation,
                                                   dd4hep::Position(dimensions.x_offset(), envYPos, 0.0)));
                envPhys.addPhysVolID("chamber", endcapIdCounter);
                dd4hep::DetElement envDE(rectangleEnvelopeDE, EndcapChamberName, endcapIdCounter);
                envDE.setPlacement(envPhys);
                envVolume.setVisAttributes(lcdd, xmlDet.visStr());

                double sliceXOffset = -dimensions.x();
                for (unsigned sliceIdx = 0; sliceIdx < numSlices; ++sliceIdx) {
                  dd4hep::xml::DetElement sliceDet = static_cast<dd4hep::xml::DetElement>(slice);
                  dd4hep::Box    sliceShape(sliceDet.x(), dimensions.y(), chamberHalfZ);
                  std::string    sliceName = dd4hep::xml::_toString(sliceIdx, "slice%d");
                  dd4hep::Volume sliceVolume(sliceName, sliceShape, lcdd.material(slice.attr<std::string>("material")));
                  dd4hep::PlacedVolume slicePV = envVolume.placeVolume(
                      sliceVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.),
                                                       dd4hep::Position(sliceXOffset + sliceDet.x(), 0.0, 0.0)));
                  if (slice.hasAttr("vis"))       sliceVolume.setVisAttributes(lcdd, sliceDet.visStr());
                  if (slice.hasAttr("sensitive") && sliceDet.isSensitive()) {
                    dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                    sensDet.setType(sdType.typeStr());
                    sliceVolume.setSensitiveDetector(sensDet);
                    slicePV.addPhysVolID("slice", sensitiveSliceIndex);
                    dd4hep::DetElement sliceDE(envDE, "slice_" + std::to_string(sensitiveSliceIndex), sensitiveSliceIndex);
                    sliceDE.setPlacement(slicePV);
                    ++sensitiveSliceIndex;
                  }
                  sliceXOffset += 2.0 * sliceDet.x();
                  slice.m_node = Slices.next();
                }
              }
            } // chamberIndex
          }   // rectangle
        }     // side

      // ---- ENDCAP RADIATOR LAYER ----------------------------------------
      } else {
        int numEndcapRadiatorLayer = endcapRadLayerCounter++;

        dd4hep::PolyhedraRegular endcapRadiatorEnvelope(numSides,
            endcapRadiatorLayerInnerRadius, endcapRadiatorLayerOuterRadius, spec.thickness);
        std::string endcapRadiatorEnvelopeName =
            name + "-EndcapRadiatorLayer" + std::to_string(numEndcapRadiatorLayer + 1);
        dd4hep::Volume endcapRadiatorEnvVol(endcapRadiatorEnvelopeName, endcapRadiatorEnvelope, mat);

        dd4hep::PlacedVolume endcapRadiatorEnvelopePhys = endcapVolume.placeVolume(
            endcapRadiatorEnvVol,
            dd4hep::Transform3D(dd4hep::RotationZ(0.), dd4hep::Position(0., 0., layerZCenter)));
        dd4hep::DetElement endcapRadiatorEnvelopeDE(
            EndcapDE, endcapRadiatorEnvelopeName + "DE", numEndcapRadiatorLayer + 1);
        endcapRadiatorEnvelopeDE.setPlacement(endcapRadiatorEnvelopePhys);
        endcapRadiatorEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

        for (int side = 0; side < numSides; ++side) {
          int sideID = (numEndcapRadiatorLayer + 1) * 100 + (side + 1);
          dd4hep::Trapezoid endcapRadiator(
              spec.thickness / 2.0, spec.thickness / 2.0,
              endcapRadiatorSideLength / 2.0, endcapRadiatorSideLength2 / 2.0,
              endcapYLength / 2.0);
          std::string    endcapRadiatorName = dd4hep::xml::_toString(sideID, "endcapRadiator%d");
          dd4hep::Volume endcapRadiatorVol(endcapRadiatorName, endcapRadiator, spec.matObj);

          double angle_degrees = shapeAngle * side;
          double angle_radians = angle_degrees * M_PI / 180.0;
          double ecRadMidRadius = endcapRadiatorLayerInnerRadius + endcapYLength / 2.0;
          double ecRXOffset     = ecRadMidRadius * std::cos(angle_radians + shapeAngle_radians);
          double ecRYOffset     = ecRadMidRadius * std::sin(angle_radians + shapeAngle_radians);

          dd4hep::RotationZ  ecRRotZ(angle_radians + shapeAngle_radians);
          dd4hep::Rotation3D ecRRot(ecRRotZ);

          dd4hep::PlacedVolume endcapRadiatorPhys = endcapRadiatorEnvVol.placeVolume(
              endcapRadiatorVol,
              dd4hep::Transform3D(ecRRot * dd4hep::RotationY(90.0 * dd4hep::degree),
                                  dd4hep::Position(ecRXOffset, ecRYOffset, 0.0)));
          dd4hep::DetElement endcapRadiatorDE(endcapRadiatorEnvelopeDE, endcapRadiatorName + "DE", sideID);
          endcapRadiatorDE.setPlacement(endcapRadiatorPhys);
          endcapRadiatorVol.setVisAttributes(lcdd.visAttributes("yoke_vis"));
        }
      }
    } // endcap sequence
  }   // endcapType

  return detElement;
}

DECLARE_DETELEMENT(muonSystem_o1_v02, createmuonSystem_o1_v02)
