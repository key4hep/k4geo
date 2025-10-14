/*
 @author: Mahmoud Ali
mahmoud.ali@cern.ch

Factory for IDEA muon system
Expected xml structure (the 'sensitive' keyword is optional and defaults to false):
<detector type="muonSystemMuRWELL_o1_v01" ...>
  <dimensions x="..." y="..." z="..." z_offset="..." x_offset="..." y_offset="...">  <!--  dimension of the local
chamber envelope. x: the half length of the thickness of the chamber, y&z: the half length of the 2D plane dimensions of
the chambers--> <sensitive type="tracker"/>

  <!-- Specify the detector parameters and the overlap // if you want exclude any component, e.g: endcap, just put
Endcap =0 // radius is put in the middle, so its not the inner neither the outer  --> <generalParameters numSides="..."
overlapY="..." overlapZ="..." clearance="..."/> <Barrel numDetectorLayers ="..."  rmin="..." length="..." numYokes="..."
yoke_Thickness="..." yoke_Material="..."/> <Endcap numDetectorLayers="..."  rmin="..." rmax="..." z_offset="..."
numYokes="..." yoke_Thickness="..." yoke_Material="..." />

  <slice x="..."  material="...">
  . . . .
  <slice x="..."  material="..." sensitive="true">
</detector>

If used with sensitive layers, the readout must contain a "slice" field
*/

#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/DocumentHandler.h"
#include "XML/Utilities.h"
#include "XML/XMLElements.h"
#include <cmath>
#include <sstream>

using namespace std;
using namespace dd4hep;

using dd4hep::rec::SurfaceType;
using dd4hep::rec::Vector3D;
using dd4hep::rec::VolCylinder;
using dd4hep::rec::VolPlane;
using dd4hep::rec::volSurfaceList;

static dd4hep::Ref_t createmuonSystemMuRWELL_o1_v01(dd4hep::Detector& lcdd, dd4hep::xml::Handle_t xmlElement,
                                                    dd4hep::SensitiveDetector sensDet) {
  dd4hep::xml::DetElement xmlDet = static_cast<dd4hep::xml::DetElement>(xmlElement);
  std::string name = xmlDet.nameStr();
  dd4hep::DetElement detElement(name, xmlDet.id());
  dd4hep::Material mat = lcdd.material("Air");
  dd4hep::Volume experimentalHall = lcdd.pickMotherVolume(detElement);

  xml_comp_t dimensions(xmlDet.dimensions());

  // ----------------------------------------------------------------------------------------------------
  //                         --- General parameters ---

  auto generalParameters = xmlElement.child(_Unicode(generalParameters));
  int numSides = generalParameters.attr<int>("numSides");
  double overlapY = generalParameters.attr<double>("overlapY");
  double overlapZ = generalParameters.attr<double>("overlapZ");
  double clearance = generalParameters.attr<double>(
      "clearance"); // it's a small distance to be used to avoid overlapping between the different volumes ~ 1 mm

  //                         --- Barrel parameters ---

  auto Barrel = xmlElement.child(_Unicode(Barrel));
  int numBarrelDetectorLayers = Barrel.attr<double>("numDetectorLayers");
  double radius = Barrel.attr<double>("rmin");
  double barrelLength = Barrel.attr<double>("length");
  int numBarrelRadiators = Barrel.attr<double>("numYokes");
  double barrelRadiatorThickness = Barrel.attr<double>("yoke_Thickness");
  dd4hep::Material barrelRadiatorMaterial = lcdd.material(Barrel.attr<std::string>("yoke_Material"));

  //                    --- Endcap parameters ---

  auto Endcap = xmlElement.child(_Unicode(Endcap));
  int numEndcapDetectorLayers = Endcap.attr<double>("numDetectorLayers");
  double endcapDetectorLayerInnerRadius = Endcap.attr<double>("rmin");
  double endcapZOffset = Endcap.attr<double>("z_offset");
  int numEndcapRadiators = Endcap.attr<double>("numYokes");
  double endcapRadiatorThickness = Endcap.attr<double>("yoke_Thickness");
  double endcapRadiatorLayerInnerRadius = endcapDetectorLayerInnerRadius;
  dd4hep::Material endcapRadiatorMaterial = lcdd.material(Endcap.attr<std::string>("yoke_Material"));

  // -----------------------------------------------------------------------------------------------------------

  double shapeAngle = 360.0 / numSides;        // it's the full angle
  double shapeAngle_radians = M_PI / numSides; // Isn't it a half angle!!
  double angle_clearance = 0.0; // 0.02 works good but needs the detector volume thick to be more than 60 mm. // it's
                                // less than 1 degree, we use the clearnce to avoid overlapping

  double sideLengthEstimate = 2 * (radius)*std::tan(shapeAngle_radians);

  double chamberAngle_rad = std::atan(dimensions.x() / dimensions.y());
  double rectangleThickness = (2 * dimensions.y()) * std::sin(chamberAngle_rad) +
                              (2 * dimensions.x()) * std::cos(chamberAngle_rad) + 1.2 * clearance;
  double rectangleAngle_rad = std::atan(rectangleThickness / dimensions.z());
  double detectorVolumeThickness;
  double endcapDetectorEnvZ; // endcap detetcor layer thickness
  if (sideLengthEstimate <= (2 * dimensions.y()) &&
      numBarrelDetectorLayers == 1) { // putting this condtition to minimize the volume thickness in case of there is
                                      // only one chamber per side(espicially for circular shape detector)
    detectorVolumeThickness =
        rectangleThickness *
        1.33; // multiplying by 1.33 fit the best thickness avoiding overlaps with smaller volumes insides (rectangles).
    endcapDetectorEnvZ = 2 * detectorVolumeThickness;
  } else {
    detectorVolumeThickness =
        (2 * dimensions.z()) * std::sin(rectangleAngle_rad) + rectangleThickness * std::cos(rectangleAngle_rad);
    endcapDetectorEnvZ = detectorVolumeThickness;
  }

  // Automation Endcap R-max to be enclosed by Barrel last layer
  double endcapDetectorLayerOuterRadius = radius + (numBarrelDetectorLayers - 1) * (2 * detectorVolumeThickness) +
                                          numBarrelRadiators * barrelRadiatorThickness;
  double endcapRadiatorLayerOuterRadius = endcapDetectorLayerOuterRadius;

  double endcapDetectorSideLength =
      (2 * (endcapDetectorLayerInnerRadius + 2 * dimensions.y()) * std::tan(shapeAngle_radians)) + 2 * clearance;
  double endcapDetectorSideTrapLength =
      (2 * (endcapDetectorLayerOuterRadius)*std::tan(shapeAngle_radians)) + 2 * clearance;

  double endcapDetectorSideTrapYLength =
      endcapDetectorLayerOuterRadius - 2 * dimensions.z() - endcapDetectorLayerInnerRadius;
  double endcapDetectorSideBoxLength = 2 * (endcapDetectorLayerOuterRadius)*std::tan(shapeAngle_radians);
  double endcapDetectorSideBoxYLength = 2 * dimensions.z();

  // endcapRadiatorSideLength and endcapRadiatorSideLength2 are the lengths of the parallel sides of the trapezoid.
  double endcapRadiatorSideLength = 2 * (endcapRadiatorLayerInnerRadius)*std::tan(
                                            shapeAngle_radians); // it is also the same for the endcap detector layers.
  double endcapRadiatorSideLength2 = 2 * (endcapRadiatorLayerOuterRadius)*std::tan(
                                             shapeAngle_radians); // it is also the same for the endcap detector layers.

  double endcapDetectorYLength = endcapDetectorLayerOuterRadius - endcapDetectorLayerInnerRadius;
  double endcapYLength = endcapRadiatorLayerOuterRadius -
                         endcapRadiatorLayerInnerRadius; // It is the distance betwwen the inner and the outer radius of
                                                         // the endcap, it can be in both Y and X dimensions //it is
                                                         // also the same for the endcap detector layers.

  double endcapRemainderZ =
      std::fmod((endcapDetectorYLength - 2 * clearance), (2 * dimensions.z() - overlapZ)) / (2 * dimensions.z()) -
      (2 * clearance / dimensions.z());

  double barrelRMax =
      radius + numBarrelDetectorLayers * (2 * detectorVolumeThickness) + numBarrelRadiators * barrelRadiatorThickness;
  double barreltotalLength =
      barrelLength + (numEndcapDetectorLayers * 2) * (2 * detectorVolumeThickness) +
      (numEndcapRadiators * 2) * endcapRadiatorThickness; // This condition to make the last barrel layer encloses all
                                                          // the endcap layers inside it.
  double EndcaptotalLength =
      numEndcapDetectorLayers * (2 * endcapDetectorEnvZ) + numEndcapRadiators * endcapRadiatorThickness;
  double endcapOffset = endcapZOffset + EndcaptotalLength / 2.0;

  int barrelIdCounter = 1;
  int endcapIdCounter = 1;

  //-------------------------// Building system envelope //----------------------------

  dd4hep::PolyhedraRegular BarrelEnv(numSides, radius, barrelRMax, barreltotalLength);
  dd4hep::PolyhedraRegular EndcapEnv(numSides, endcapDetectorLayerInnerRadius, endcapDetectorLayerOuterRadius,
                                     EndcaptotalLength);

  double unionOffsetZpositive = endcapOffset;
  double unionOffsetZnegative = -endcapOffset;

  dd4hep::Position unionPos(0.0, 0.0, unionOffsetZpositive);
  dd4hep::Position unionPos2(0.0, 0.0, unionOffsetZnegative);
  dd4hep::Rotation3D unionRot(dd4hep::RotationY(0.0 * dd4hep::degree));

  dd4hep::Transform3D unionTransform(unionRot, unionPos);
  dd4hep::Transform3D unionTransform2(unionRot, unionPos2);

  // Combining two shapes by UnionSolid: the first shape is centralized and the second transform around the first..
  dd4hep::UnionSolid barrelAndPositiveEndcap(BarrelEnv, EndcapEnv, unionTransform);
  dd4hep::UnionSolid systemEnvelope(barrelAndPositiveEndcap, EndcapEnv, unionTransform2);
  dd4hep::Volume detectorVolume(name, systemEnvelope, mat);

  dd4hep::Position detectorTrans(0., 0., 0.);
  dd4hep::PlacedVolume detectorPhys = experimentalHall.placeVolume(
      detectorVolume, dd4hep::Transform3D(dd4hep::RotationZ(shapeAngle_radians), detectorTrans));
  detectorPhys.addPhysVolID("system", xmlDet.id());
  detElement.setPlacement(detectorPhys);
  detectorVolume.setVisAttributes(lcdd.visAttributes("no_vis"));

  // ----------------------------------------------------------------------------------------------------
  // ------------------------------// B A R R E L // ----------------------------------------------------

  dd4hep::PolyhedraRegular BarrelEnvWithoutLastLayer(numSides, radius, barrelRMax, barrelLength);
  dd4hep::PolyhedraRegular BarrelLastLayerEnv(numSides, (barrelRMax - 2 * detectorVolumeThickness), barrelRMax,
                                              barreltotalLength);
  std::string barrelName = name + "-Barrel" + std::to_string(0);

  dd4hep::Position barrelUnionPos(0.0, 0.0, 0.0);
  dd4hep::Rotation3D barrelUnionRot(dd4hep::RotationY(0.0 * dd4hep::degree));
  dd4hep::Transform3D barrelUnionTransform(barrelUnionRot, barrelUnionPos);

  dd4hep::UnionSolid barrelUnion(BarrelEnvWithoutLastLayer, BarrelLastLayerEnv, barrelUnionTransform);
  dd4hep::Volume BarrelVolume(barrelName, barrelUnion, mat);

  dd4hep::Position barrelTrans(0., 0., 0.);
  dd4hep::PlacedVolume barrelPhys =
      detectorVolume.placeVolume(BarrelVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), barrelTrans));
  barrelPhys.addPhysVolID("type", 0);
  dd4hep::DetElement BarrelDE(detElement, name + "-BarrelDE", 0);
  BarrelDE.setPlacement(barrelPhys);
  BarrelVolume.setVisAttributes(lcdd.visAttributes("no_vis"));

  // ---------- loop to creat the barrel layers -----------

  for (int numBarrelLayer = 0; numBarrelLayer < numBarrelDetectorLayers; ++numBarrelLayer) {

    double barrelLayerRMin =
        radius + numBarrelLayer * (2 * detectorVolumeThickness) +
        numBarrelLayer * barrelRadiatorThickness; // Automation of inner radius of different layers, taking into account
                                                  // that every detector layer is followed by a yoke(radiator) layer.
    double barrelLayerRMax = barrelLayerRMin + 2 * detectorVolumeThickness;
    double barrelLayerRMid = (barrelLayerRMin + barrelLayerRMax) / 2.0;
    double barrelLayerLength;
    if (numBarrelLayer == numBarrelDetectorLayers - 1) {
      barrelLayerLength =
          barrelLength + (numEndcapDetectorLayers * 2) * (2 * detectorVolumeThickness) +
          (numEndcapRadiators * 2) * endcapRadiatorThickness; // This condition to make the last barrel layer encloses
                                                              // all the endcap layers inside it.
    } else {
      barrelLayerLength = barrelLength;
    }

    dd4hep::PolyhedraRegular BarrelDetectorLayer(numSides, barrelLayerRMin, barrelLayerRMax, barrelLayerLength);
    std::string barrelDetectorName = name + "-BarrelDetectorLayer" + std::to_string(numBarrelLayer + 1);
    dd4hep::Volume BarrelDetectorLayerVolume(barrelDetectorName, BarrelDetectorLayer, mat);

    // sideLength and sideLength2 are the lengths of the parallel sides of the trapezoid.
    double sideLength = 2 * (barrelLayerRMid - detectorVolumeThickness / 2.0) * std::tan(shapeAngle_radians);
    double sideLength2 = 2 * (barrelLayerRMid + detectorVolumeThickness / 2.0) * std::tan(shapeAngle_radians);
    double sideEnvX = detectorVolumeThickness / 2.0;
    double sideEnvY = (sideLength / 2.0);
    double sideEnvY2 = (sideLength2 / 2.0);
    double sideEnvZ = (barrelLayerLength / 2.0);
    double remainderZ =
        std::fmod((barrelLayerLength - 2 * clearance), (2 * dimensions.z() - overlapZ)) / (2 * dimensions.z()) -
        (2 * clearance / dimensions.z());

    dd4hep::Position detectorLayerTrans(0., 0., 0.);
    dd4hep::PlacedVolume detectorLayerPhys = BarrelVolume.placeVolume(
        BarrelDetectorLayerVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), detectorLayerTrans));
    detectorLayerPhys.addPhysVolID("layer", numBarrelLayer);
    dd4hep::DetElement BarrelDetectorLayerDE(
        BarrelDE, name + "-Barrel_DetectorLayerDE_" + std::to_string(numBarrelLayer + 1), numBarrelLayer + 1);
    BarrelDetectorLayerDE.setPlacement(detectorLayerPhys);
    BarrelDetectorLayerVolume.setVisAttributes(lcdd.visAttributes("no_vis"));

    // ------- Dividing every layer into # of sides = # of polyhedron sides -----------
    for (int side = 0; side < numSides; ++side) {

      int sideID = (numBarrelLayer + 1) * 100 + (side + 1); // to differentiated with the same side in different layers.
      dd4hep::Box sideEnvelope(sideEnvX, sideEnvY, sideEnvZ);
      dd4hep::Box sideEnvelope2(sideEnvX, sideEnvY2, sideEnvZ);
      std::string sideName = dd4hep::xml::_toString(sideID, "side%d");
      dd4hep::Volume sideVol(sideName, sideEnvelope, mat);
      dd4hep::Volume sideVol2(sideName, sideEnvelope2, mat);

      double angle_degrees = shapeAngle * side; // Calculate the angle for each side
      double angle_radians = angle_degrees * M_PI / 180.0;

      double sideXOffset =
          (barrelLayerRMid - detectorVolumeThickness / 2.0) * std::cos(angle_radians + shapeAngle_radians);
      double sideYOffset =
          (barrelLayerRMid - detectorVolumeThickness / 2.0) * std::sin(angle_radians + shapeAngle_radians);

      double sideXOffset2 =
          (barrelLayerRMid + detectorVolumeThickness / 2.0) * std::cos(angle_radians + shapeAngle_radians);
      double sideYOffset2 =
          (barrelLayerRMid + detectorVolumeThickness / 2.0) * std::sin(angle_radians + shapeAngle_radians);

      dd4hep::RotationZ sideRotationZ(angle_radians + shapeAngle_radians + angle_clearance);
      dd4hep::Rotation3D sideRotation = dd4hep::Rotation3D(sideRotationZ);

      double sideXPos = sideXOffset;
      double sideYPos = sideYOffset;
      double sideZPos = 0.0;
      double sideXPos2 = sideXOffset2;
      double sideYPos2 = sideYOffset2;

      dd4hep::Position sideTrans;
      dd4hep::PlacedVolume sidePhys;
      dd4hep::DetElement sideDE;

      if (side % 2 == 0) {
        sideTrans = dd4hep::Position(sideXPos, sideYPos, sideZPos);
        sidePhys = BarrelDetectorLayerVolume.placeVolume(sideVol, dd4hep::Transform3D(sideRotation, sideTrans));
      } else {
        sideTrans = dd4hep::Position(sideXPos2, sideYPos2, sideZPos);
        sidePhys = BarrelDetectorLayerVolume.placeVolume(sideVol2, dd4hep::Transform3D(sideRotation, sideTrans));
      }

      sideDE = dd4hep::DetElement(BarrelDetectorLayerDE, sideName + "DE", sideID);
      sideDE.setPlacement(sidePhys);
      sideVol.setVisAttributes(lcdd, xmlDet.visStr());
      sideVol2.setVisAttributes(lcdd, xmlDet.visStr());

      // -------------------------------------------------------------------------------------------------------
      //  Dividing every side to small rectangles  //
      // -------------------------------------------------------------------------------------------------------

      int numRectangles =
          barrelLayerLength /
          (2 * dimensions.z() - overlapZ); // numbers of the rectangles in each side.. depends on the number of the
                                           // chambers that can be overlapped in Z direction

      for (int rectangle = 0; rectangle < (numRectangles + 1); ++rectangle) {

        double rectangleEnvY;
        if (side % 2 == 0) {
          rectangleEnvY = sideLength / 2.0; // + clearance;
        } else {
          rectangleEnvY = sideLength2 / 2.0; // + clearance;
        }
        double rectangleEnvZ = dimensions.z() + clearance / 2.0;
        double rectangleRemainderEnvZ =
            (remainderZ * dimensions.z()) +
            clearance / 4.0; // this is the last rectangle in the side, and it is smaller than the usual one, with
                             // larger rotation angle. so it need to be shorter in length to avoid overlap.
        double rectangleEnvX;
        if (rectangleEnvY <= dimensions.y()) {
          rectangleEnvX = dimensions.x(); // in case if there is only one chamber inside the rectangle.
        } else {
          rectangleEnvX = detectorVolumeThickness /
                          4.5; // dividing by 4.5 gets the best thickness for the recatngle to avoid any overlap ~ in
                               // our case the uRWELL the best rectangle thickness is 26.667 mm which is 120/4.5.
        }

        dd4hep::Box rectangleEnvelope(rectangleEnvX, rectangleEnvY, rectangleEnvZ);
        std::string rectangleEnvelopeName = dd4hep::xml::_toString(rectangle, "rectangleEnvelope%d");
        dd4hep::Volume rectangleEnvVol(rectangleEnvelopeName, rectangleEnvelope, mat);

        dd4hep::Box rectangleRemainderEnvelope(rectangleEnvX, rectangleEnvY, rectangleRemainderEnvZ);
        std::string rectangleRemainderEnvelopeName = dd4hep::xml::_toString(rectangle, "rectangleRemainderEnvelope%d");
        dd4hep::Volume rectangleRemainderEnvVol(rectangleRemainderEnvelopeName, rectangleRemainderEnvelope, mat);

        double rectangleEnvXPos = 0.0;
        double rectangleEnvYPos = 0.0;
        double rectangleEnvZPos =
            barrelLayerLength / 2.0 - dimensions.z() - (rectangle * (2 * dimensions.z() - overlapZ)) - clearance;
        double rectangleRemainderEnvXPos = sideEnvX / 2.0;
        double rectangleRemainderEnvZPos = barrelLayerLength / 2.0 - (remainderZ * dimensions.z()) -
                                           (rectangle * (2 * dimensions.z() - overlapZ)) - clearance;

        double yRotation = std::atan(rectangleEnvX / (rectangleEnvZ - (2 * overlapZ)));
        double yRemainderRotation = std::atan(rectangleEnvX / (rectangleRemainderEnvZ - (2 * overlapZ)));
        dd4hep::RotationY rotationY(yRotation);
        dd4hep::RotationY remainderRotationY;
        if (rectangleRemainderEnvZ <= rectangleEnvZ / 10.0) {
          remainderRotationY = dd4hep::RotationY(yRotation);
        } else {
          remainderRotationY = dd4hep::RotationY(yRemainderRotation);
        }

        // _________________________________________________________________________________________________
        if (rectangle == numRectangles) {

          dd4hep::Position rectangeEnvelopeTrans;
          if (rectangleRemainderEnvZ <= rectangleEnvZ / 10.0) {
            rectangeEnvelopeTrans =
                dd4hep::Position(rectangleRemainderEnvXPos, rectangleEnvYPos, rectangleRemainderEnvZPos);
          } else {
            rectangeEnvelopeTrans = dd4hep::Position(rectangleEnvXPos, rectangleEnvYPos, rectangleRemainderEnvZPos);
          }
          dd4hep::PlacedVolume rectangleEnvelopePhys;
          if (side % 2 == 0) {
            rectangleEnvelopePhys = sideVol.placeVolume(rectangleRemainderEnvVol,
                                                        dd4hep::Transform3D(remainderRotationY, rectangeEnvelopeTrans));
          } else {
            rectangleEnvelopePhys = sideVol2.placeVolume(
                rectangleRemainderEnvVol, dd4hep::Transform3D(remainderRotationY, rectangeEnvelopeTrans));
          }
          dd4hep::DetElement rectangleEnvelopeDE(sideDE, rectangleRemainderEnvelopeName + "DE", rectangle);
          rectangleEnvelopeDE.setPlacement(rectangleEnvelopePhys);
          rectangleEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

          // ------------------------ start to build the chamber envelopes -------------------

          int numChambersInRectangle;
          if (rectangleEnvY <= dimensions.y()) {
            numChambersInRectangle = 0; // number of the chambers in each rectangle
          } else {
            numChambersInRectangle =
                2 * rectangleEnvY / (2 * dimensions.y() - overlapY); // number of the chambers in each rectangle
          }

          for (int chamberIndex = 0; chamberIndex < (numChambersInRectangle + 1); chamberIndex++) {

            std::stringstream barrelNameStream;
            barrelNameStream << "-MuRWELL_Barrel_" << barrelIdCounter++;
            std::string BarrelChamberName = name + barrelNameStream.str();

            dd4hep::Box envelope(dimensions.x(), dimensions.y(), remainderZ * dimensions.z());
            dd4hep::Volume envVolume(BarrelChamberName, envelope, lcdd.material(dimensions.materialStr()));

            double rectangleRemainderY;
            double rectangleRemainderREnvYPos;
            if (numChambersInRectangle == 0) {
              rectangleRemainderY =
                  std::abs(std::fmod((2 * rectangleEnvY - clearance), (2 * dimensions.y() - clearance))) /
                  (2 * dimensions.y());
              rectangleRemainderREnvYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) +
                                           dimensions.y_offset() - rectangleEnvY +
                                           rectangleRemainderY * dimensions.y() + clearance / 20;
            } else {
              rectangleRemainderY =
                  std::fmod(2 * (rectangleEnvY - clearance), (2 * dimensions.y() - overlapY)) / (2 * dimensions.y());
              rectangleRemainderREnvYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) +
                                           dimensions.y_offset() - rectangleEnvY +
                                           rectangleRemainderY * dimensions.y() + 1.5 * clearance;
            }

            dd4hep::Box rectangleRemainderYEnvelope(dimensions.x(), rectangleRemainderY * dimensions.y(),
                                                    remainderZ * dimensions.z());
            dd4hep::Volume rectangleRemainderYEnvVolume(BarrelChamberName + "rectangleRemainderY",
                                                        rectangleRemainderYEnvelope,
                                                        lcdd.material(dimensions.materialStr()));

            double envYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() -
                             rectangleEnvY + dimensions.y() +
                             clearance / 20.0; // found that the positioning of the chambers inside the rectangle had an
                                               // overlap with the mother volume ~ 45 um.
            // double rectangleRemainderREnvYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) +
            // dimensions.y_offset() - rectangleEnvY + rectangleRemainderY * dimensions.y() + 1.5 * clearance;

            double zRotation;
            double rectangleRemainderZRotation;
            if (numChambersInRectangle == 0) {
              zRotation = 0.0;
              rectangleRemainderZRotation = 0.0;
            } else {
              zRotation = std::atan(dimensions.x() / (dimensions.y() - (2 * overlapY)));
              rectangleRemainderZRotation =
                  std::atan(dimensions.x() / (rectangleRemainderY * dimensions.y() -
                                              overlapY * 1.5)); // Y and Z are reversed in local remainder
            }

            dd4hep::RotationZ chamberRotation(zRotation);
            dd4hep::RotationZ rectangleRemainderRotationZ(rectangleRemainderZRotation);

            auto Slices = xmlElement.children(_Unicode(slice));
            auto numSlices = xmlElement.numChildren(_Unicode(slice), true);
            dd4hep::xml::Handle_t slice(Slices.reset());
            int sensitiveSliceIndex = 0;

            // --- two cases : one for remainder chambers ----------------------------------------------

            if (chamberIndex == numChambersInRectangle) {

              dd4hep::Position rectangleRemainderTrans(dimensions.x_offset(), rectangleRemainderREnvYPos, 0.0);
              dd4hep::PlacedVolume rectangleRemainderEnvPhys = rectangleRemainderEnvVol.placeVolume(
                  rectangleRemainderYEnvVolume,
                  dd4hep::Transform3D(rectangleRemainderRotationZ, rectangleRemainderTrans));
              rectangleRemainderEnvPhys.addPhysVolID("chamber", barrelIdCounter);
              dd4hep::DetElement rectangleRemainderEnvDE(rectangleEnvelopeDE, BarrelChamberName, barrelIdCounter);
              rectangleRemainderEnvDE.setPlacement(rectangleRemainderEnvPhys);
              rectangleRemainderYEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());
              // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() <<
              // std::endl;

              double sliceXOffset = -dimensions.x();
              for (unsigned sliceIdx = 0; sliceIdx < numSlices; ++sliceIdx) {
                dd4hep::xml::DetElement sliceDet = static_cast<dd4hep::xml::DetElement>(slice);
                dd4hep::Box sliceShape(
                    sliceDet.x(), rectangleRemainderY * dimensions.y(),
                    remainderZ *
                        dimensions.z()); // I made the y-z dimensions of the slices are the same like the chamber (which
                                         // is normal case in most of the detectors),  if you need change any slice y-z
                                         // dimension replace dimensions.y() with sliceDet.y()
                std::string sliceName = dd4hep::xml::_toString(sliceIdx, "slice%d");
                dd4hep::Volume sliceVolume(sliceName, sliceShape, lcdd.material(slice.attr<std::string>("material")));

                dd4hep::Position transSlice(
                    sliceXOffset + sliceDet.x(), 0.0,
                    0.0); // all the slices are centered in the chamber except in x-axis they are accumalted over each
                          // other respectivley.. If you want to add an offset to any direction by the user you can add
                          // sliceDet.y_offset() instead of 0.
                dd4hep::PlacedVolume slicePlacedVolume = rectangleRemainderYEnvVolume.placeVolume(
                    sliceVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), transSlice));

                if (slice.hasAttr("vis")) {
                  sliceVolume.setVisAttributes(lcdd, sliceDet.visStr());
                }
                if (slice.hasAttr("sensitive") && sliceDet.isSensitive()) {
                  dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                  sensDet.setType(sdType.typeStr());
                  sliceVolume.setSensitiveDetector(sensDet);
                  slicePlacedVolume.addPhysVolID("slice", sensitiveSliceIndex);
                  dd4hep::DetElement sliceDE(rectangleRemainderEnvDE, "slice_" + std::to_string(sensitiveSliceIndex),
                                             sensitiveSliceIndex);
                  sliceDE.setPlacement(slicePlacedVolume);
                  //  dd4hep::printout(dd4hep::INFO,"Sensitive Slice has been created at", name, BarrelChamberName);
                  sensitiveSliceIndex++;
                }
                // Increment the current x-offset by the width of the current slice
                sliceXOffset += (2 * sliceDet.x());
                slice.m_node = Slices.next();
              }

              // ---------------- Second case: for the full chambers
            } else {

              dd4hep::Position trans(dimensions.x_offset(), envYPos, 0.0);
              dd4hep::PlacedVolume envPhys =
                  rectangleRemainderEnvVol.placeVolume(envVolume, dd4hep::Transform3D(chamberRotation, trans));
              envPhys.addPhysVolID("chamber", barrelIdCounter);
              dd4hep::DetElement envDE(rectangleEnvelopeDE, BarrelChamberName, barrelIdCounter);
              envDE.setPlacement(envPhys);
              envVolume.setVisAttributes(lcdd, xmlDet.visStr());

              // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() <<
              // std::endl;
              double sliceXOffset = -dimensions.x();
              for (unsigned sliceIdx = 0; sliceIdx < numSlices; ++sliceIdx) {
                dd4hep::xml::DetElement sliceDet = static_cast<dd4hep::xml::DetElement>(slice);
                dd4hep::Box sliceShape(sliceDet.x(), dimensions.y(), remainderZ * dimensions.z());
                std::string sliceName = dd4hep::xml::_toString(sliceIdx, "slice%d");
                dd4hep::Volume sliceVolume(sliceName, sliceShape, lcdd.material(slice.attr<std::string>("material")));
                dd4hep::Position transSlice(sliceXOffset + sliceDet.x(), 0.0, 0.0);
                dd4hep::PlacedVolume slicePlacedVolume =
                    envVolume.placeVolume(sliceVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), transSlice));

                if (slice.hasAttr("vis")) {
                  sliceVolume.setVisAttributes(lcdd, sliceDet.visStr());
                }
                if (slice.hasAttr("sensitive") && sliceDet.isSensitive()) {
                  dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                  sensDet.setType(sdType.typeStr());
                  sliceVolume.setSensitiveDetector(sensDet);
                  slicePlacedVolume.addPhysVolID("slice", sensitiveSliceIndex);
                  dd4hep::DetElement sliceDE(envDE, "slice_" + std::to_string(sensitiveSliceIndex),
                                             sensitiveSliceIndex);
                  sliceDE.setPlacement(slicePlacedVolume);
                  // dd4hep::printout(dd4hep::INFO,"Sensitive slice has been created at", name, BarrelChamberName);
                  sensitiveSliceIndex++;
                }
                // Increment the current x-offset by the width of the current slice
                sliceXOffset += (2 * sliceDet.x());
                slice.m_node = Slices.next();
              }
            }
          }

        } else {
          // ......_____________________________________________________________________________________________

          dd4hep::Position rectangeEnvelopeTrans(rectangleEnvXPos, rectangleEnvYPos, rectangleEnvZPos);
          dd4hep::PlacedVolume rectangleEnvelopePhys;
          if (side % 2 == 0) {
            rectangleEnvelopePhys =
                sideVol.placeVolume(rectangleEnvVol, dd4hep::Transform3D(rotationY, rectangeEnvelopeTrans));
          } else {
            rectangleEnvelopePhys =
                sideVol2.placeVolume(rectangleEnvVol, dd4hep::Transform3D(rotationY, rectangeEnvelopeTrans));
          }
          dd4hep::DetElement rectangleEnvelopeDE(sideDE, rectangleEnvelopeName + "DE", rectangle);
          rectangleEnvelopeDE.setPlacement(rectangleEnvelopePhys);
          rectangleEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

          //  ------------------------ start to build the chamber envelopes -------------------

          int numChambersInRectangle;
          if (rectangleEnvY <= dimensions.y()) {
            numChambersInRectangle =
                0; // number of the chambers in each rectangle, in that case it will create just 1 chamber.
          } else {
            numChambersInRectangle =
                2 * rectangleEnvY / (2 * dimensions.y() - overlapY); // number of the chambers in each rectangle
          }

          for (int chamberIndex = 0; chamberIndex < (numChambersInRectangle + 1); chamberIndex++) {

            std::stringstream barrelNameStream;
            barrelNameStream << "-MuRWELL_Barrel_" << barrelIdCounter++;
            std::string BarrelChamberName = name + barrelNameStream.str();

            dd4hep::Box envelope(dimensions.x(), dimensions.y(), dimensions.z());
            dd4hep::Volume envVolume(BarrelChamberName, envelope, lcdd.material(dimensions.materialStr()));

            double rectangleRemainderY;
            double rectangleRemainderREnvYPos;
            if (numChambersInRectangle == 0) {
              rectangleRemainderY =
                  std::abs(std::fmod((2 * rectangleEnvY - clearance), (2 * dimensions.y() - clearance))) /
                  (2 * dimensions.y());
              rectangleRemainderREnvYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) +
                                           dimensions.y_offset() - rectangleEnvY +
                                           rectangleRemainderY * dimensions.y() + clearance / 20;
            } else {
              rectangleRemainderY =
                  std::fmod(2 * (rectangleEnvY - clearance), (2 * dimensions.y() - overlapY)) / (2 * dimensions.y());
              rectangleRemainderREnvYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) +
                                           dimensions.y_offset() - rectangleEnvY +
                                           rectangleRemainderY * dimensions.y() + 1.5 * clearance;
            }

            dd4hep::Box rectangleRemainderYEnvelope(dimensions.x(), rectangleRemainderY * dimensions.y(),
                                                    dimensions.z());
            dd4hep::Volume rectangleRemainderYEnvVolume(BarrelChamberName + "rectangleRemainderY",
                                                        rectangleRemainderYEnvelope,
                                                        lcdd.material(dimensions.materialStr()));

            double envYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() -
                             rectangleEnvY + dimensions.y() +
                             clearance / 20.0; // found that the positioning of the chambers inside the rectangle had an
                                               // overlap with the mother volume ~ 45 um.
            // double rectangleRemainderREnvYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) +
            // dimensions.y_offset() - rectangleEnvY + rectangleRemainderY * dimensions.y() + 1.5 * clearance;

            double zRotation;
            double rectangleRemainderZRotation;
            if (numChambersInRectangle == 0) {
              zRotation = 0.0;
              rectangleRemainderZRotation = 0.0;
            } else {
              zRotation = std::atan(dimensions.x() / (dimensions.y() - (2 * overlapY)));
              rectangleRemainderZRotation =
                  std::atan(dimensions.x() / (rectangleRemainderY * dimensions.y() -
                                              overlapY * 1.5)); // Y and Z are reversed in local remainder
            }

            dd4hep::RotationZ chamberRotation(zRotation);
            dd4hep::RotationZ rectangleRemainderRotationZ(rectangleRemainderZRotation);

            auto Slices = xmlElement.children(_Unicode(slice));
            auto numSlices = xmlElement.numChildren(_Unicode(slice), true);
            dd4hep::xml::Handle_t slice(Slices.reset());
            int sensitiveSliceIndex = 0;

            // --- two cases : one for remainder chambers ----------------------------------------------

            if (chamberIndex == numChambersInRectangle) {

              dd4hep::Position rectangleRemainderTrans(dimensions.x_offset(), rectangleRemainderREnvYPos, 0.0);
              dd4hep::PlacedVolume rectangleRemainderEnvPhys = rectangleEnvVol.placeVolume(
                  rectangleRemainderYEnvVolume,
                  dd4hep::Transform3D(rectangleRemainderRotationZ, rectangleRemainderTrans));
              rectangleRemainderEnvPhys.addPhysVolID("chamber", barrelIdCounter);
              dd4hep::DetElement rectangleRemainderEnvDE(rectangleEnvelopeDE, BarrelChamberName, barrelIdCounter);
              rectangleRemainderEnvDE.setPlacement(rectangleRemainderEnvPhys);
              rectangleRemainderYEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());

              // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() <<
              // std::endl;
              double sliceXOffset = -dimensions.x();
              for (unsigned sliceIdx = 0; sliceIdx < numSlices; ++sliceIdx) {
                dd4hep::xml::DetElement sliceDet = static_cast<dd4hep::xml::DetElement>(slice);
                dd4hep::Box sliceShape(sliceDet.x(), rectangleRemainderY * dimensions.y(), dimensions.z());
                std::string sliceName = dd4hep::xml::_toString(sliceIdx, "slice%d");
                dd4hep::Volume sliceVolume(sliceName, sliceShape, lcdd.material(slice.attr<std::string>("material")));
                dd4hep::Position transSlice(sliceXOffset + sliceDet.x(), 0.0, 0.0);
                dd4hep::PlacedVolume slicePlacedVolume = rectangleRemainderYEnvVolume.placeVolume(
                    sliceVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), transSlice));

                if (slice.hasAttr("vis")) {
                  sliceVolume.setVisAttributes(lcdd, sliceDet.visStr());
                }
                if (slice.hasAttr("sensitive") && sliceDet.isSensitive()) {
                  dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                  sensDet.setType(sdType.typeStr());
                  sliceVolume.setSensitiveDetector(sensDet);
                  slicePlacedVolume.addPhysVolID("slice", sensitiveSliceIndex);
                  dd4hep::DetElement sliceDE(rectangleRemainderEnvDE, "slice_" + std::to_string(sensitiveSliceIndex),
                                             sensitiveSliceIndex);
                  sliceDE.setPlacement(slicePlacedVolume);
                  // dd4hep::printout(dd4hep::INFO,"Sensitive slice has been created at", name, BarrelChamberName);
                  sensitiveSliceIndex++;
                }
                // Increment the current x-offset by the width of the current slice
                sliceXOffset += (2 * sliceDet.x());
                slice.m_node = Slices.next();
              }

              // ---------------- Second case: for the full chambers
            } else {

              dd4hep::Position trans(dimensions.x_offset(), envYPos, 0.0);
              dd4hep::PlacedVolume envPhys =
                  rectangleEnvVol.placeVolume(envVolume, dd4hep::Transform3D(chamberRotation, trans));
              envPhys.addPhysVolID("chamber", barrelIdCounter);
              dd4hep::DetElement envDE(rectangleEnvelopeDE, BarrelChamberName, barrelIdCounter);
              envDE.setPlacement(envPhys);
              envVolume.setVisAttributes(lcdd, xmlDet.visStr());

              // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() <<
              // std::endl;
              double sliceXOffset = -dimensions.x();
              for (unsigned sliceIdx = 0; sliceIdx < numSlices; ++sliceIdx) {
                dd4hep::xml::DetElement sliceDet = static_cast<dd4hep::xml::DetElement>(slice);
                dd4hep::Box sliceShape(sliceDet.x(), dimensions.y(), dimensions.z());
                std::string sliceName = dd4hep::xml::_toString(sliceIdx, "slice%d");
                dd4hep::Volume sliceVolume(sliceName, sliceShape, lcdd.material(slice.attr<std::string>("material")));
                dd4hep::Position transSlice(sliceXOffset + sliceDet.x(), 0.0, 0.0);
                dd4hep::PlacedVolume slicePlacedVolume =
                    envVolume.placeVolume(sliceVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), transSlice));

                if (slice.hasAttr("vis")) {
                  sliceVolume.setVisAttributes(lcdd, sliceDet.visStr());
                }
                if (slice.hasAttr("sensitive") && sliceDet.isSensitive()) {
                  dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                  sensDet.setType(sdType.typeStr());
                  sliceVolume.setSensitiveDetector(sensDet);
                  slicePlacedVolume.addPhysVolID("slice", sensitiveSliceIndex);
                  dd4hep::DetElement sliceDE(envDE, "slice_" + std::to_string(sensitiveSliceIndex),
                                             sensitiveSliceIndex);
                  sliceDE.setPlacement(slicePlacedVolume);
                  // dd4hep::printout(dd4hep::INFO,"Sensitive slice has been created at", name, BarrelChamberName);
                  sensitiveSliceIndex++;
                }
                // Increment the current x-offset by the width of the current slice
                sliceXOffset += (2 * sliceDet.x());
                slice.m_node = Slices.next();
              }
            }
          }
        }
      }
      //-----------------------------------------------------------------------------------
      //-----------------------------------------------------------------------------------
    }
  }

  // ---------------------------------------// B A R R E L // R A D I A T O R S
  // //------------------------------------------------

  for (int numBarrelRadiatorLayer = 0; numBarrelRadiatorLayer < numBarrelRadiators; ++numBarrelRadiatorLayer) {

    double barrelRadiatorLayerRMin =
        radius + (numBarrelRadiatorLayer + 1) * (2 * detectorVolumeThickness) +
        numBarrelRadiatorLayer * barrelRadiatorThickness; // I'm assuming that the yoke is following the detctor layer;
                                                          // one detector layer then -> one yoke layer.
    double barrelRadiatorLayerRMax = barrelRadiatorLayerRMin + barrelRadiatorThickness;
    double barrelRadiatorLayerRMid = (barrelRadiatorLayerRMin + barrelRadiatorLayerRMax) / 2.0;

    dd4hep::PolyhedraRegular BarrelRadiatorLayer(numSides, barrelRadiatorLayerRMin, barrelRadiatorLayerRMax,
                                                 barrelLength);
    std::string barrelRadiatorEnvName = name + "-BarrelRadiatorLayer" + std::to_string(numBarrelRadiatorLayer + 1);
    dd4hep::Volume BarrelRadiatorLayerVolume(barrelRadiatorEnvName, BarrelRadiatorLayer, mat);

    double barrelRadiatorSideLength =
        2 * (barrelRadiatorLayerRMid - barrelRadiatorThickness / 2.0) * std::tan(shapeAngle_radians);
    double barrelRadiatorSideLength2 =
        2 * (barrelRadiatorLayerRMid + barrelRadiatorThickness / 2.0) * std::tan(shapeAngle_radians);

    for (int side = 0; side < numSides; ++side) {
      int sideID =
          (numBarrelRadiatorLayer + 1) * 100 + (side + 1); // to differentiated with the same side in different layers.
      dd4hep::Trapezoid barrelRadiatorEnvelope(barrelLength / 2.0, barrelLength / 2.0, barrelRadiatorSideLength / 2.0,
                                               barrelRadiatorSideLength2 / 2.0, barrelRadiatorThickness / 2.0);
      std::string barrelRadiatorEnvelopeName = name + "-BarrelRadiatorSide" + std::to_string(sideID);
      dd4hep::Volume barrelRadiatorEnvVol(barrelRadiatorEnvelopeName, barrelRadiatorEnvelope, mat);

      double angle_degrees = shapeAngle * side; // Calculate the angle for each side
      double angle_radians = angle_degrees * M_PI / 180.0;

      double barrelRadiatorXOffset = barrelRadiatorLayerRMid * std::cos(angle_radians + shapeAngle_radians);
      double barrelRadiatorYOffset = barrelRadiatorLayerRMid * std::sin(angle_radians + shapeAngle_radians);

      dd4hep::RotationZ barrelRadiatorRotationZ(angle_radians + shapeAngle_radians);
      dd4hep::RotationY barrelRadiatorRotationY(90. * dd4hep::degree);

      dd4hep::Rotation3D barrelRadiatorRotation = barrelRadiatorRotationZ * barrelRadiatorRotationY;

      double barrelRadiatorXPos = barrelRadiatorXOffset;
      double barrelRadiatorYPos = barrelRadiatorYOffset;
      double barrelRadiatorZPos = 0.0;

      dd4hep::Position barrelRadiatorEnvelopeTrans(barrelRadiatorXPos, barrelRadiatorYPos, barrelRadiatorZPos);
      dd4hep::PlacedVolume barrelRadiatorEnvelopePhys = BarrelRadiatorLayerVolume.placeVolume(
          barrelRadiatorEnvVol, dd4hep::Transform3D(barrelRadiatorRotation, barrelRadiatorEnvelopeTrans));
      dd4hep::DetElement barrelRadiatorEnvelopeDE(detElement, barrelRadiatorEnvelopeName + "DE", sideID);
      barrelRadiatorEnvelopeDE.setPlacement(barrelRadiatorEnvelopePhys);
      barrelRadiatorEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

      // ----------------------------------------- Building the yokes ------------------

      dd4hep::Trapezoid barrelRadiator(barrelLength / 2.0, barrelLength / 2.0, barrelRadiatorSideLength / 2.0,
                                       barrelRadiatorSideLength2 / 2.0, barrelRadiatorThickness / 2.0);
      std::string barrelRadiatorName = dd4hep::xml::_toString(side, "barrelRadiator%d");
      dd4hep::Volume barrelRadiatorVol(barrelRadiatorName, barrelRadiator, barrelRadiatorMaterial);

      dd4hep::Position barrelRadiatorTrans(0, 0, 0);
      dd4hep::PlacedVolume barrelRadiatorPhys = barrelRadiatorEnvVol.placeVolume(
          barrelRadiatorVol, dd4hep::Transform3D(dd4hep::RotationY(0.), barrelRadiatorTrans));
      dd4hep::DetElement barrelRadiatorDE(barrelRadiatorEnvelopeDE, barrelRadiatorName + "DE", sideID);
      barrelRadiatorDE.setPlacement(barrelRadiatorPhys);
      barrelRadiatorVol.setVisAttributes(lcdd.visAttributes("yoke_vis"));
    }

    dd4hep::Position radiatorLayerTrans(0., 0., 0.);
    dd4hep::PlacedVolume radiatorLayerPhys = BarrelVolume.placeVolume(
        BarrelRadiatorLayerVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), radiatorLayerTrans));
    dd4hep::DetElement BarrelRadiatorLayerDE(
        BarrelDE, name + "-Barrel_RadiatorLayerDE_" + std::to_string(numBarrelRadiatorLayer + 1),
        numBarrelRadiatorLayer + 1);
    BarrelRadiatorLayerDE.setPlacement(radiatorLayerPhys);
    BarrelRadiatorLayerVolume.setVisAttributes(lcdd.visAttributes("no_vis"));
  }

  // -------------------------------------------------------------------------------------------
  //  ----------------------------------// E N D C A P //---------------------------------------
  //--------------------------------------- Endcap Detectors------------------------------------

  std::string EndcapName;
  dd4hep::Volume endcapVolume;
  dd4hep::Position endcapTrans;
  dd4hep::PlacedVolume endcapPhys;
  dd4hep::DetElement EndcapDE;

  for (int endcapType = -1; endcapType < 2; ++endcapType) {
    if (endcapType == 0) {
      continue; // Skip the iteration when endcapType is 0
    }
    EndcapName = name + "-Endcap" + std::to_string(endcapType);
    endcapVolume = dd4hep::Volume(EndcapName, EndcapEnv, mat);
    endcapTrans = dd4hep::Position(0., 0., endcapType * endcapOffset);
    endcapPhys = detectorVolume.placeVolume(endcapVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), endcapTrans));
    endcapPhys.addPhysVolID("type", endcapType);
    EndcapDE = dd4hep::DetElement(detElement, name + "EndcapDE_" + std::to_string(endcapType), endcapType);
    EndcapDE.setPlacement(endcapPhys);
    endcapVolume.setVisAttributes(lcdd.visAttributes("no_vis"));

    // ------------------ loop to create endcap layers ------------------

    for (int numEndcapLayer = 0; numEndcapLayer < numEndcapDetectorLayers; ++numEndcapLayer) {

      double endcapLayerZOffset;
      std::string endcapDetectorEnvelopeName;
      double endcapDetectorEnvZPos;

      dd4hep::PolyhedraRegular endcapDetectorEnvelope(numSides, endcapDetectorLayerInnerRadius,
                                                      endcapDetectorLayerOuterRadius, 2 * endcapDetectorEnvZ);

      // Automation of inner Z-Offset of different layers, taking into account that
      // every detector layer is followed by a yoke(radiator) layer
      if (endcapType == -1) {
        // For negative endcap: reverse the layer positioning so layer 0 is closest to IP
        int reversedLayer = numEndcapDetectorLayers - 1 - numEndcapLayer;
        endcapLayerZOffset = -EndcaptotalLength / 2.0 + endcapDetectorEnvZ + reversedLayer * (2 * endcapDetectorEnvZ) +
                             reversedLayer * endcapRadiatorThickness;
      } else {
        // For positive endcap: use normal positioning
        endcapLayerZOffset = -EndcaptotalLength / 2.0 + endcapDetectorEnvZ + numEndcapLayer * (2 * endcapDetectorEnvZ) +
                             numEndcapLayer * endcapRadiatorThickness;
      }

      endcapDetectorEnvelopeName = name + "-EndcapDetectorLayer" + std::to_string(numEndcapLayer + 1);
      endcapDetectorEnvZPos = endcapLayerZOffset;

      dd4hep::Volume endcapDetectorEnvVol(endcapDetectorEnvelopeName, endcapDetectorEnvelope, mat);

      double endcapDetectorEnvXPos = 0.0;
      double endcapDetectorEnvYPos = 0.0;

      int detElementID = (numEndcapLayer < 0) ? numEndcapLayer + numEndcapDetectorLayers
                                              : numEndcapLayer + numEndcapDetectorLayers + 1;

      dd4hep::Position endcapDetectorEnvelopeTrans(endcapDetectorEnvXPos, endcapDetectorEnvYPos, endcapDetectorEnvZPos);
      dd4hep::PlacedVolume endcapDetectorEnvelopePhys = endcapVolume.placeVolume(
          endcapDetectorEnvVol, dd4hep::Transform3D(dd4hep::RotationZ(0.), endcapDetectorEnvelopeTrans));
      dd4hep::DetElement endcapDetectorEnvelopeDE(EndcapDE, endcapDetectorEnvelopeName + "DE", detElementID);
      endcapDetectorEnvelopeDE.setPlacement(endcapDetectorEnvelopePhys);
      endcapDetectorEnvelopePhys.addPhysVolID("layer", numEndcapLayer);
      endcapDetectorEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

      // -------------------------- Building detector sides ---------------------

      for (int side = 0; side < numSides; ++side) {

        int sideID =
            (numEndcapLayer + 1) * 100 + (side + 1); // to differentiated with the same side in different layers.
        dd4hep::Trapezoid endcapDetectorSideTrap(detectorVolumeThickness / 2.0, detectorVolumeThickness / 2.0,
                                                 endcapDetectorSideLength / 2.0, endcapDetectorSideTrapLength / 2.0,
                                                 endcapDetectorSideTrapYLength / 2.0);
        dd4hep::Box endcapDetectorSideBox(detectorVolumeThickness / 2.0, endcapDetectorSideBoxLength / 2.0,
                                          endcapDetectorSideBoxYLength / 2.0);

        double boxOffsetZ = endcapDetectorYLength / 2.0;

        dd4hep::Position boxPos(0.0, 0.0, boxOffsetZ);
        dd4hep::Rotation3D boxRot(dd4hep::RotationY(0.0 * dd4hep::degree));

        dd4hep::Transform3D boxTransform(boxRot, boxPos);

        // Combining two shapes by UnionSolid: the first shape is centralized and the second transform around the
        // first..
        dd4hep::UnionSolid endcapDetectorSideEnvelope(endcapDetectorSideTrap, endcapDetectorSideBox, boxTransform);
        std::string endcapDetectorSideEnvName = dd4hep::xml::_toString(sideID, "endcapDetectorSideEnv%d");
        dd4hep::Volume endcapDetectorSideEnvVol(endcapDetectorSideEnvName, endcapDetectorSideEnvelope, mat);

        double angle_degrees = shapeAngle * side; // Calculate the angle for each side
        double angle_radians = angle_degrees * M_PI / 180.0;

        // double endcapDetectorMidRadius = endcapDetectorLayerInnerRadius + (endcapDetectorYLength /2.0);
        double endcapDetectorTrapCenterRadius = endcapDetectorLayerInnerRadius + (endcapDetectorSideTrapYLength / 2.0);

        double endcapDetectorXOffset = endcapDetectorTrapCenterRadius * std::cos(angle_radians + shapeAngle_radians);
        double endcapDetectorYOffset = endcapDetectorTrapCenterRadius * std::sin(angle_radians + shapeAngle_radians);

        dd4hep::RotationZ endcapDetectorRotationZ(angle_radians + shapeAngle_radians);
        dd4hep::Rotation3D endcapDetectorRotation = dd4hep::Rotation3D(endcapDetectorRotationZ);

        double endcapDetectorSideEnvXPos = endcapDetectorXOffset;
        double endcapDetectorSideEnvYPos = endcapDetectorYOffset;
        double endcapDetectorSideEnvZPos;
        double endcapDetectorSideEnvZPos2;

        dd4hep::Position endcapDetectorSideEnvTrans;
        dd4hep::PlacedVolume endcapDetectorSideEnvPhys;
        dd4hep::DetElement endcapDetectorSideEnvDE;

        if (sideLengthEstimate <=
            (2 * dimensions.y())) { // in case of small side's length, put the sides in 4 differnet positions, to avoid
                                    // overlaps between different sides, and the probability increases for small sides
          endcapDetectorSideEnvZPos = -endcapDetectorEnvZ / 4.0;
          endcapDetectorSideEnvZPos2 = -endcapDetectorEnvZ * 3.0 / 4.0;
          double endcapDetectorSideEnvZPos3 = endcapDetectorEnvZ / 4.0;
          double endcapDetectorSideEnvZPos4 = endcapDetectorEnvZ * 3.0 / 4.0;

          if (side % 4 == 0) {
            endcapDetectorSideEnvTrans =
                dd4hep::Position(endcapDetectorSideEnvXPos, endcapDetectorSideEnvYPos, endcapDetectorSideEnvZPos);
          } else if (side % 4 == 1) {
            endcapDetectorSideEnvTrans =
                dd4hep::Position(endcapDetectorSideEnvXPos, endcapDetectorSideEnvYPos, endcapDetectorSideEnvZPos2);
          } else if (side % 4 == 2) {
            endcapDetectorSideEnvTrans =
                dd4hep::Position(endcapDetectorSideEnvXPos, endcapDetectorSideEnvYPos, endcapDetectorSideEnvZPos3);
          } else {
            endcapDetectorSideEnvTrans =
                dd4hep::Position(endcapDetectorSideEnvXPos, endcapDetectorSideEnvYPos, endcapDetectorSideEnvZPos4);
          }

        } else { // else, put the sides in 2 differnet positions.
          endcapDetectorSideEnvZPos = -endcapDetectorEnvZ / 2.0;
          endcapDetectorSideEnvZPos2 = endcapDetectorEnvZ / 2.0;

          if (side % 2 == 0) {
            endcapDetectorSideEnvTrans =
                dd4hep::Position(endcapDetectorSideEnvXPos, endcapDetectorSideEnvYPos, endcapDetectorSideEnvZPos);
          } else {
            endcapDetectorSideEnvTrans =
                dd4hep::Position(endcapDetectorSideEnvXPos, endcapDetectorSideEnvYPos, endcapDetectorSideEnvZPos2);
          }
        }

        endcapDetectorSideEnvPhys = endcapDetectorEnvVol.placeVolume(
            endcapDetectorSideEnvVol,
            dd4hep::Transform3D(endcapDetectorRotation * dd4hep::RotationY(90.0 * dd4hep::degree),
                                endcapDetectorSideEnvTrans));
        endcapDetectorSideEnvDE =
            dd4hep::DetElement(endcapDetectorEnvelopeDE, endcapDetectorSideEnvName + "DE", sideID);
        endcapDetectorSideEnvDE.setPlacement(endcapDetectorSideEnvPhys);
        endcapDetectorSideEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

        // ----- dividing the trapezoid envelope to smaller pieces (rectangles)

        int numRectangles =
            endcapDetectorYLength /
            (2 * dimensions.z() - overlapZ); // numbers of the rectangles in each teapezoids.. depends on the number of
                                             // the chambers that can be overlapped in Y direction

        for (int rectangle = 0; rectangle < (numRectangles + 1); ++rectangle) {

          double rectangleEnvY =
              (endcapDetectorLayerOuterRadius - rectangle * (2 * dimensions.y() - overlapY)) *
              std::tan(shapeAngle_radians); // without multiplying by 2 .. because it is the half length // it should be
                                            // dimensions.x() instead of z, but in the endcap its perpendicular to usual
                                            // dimension set
          double rectangleEnvZ;
          if (rectangle == numRectangles) {
            rectangleEnvZ = (endcapRemainderZ * dimensions.z()) + clearance / 4.0;
          } else {
            rectangleEnvZ = dimensions.z() + clearance / 2.0;
          }
          double rectangleEnvX;
          if (rectangleEnvY <= dimensions.y()) {
            rectangleEnvX = dimensions.x(); // in case if there is only one chamber inside the rectangle.
          } else {
            rectangleEnvX = detectorVolumeThickness /
                            4.5; // dividing by 4.5 gets the best thickness for the recatngle to avoid any overlap ~ in
                                 // our case the uRWELL the best rectangle thickness is 26.667 mm which is 120/4.5.
          }

          dd4hep::Box rectangleEnvelope(rectangleEnvX, rectangleEnvY, rectangleEnvZ);
          std::string rectangleEnvelopeName;
          if (rectangle == numRectangles) {
            rectangleEnvelopeName = dd4hep::xml::_toString(rectangle, "rectangleRemainderEnvelope%d");
          } else {
            rectangleEnvelopeName = dd4hep::xml::_toString(rectangle, "rectangleEnvelope%d");
          }
          dd4hep::Volume rectangleEnvVol(rectangleEnvelopeName, rectangleEnvelope, mat);

          double rectangleEnvXPos = 0.0;
          double rectangleEnvYPos = 0.0;
          double rectangleEnvZPos;
          if (rectangle == numRectangles) {
            rectangleEnvZPos = endcapDetectorYLength / 2.0 + ((1 - endcapRemainderZ) * dimensions.z()) -
                               (rectangle * (2 * dimensions.z() - overlapZ)) - clearance;
          } else {
            rectangleEnvZPos = endcapDetectorYLength / 2.0 - (rectangle * (2 * dimensions.y() - overlapY)) - clearance;
          }

          double yRotation = std::atan(rectangleEnvX / (rectangleEnvZ - (2 * overlapY)));
          dd4hep::RotationY rotationY(yRotation);

          dd4hep::Position rectangeEnvelopeTrans(rectangleEnvXPos, rectangleEnvYPos, rectangleEnvZPos);
          dd4hep::PlacedVolume rectangleEnvelopePhys = endcapDetectorSideEnvVol.placeVolume(
              rectangleEnvVol, dd4hep::Transform3D(rotationY, rectangeEnvelopeTrans));
          dd4hep::DetElement rectangleEnvelopeDE(endcapDetectorSideEnvDE, rectangleEnvelopeName + "DE", rectangle);
          rectangleEnvelopeDE.setPlacement(rectangleEnvelopePhys);
          rectangleEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

          // ------------------------ start to build the chamber envelopes -------------------

          int numChambersInRectangle;
          if (rectangleEnvY <= dimensions.y()) {
            numChambersInRectangle = 0; // number of the chambers in each rectangle
          } else {
            numChambersInRectangle =
                2 * rectangleEnvY / (2 * dimensions.y() - overlapY); // number of the chambers in each rectangle
          }

          for (int chamberIndex = 0; chamberIndex < (numChambersInRectangle + 1); chamberIndex++) {

            std::stringstream endcapNameStream;
            endcapNameStream << "-MuRWELL_Endcap_" << endcapIdCounter++;
            std::string EndcapChamberName = name + endcapNameStream.str();

            dd4hep::Box envelope;
            if (rectangle == numRectangles) {
              envelope = dd4hep::Box(dimensions.x(), dimensions.y(), endcapRemainderZ * dimensions.z());
            } else {
              envelope = dd4hep::Box(dimensions.x(), dimensions.y(), dimensions.z());
            }

            dd4hep::Volume envVolume(EndcapChamberName, envelope, lcdd.material(dimensions.materialStr()));

            double rectangleRemainderY;
            if (numChambersInRectangle == 0) {
              rectangleRemainderY =
                  std::abs(std::fmod((2 * rectangleEnvY - clearance), (2 * dimensions.y() - clearance))) /
                  (2 * dimensions.y());
            } else {
              rectangleRemainderY =
                  std::fmod(2 * (rectangleEnvY - clearance), (2 * dimensions.y() - overlapY)) / (2 * dimensions.y());
            }

            dd4hep::Box rectangleRemainderYEnvelope;
            if (rectangle == numRectangles) {
              rectangleRemainderYEnvelope =
                  dd4hep::Box(dimensions.x(), rectangleRemainderY * dimensions.y(), endcapRemainderZ * dimensions.z());
            } else {
              rectangleRemainderYEnvelope =
                  dd4hep::Box(dimensions.x(), rectangleRemainderY * dimensions.y(), dimensions.z());
            }

            dd4hep::Volume rectangleRemainderYEnvVolume(EndcapChamberName + "rectangleRemainderY",
                                                        rectangleRemainderYEnvelope,
                                                        lcdd.material(dimensions.materialStr()));

            double envYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() -
                             rectangleEnvY + dimensions.y() +
                             0.005; // found that the positioning of the chambers inside the rectangle had an overlap
                                    // with the mother volume ~ 45 um.
            double rectangleRemainderREnvYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) +
                                                dimensions.y_offset() - rectangleEnvY +
                                                rectangleRemainderY * dimensions.y() + 0.005;

            double zRotation;
            double rectangleRemainderZRotation;
            if (numChambersInRectangle == 0) {
              zRotation = 0.0;
              rectangleRemainderZRotation = 0.0;
            } else {
              zRotation = std::atan(dimensions.x() / (dimensions.y() - (2 * overlapY)));
              rectangleRemainderZRotation =
                  std::atan(dimensions.x() / (rectangleRemainderY * dimensions.y() -
                                              (2 * overlapY))); // Y and Z are reversed in local remainder
            }

            dd4hep::RotationZ chamberRotation(zRotation);
            dd4hep::RotationZ rectangleRemainderRotationZ(rectangleRemainderZRotation);

            auto Slices = xmlElement.children(_Unicode(slice));
            auto numSlices = xmlElement.numChildren(_Unicode(slice), true);
            dd4hep::xml::Handle_t slice(Slices.reset());
            int sensitiveSliceIndex = 0;

            // --- two cases : one for full chambers ----------------------------------------------

            if (chamberIndex == numChambersInRectangle) {

              dd4hep::Position rectangleRemainderTrans(dimensions.x_offset(), rectangleRemainderREnvYPos, 0.0);
              dd4hep::PlacedVolume rectangleRemainderEnvPhys = rectangleEnvVol.placeVolume(
                  rectangleRemainderYEnvVolume,
                  dd4hep::Transform3D(rectangleRemainderRotationZ, rectangleRemainderTrans));
              rectangleRemainderEnvPhys.addPhysVolID("chamber", endcapIdCounter);
              dd4hep::DetElement rectangleRemainderEnvDE(rectangleEnvelopeDE, EndcapChamberName, endcapIdCounter);
              rectangleRemainderEnvDE.setPlacement(rectangleRemainderEnvPhys);
              rectangleRemainderYEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());

              // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() <<
              // std::endl;
              double sliceXOffset = -dimensions.x();
              for (unsigned sliceIdx = 0; sliceIdx < numSlices; ++sliceIdx) {
                dd4hep::xml::DetElement sliceDet = static_cast<dd4hep::xml::DetElement>(slice);
                dd4hep::Box sliceShape;
                if (rectangle == numRectangles) {
                  sliceShape = dd4hep::Box(sliceDet.x(), rectangleRemainderY * dimensions.y(),
                                           endcapRemainderZ * dimensions.z());
                } else {
                  sliceShape = dd4hep::Box(sliceDet.x(), rectangleRemainderY * dimensions.y(), dimensions.z());
                }
                std::string sliceName = dd4hep::xml::_toString(sliceIdx, "slice%d");
                dd4hep::Volume sliceVolume(sliceName, sliceShape, lcdd.material(slice.attr<std::string>("material")));
                dd4hep::Position transSlice(sliceXOffset + sliceDet.x(), 0.0, 0.0);
                dd4hep::PlacedVolume slicePlacedVolume = rectangleRemainderYEnvVolume.placeVolume(
                    sliceVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), transSlice));

                if (slice.hasAttr("vis")) {
                  sliceVolume.setVisAttributes(lcdd, sliceDet.visStr());
                }
                if (slice.hasAttr("sensitive") && sliceDet.isSensitive()) {
                  dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                  sensDet.setType(sdType.typeStr());
                  sliceVolume.setSensitiveDetector(sensDet);
                  slicePlacedVolume.addPhysVolID("slice", sensitiveSliceIndex);
                  dd4hep::DetElement sliceDE(rectangleRemainderEnvDE, "slice_" + std::to_string(sensitiveSliceIndex),
                                             sensitiveSliceIndex);
                  sliceDE.setPlacement(slicePlacedVolume);
                  //  dd4hep::printout(dd4hep::INFO,"Sensitive slice has been created at", name,EndcapChamberName);
                  sensitiveSliceIndex++;
                }
                // Increment the current x-offset by the width of the current slice
                sliceXOffset += (2 * sliceDet.x());
                slice.m_node = Slices.next();
              }

              // ----------------
            } else {

              dd4hep::Position trans(dimensions.x_offset(), envYPos, 0.0);
              dd4hep::PlacedVolume envPhys =
                  rectangleEnvVol.placeVolume(envVolume, dd4hep::Transform3D(chamberRotation, trans));
              envPhys.addPhysVolID("chamber", endcapIdCounter);
              dd4hep::DetElement envDE(rectangleEnvelopeDE, EndcapChamberName, endcapIdCounter);
              envDE.setPlacement(envPhys);
              envVolume.setVisAttributes(lcdd, xmlDet.visStr());

              // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() <<
              // std::endl;
              double sliceXOffset = -dimensions.x();
              for (unsigned sliceIdx = 0; sliceIdx < numSlices; ++sliceIdx) {
                dd4hep::xml::DetElement sliceDet = static_cast<dd4hep::xml::DetElement>(slice);
                dd4hep::Box sliceShape;
                if (rectangle == numRectangles) {
                  sliceShape = dd4hep::Box(sliceDet.x(), dimensions.y(), endcapRemainderZ * dimensions.z());
                } else {
                  sliceShape = dd4hep::Box(sliceDet.x(), dimensions.y(), dimensions.z());
                }
                std::string sliceName = dd4hep::xml::_toString(sliceIdx, "slice%d");
                dd4hep::Volume sliceVolume(sliceName, sliceShape, lcdd.material(slice.attr<std::string>("material")));
                dd4hep::Position transSlice(sliceXOffset + sliceDet.x(), 0.0, 0.0);
                dd4hep::PlacedVolume slicePlacedVolume =
                    envVolume.placeVolume(sliceVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), transSlice));

                if (slice.hasAttr("vis")) {
                  sliceVolume.setVisAttributes(lcdd, sliceDet.visStr());
                }
                if (slice.hasAttr("sensitive") && sliceDet.isSensitive()) {
                  dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                  sensDet.setType(sdType.typeStr());
                  sliceVolume.setSensitiveDetector(sensDet);
                  slicePlacedVolume.addPhysVolID("slice", sensitiveSliceIndex);
                  dd4hep::DetElement sliceDE(envDE, "slice_" + std::to_string(sensitiveSliceIndex),
                                             sensitiveSliceIndex);
                  sliceDE.setPlacement(slicePlacedVolume);
                  // dd4hep::printout(dd4hep::INFO,"Sensitive slice has been created at", name, EndcapChamberName);
                  sensitiveSliceIndex++;
                }
                // Increment the current x-offset by the width of the current slice
                sliceXOffset += (2 * sliceDet.x());
                slice.m_node = Slices.next();
              }
            }
          }
        }
      }
    }
    // --------------------------------------Radiators--------------------------------------------

    for (int numEndcapRadiatorLayer = 0; numEndcapRadiatorLayer < numEndcapRadiators; ++numEndcapRadiatorLayer) {

      // Automation of inner Z-Offset of different layers, taking into account that
      // every detector layer is followed by a yoke(radiator) layer
      double endcapRadiatorLayerZOffset;
      if (endcapType == -1) {
        // For negative endcap: reverse the radiator layer positioning
        int reversedRadiatorLayer = numEndcapRadiators - 1 - numEndcapRadiatorLayer;
        endcapRadiatorLayerZOffset = -EndcaptotalLength / 2.0 + (endcapRadiatorThickness / 2.0) +
                                     (reversedRadiatorLayer + 1) * (2 * endcapDetectorEnvZ) +
                                     reversedRadiatorLayer * endcapRadiatorThickness;
      } else {
        // For positive endcap: use normal positioning
        endcapRadiatorLayerZOffset = -EndcaptotalLength / 2.0 + (endcapRadiatorThickness / 2.0) +
                                     (numEndcapRadiatorLayer + 1) * (2 * endcapDetectorEnvZ) +
                                     numEndcapRadiatorLayer * endcapRadiatorThickness;
      }

      dd4hep::PolyhedraRegular endcapRadiatorEnvelope(numSides, endcapRadiatorLayerInnerRadius,
                                                      endcapRadiatorLayerOuterRadius, endcapRadiatorThickness);
      std::string endcapRadiatorEnvelopeName =
          name + "-EndcapRadiatorLayer" + std::to_string(numEndcapRadiatorLayer + 1);
      dd4hep::Volume endcapRadiatorEnvVol(endcapRadiatorEnvelopeName, endcapRadiatorEnvelope, mat);

      double endcapRadiatorEnvXPos = 0.0;
      double endcapRadiatorEnvYPos = 0.0;
      double endcapRadiatorEnvZPos = endcapRadiatorLayerZOffset;

      dd4hep::Position endcapRadiatorEnvelopeTrans(endcapRadiatorEnvXPos, endcapRadiatorEnvYPos, endcapRadiatorEnvZPos);
      dd4hep::PlacedVolume endcapRadiatorEnvelopePhys = endcapVolume.placeVolume(
          endcapRadiatorEnvVol, dd4hep::Transform3D(dd4hep::RotationZ(0.), endcapRadiatorEnvelopeTrans));
      dd4hep::DetElement endcapRadiatorEnvelopeDE(EndcapDE, endcapRadiatorEnvelopeName + "DE",
                                                  numEndcapRadiatorLayer + 1);
      endcapRadiatorEnvelopeDE.setPlacement(endcapRadiatorEnvelopePhys);
      endcapRadiatorEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

      // ----------------Building radiator sides------------------

      for (int side = 0; side < numSides; ++side) {

        int sideID = (numEndcapRadiatorLayer + 1) * 100 +
                     (side + 1); // to differentiated with the same side in different layers.
        dd4hep::Trapezoid endcapRadiator(endcapRadiatorThickness / 2.0, endcapRadiatorThickness / 2.0,
                                         endcapRadiatorSideLength / 2.0, endcapRadiatorSideLength2 / 2.0,
                                         endcapYLength / 2.0);
        std::string endcapRadiatorName = dd4hep::xml::_toString(sideID, "endcapRadiator%d");
        dd4hep::Volume endcapRadiatorVol(endcapRadiatorName, endcapRadiator, endcapRadiatorMaterial);

        double angle_degrees = shapeAngle * side;            // Calculate the angle for each side
        double angle_radians = angle_degrees * M_PI / 180.0; // radian angle for each side (rotating angle)

        double endcapRadiatorMidRadius = endcapRadiatorLayerInnerRadius + (endcapYLength / 2.0);

        double endcapRadiatorXOffset = endcapRadiatorMidRadius * std::cos(angle_radians + shapeAngle_radians);
        double endcapRadiatorYOffset = endcapRadiatorMidRadius * std::sin(angle_radians + shapeAngle_radians);

        dd4hep::RotationZ endcapRadiatorRotationZ(angle_radians + shapeAngle_radians);
        dd4hep::Rotation3D endcapRadiatorRotation = dd4hep::Rotation3D(endcapRadiatorRotationZ);

        double endcapRadiatorXPos = endcapRadiatorXOffset;
        double endcapRadiatorYPos = endcapRadiatorYOffset;
        double endcapRadiatorZPos = 0.0;

        dd4hep::Position endcapRadiatorTrans(endcapRadiatorXPos, endcapRadiatorYPos, endcapRadiatorZPos);
        dd4hep::PlacedVolume endcapRadiatorPhys = endcapRadiatorEnvVol.placeVolume(
            endcapRadiatorVol,
            dd4hep::Transform3D(endcapRadiatorRotation * dd4hep::RotationY(90. * dd4hep::degree), endcapRadiatorTrans));
        dd4hep::DetElement endcapRadiatorDE(endcapRadiatorEnvelopeDE, endcapRadiatorName + "DE", sideID);
        endcapRadiatorDE.setPlacement(endcapRadiatorPhys);
        endcapRadiatorVol.setVisAttributes(lcdd.visAttributes("yoke_vis"));
      }
    }
  }

  // -------------------------------------------------------------------------------------------
  return detElement;
}

DECLARE_DETELEMENT(muonSystemMuRWELL_o1_v01, createmuonSystemMuRWELL_o1_v01)