/*
 @author: Mahmoud Ali
mahmoud.ali@cern.ch

Factory for IDEA muon system
Expected xml structure (the 'sensitive' keyword is optional and defaults to false):
<detector type="muonSystem_o1_v00" ...>
  <dimensions x="..." y="..." z="..." z_offset="..." x_offset="..." y_offset="...">
  <sensitive type="SimpleTrackerSD"/>

  <!-- Specify the detector paramenters and the overlap // if you want exclude any component (e.g: endcap, radiators ..), just but its num =0 // radius is put in the middle, so its not the inner neither the outer  -->
  <barrelDetectorParameters numBarrelDetectorLayers ="... radius="... barrelLength="... numSides="... overlapY="... overlapZ="... />
  <barrelRadiatorParameters numBarrelRadiators="... barrelRadiatorThickness="... barrelRadiatorLayerRadius="... material=".../> 
  <endcapRadiatorParameters numEndcapRadiators="... endcapRadiatorThickness="...   endcapRadiatorLayerInnerRadius="... endcapRadiatorLayerOuterRadius="... endcapRadiatorZOffset="... material=".../>
  <endcapDetectorParameters numEndcapDetectorLayers="... endcapDetectorLayerInnerRadius="... endcapDetectorLayerOuterRadius="... endcapDetectorZOffset="... />


  <layer x="..." y="..." z="..." z_offset="..." x_offset="..." y_offset="..." material="...">
  . . . .
  <layer x="..." y="..." z="..." z_offset="..." x_offset="..." y_offset="..." material="..." sensitive="true">
</detector>

If used with sensitive layers, the readout must contain a "layer" field
*/


#include "DD4hep/DetFactoryHelper.h"
#include "XML/XMLElements.h"
#include <sstream> 
#include <iostream>
#include <cmath>

namespace det {


static dd4hep::Ref_t createmuonSystem_o1_v00(dd4hep::Detector& lcdd,
                                                           dd4hep::xml::Handle_t xmlElement,
                                                           dd4hep::SensitiveDetector sensDet) {
  dd4hep::xml::DetElement xmlDet = static_cast<dd4hep::xml::DetElement>(xmlElement);
  std::string name = xmlDet.nameStr();
  dd4hep::DetElement detElement(name, xmlDet.id());
  dd4hep::Volume experimentalHall = lcdd.pickMotherVolume(detElement);
  

  xml_comp_t dimensions(xmlDet.dimensions());

  dd4hep::Box barrelEnvelope(600 , 600 , 600);
  dd4hep::Material mat = lcdd.material("Air");
  dd4hep::Volume barrelVolume("Barrel", barrelEnvelope, mat);

// ----------------------------------------------------------------------------------------------------
//                         --- Barrel parameters ---

  auto barrelDetectorParameters = xmlElement.child(_Unicode(barrelDetectorParameters));
  int numBarrelDetectorLayers = barrelDetectorParameters.attr<double>("numBarrelDetectorLayers");
  double radius = barrelDetectorParameters.attr<double>("radius");
  double barrelLength = barrelDetectorParameters.attr<double>("barrelLength");
  int numSides = barrelDetectorParameters.attr<int>("numSides");
  double overlapY = barrelDetectorParameters.attr<double>("overlapY");
  double overlapZ = barrelDetectorParameters.attr<double>("overlapZ");
 // int numLayers = barrelDetectorParameters.attr<double>("numLayers");


  auto barrelRadiatorParameters = xmlElement.child(_Unicode(barrelRadiatorParameters));
  int numBarrelRadiators = barrelRadiatorParameters.attr<double>("numBarrelRadiators");
  double barrelRadiatorThickness = barrelRadiatorParameters.attr<double>("barrelRadiatorThickness");
  double barrelRadiatorLayerRadius = barrelRadiatorParameters.attr<double>("barrelRadiatorLayerRadius");
  dd4hep::Material barrelRadiatorMaterial = lcdd.material(barrelRadiatorParameters.attr<std::string>("material"));
 // double barrelRadiatorLayer2Radius = barrelDetectorParameters.attr<double>("barrelRadiatorLayer2Radius");


//                    --- Endcap parameters ---

  auto endcapDetectorParameters = xmlElement.child(_Unicode(endcapDetectorParameters));
  int numEndcapDetectorLayers = endcapDetectorParameters.attr<double>("numEndcapDetectorLayers");
  double endcapDetectorLayerInnerRadius = endcapDetectorParameters.attr<double>("endcapDetectorLayerInnerRadius");
  double endcapDetectorLayerOuterRadius = endcapDetectorParameters.attr<double>("endcapDetectorLayerOuterRadius");
  double endcapDetectorZOffset = endcapDetectorParameters.attr<double>("endcapDetectorZOffset");


  auto endcapRadiatorParameters = xmlElement.child(_Unicode(endcapRadiatorParameters));
  int numEndcapRadiators = endcapRadiatorParameters.attr<double>("numEndcapRadiators");
  double endcapRadiatorThickness = endcapRadiatorParameters.attr<double>("endcapRadiatorThickness");
  double endcapRadiatorLayerInnerRadius = endcapRadiatorParameters.attr<double>("endcapRadiatorLayerInnerRadius");
  double endcapRadiatorLayerOuterRadius = endcapRadiatorParameters.attr<double>("endcapRadiatorLayerOuterRadius");
  double endcapRadiatorZOffset = endcapRadiatorParameters.attr<double>("endcapRadiatorZOffset");
  dd4hep::Material endcapRadiatorMaterial = lcdd.material(endcapRadiatorParameters.attr<std::string>("material"));

  // Create layer boxes with their respective material, etc
//  auto layers = xmlElement.children(_Unicode(layer));
//  auto numLayers = xmlElement.numChildren(_Unicode(layer), true);
//  dd4hep::xml::Handle_t layer(layers.reset());
//  int sensitiveLayerIndex = 0;

  // -----------------------------------------------------------------------------------------------------------

   //   int numSides = 8; // Octagon has 8 sides
  double shapeAngle = 360.0 / numSides;
  double shapeAngle_radians = M_PI / numSides;  // Isn't it a half angle!!
  double sideLength = 2 * radius * std::tan(shapeAngle_radians);
  double barrelRadiatorSideLength = 2 * (barrelRadiatorLayerRadius - barrelRadiatorThickness/2.0) * std::tan(shapeAngle_radians);
  double barrelRadiatorSideLength2 = 2 * (barrelRadiatorLayerRadius + barrelRadiatorThickness/2.0) * std::tan(shapeAngle_radians);

  double endcapDetectorSideLength = 2 * (endcapDetectorLayerInnerRadius) * std::tan(shapeAngle_radians); 
  double endcapDetectorSideLength2 = 2 * (endcapDetectorLayerOuterRadius) * std::tan(shapeAngle_radians);

  double endcapRadiatorSideLength = 2 * (endcapRadiatorLayerInnerRadius) * std::tan(shapeAngle_radians); //it is also the same for the endcap detector layers.
  double endcapRadiatorSideLength2 = 2 * (endcapRadiatorLayerOuterRadius) * std::tan(shapeAngle_radians); //it is also the same for the endcap detector layers.

  double endcapDetectorYLength = endcapDetectorLayerOuterRadius - endcapDetectorLayerInnerRadius;
  double endcapYLength = endcapRadiatorLayerOuterRadius - endcapRadiatorLayerInnerRadius; // It is the distance betwwen the inner and the outer radius of the endcap, it can be in both Y and X dimensions //it is also the same for the endcap detector layers.

  int numCopiesY = sideLength / (2 * dimensions.y() - overlapY) + 1;
  int numCopiesZ = barrelLength / (2 * dimensions.z() - overlapZ) + 1;

  double remainderY = std::fmod(sideLength, (2 * dimensions.y() - overlapY)) / (2 * dimensions.y());
  double remainderZ = std::fmod(barrelLength, (2 * dimensions.z() - overlapZ)) / (2 * dimensions.z());

  double remainderYLength = remainderY * 2 * dimensions.y();
  double remainderZLength = remainderZ * 2 * dimensions.z();

  double totalOverlapValueY = (numCopiesY - 1) * overlapY;
  double totalOverlapValueZ = (numCopiesZ - 1) * overlapZ;

  double sideEnvX = 2 * dimensions.x();
  double sideEnvY = sideLength / 2.0;
  double sideEnvZ = barrelLength / 2.0;

  double barrelRadiatorEnvX = barrelRadiatorThickness/2.0;
  double barrelRadiatorEnvY = barrelRadiatorSideLength2 / 2.0;
  double barrelRadiatorEnvZ = barrelLength / 2.0;

  double endcapRadiatorEnvX = endcapRadiatorLayerOuterRadius;  // it depends on outer radius
  double endcapRadiatorEnvY = endcapRadiatorLayerOuterRadius;  // outer radius too
  double endcapRadiatorEnvZ = endcapRadiatorThickness / 2.0;  // layer thickness

  double endcapDetectorEnvX = endcapDetectorLayerOuterRadius;  // it depends on outer radius
  double endcapDetectorEnvY = endcapDetectorLayerOuterRadius;  // outer radius too
  double endcapDetectorEnvZ = 4 * dimensions.x();  // layer thickness

  double sideY = (numCopiesY -1) * 2 * dimensions.y() + remainderYLength - totalOverlapValueY;
  double sideZ = (numCopiesZ -1) * 2 * dimensions.z() + remainderZLength - totalOverlapValueZ;


  int nameCounter = 0;
  int idCounter = 0;

// ----------------------------------------------------------------------------------------------------
//                               // B A R R E L //
// ----------------------------------------------------------------------------------------------------

  if (numBarrelDetectorLayers > 0){

  for (int side = 0; side < numSides; ++side) {

      
      dd4hep::Box sideEnvelope(sideEnvX, sideEnvY, sideEnvZ);
      std::string sideName = dd4hep::xml::_toString(side, "side%d");
      dd4hep::Volume sideVol(sideName, sideEnvelope, mat);

      double angle_degrees = shapeAngle * side; // Calculate the angle for each side
      double angle_radians = angle_degrees * M_PI / 180.0;

      double sideXOffset = radius * std::cos(angle_radians);
      double sideYOffset = radius * std::sin(angle_radians);

      dd4hep::RotationZ sideRotationZ(angle_radians);
     //  dd4hep::Rotation3D sideRotation = rotationX;
      dd4hep::Rotation3D sideRotation = dd4hep::Rotation3D(sideRotationZ);

      double sideXPos = sideXOffset ; 
      double sideYPos = sideYOffset ;   
      double sideZPos = 0.0; 

      dd4hep::Position sideTrans(sideXPos, sideYPos, sideZPos);
      dd4hep::PlacedVolume sidePhys = experimentalHall.placeVolume(sideVol, dd4hep::Transform3D(sideRotation, sideTrans));
      sidePhys.addPhysVolID("side", side);
      dd4hep::DetElement sideDE(detElement, sideName + "DE", side);
      sideDE.setPlacement(sidePhys);
      sideVol.setVisAttributes(lcdd, xmlDet.visStr()); 

// -------------------------------------------------------------------------------------------------------      

           // Loop to create copies in the Y-direction
      for (int copyIndexY = 0; copyIndexY < numCopiesY; ++copyIndexY) {
 
      // Loop to create copies in the Z-direction
      for (int copyIndexZ = 0; copyIndexZ < numCopiesZ; ++copyIndexZ) {


        double yRotation = std::atan(dimensions.x() / dimensions.z());
        double zRotation = std::atan(dimensions.x() / dimensions.y());
        dd4hep::RotationY rotationY(yRotation);
        dd4hep::RotationZ rotationZ(zRotation);
        // Combine the Y and Z rotations to get a 3D rotation
        dd4hep::Rotation3D copyRotation = rotationY * rotationZ;

        // Calculate overlap values based on copyIndexY and copyIndexZ
        double overlapYValue = copyIndexY * overlapY;
        double overlapZValue = copyIndexZ * overlapZ;

        // Create layer boxes with their respective material, etc
        auto layers = xmlElement.children(_Unicode(layer));
        auto numLayers = xmlElement.numChildren(_Unicode(layer), true);
        int sensitiveLayerIndex = 0;
        // Loop through the layers
        dd4hep::xml::Handle_t layer(layers.reset());

        int chamberID = idCounter++;
            
        std::stringstream nameStream;
        nameStream << "envDE_" << nameCounter++;
        std::string chamberName = nameStream.str();


// --------------------------------------------------------------------------------------------

        
        if (copyIndexY == (numCopiesY - 1) || copyIndexZ == (numCopiesZ - 1)) {


            if (copyIndexY == (numCopiesY - 1) && copyIndexZ != (numCopiesZ - 1)) {            

              dd4hep::Box remainderYEnvelope(dimensions.x(), remainderY * dimensions.y(), dimensions.z());
              dd4hep::Volume remainderYEnvVolume(name + "remainderY", remainderYEnvelope, lcdd.material(dimensions.materialStr()));


              double remainderYEnvYPos = (copyIndexY * 2 * dimensions.y()) + (remainderY * dimensions.y()) - overlapYValue + dimensions.y_offset() - sideY/2.0;   
              double remainderYEnvZPos = (copyIndexZ * 2 * dimensions.z()) - overlapZValue + dimensions.z_offset() - sideZ/2.0 + dimensions.z(); 
              dd4hep::Position remainderYTrans(dimensions.x_offset(), remainderYEnvYPos, remainderYEnvZPos);
              dd4hep::PlacedVolume remainderYEnvPhys = sideVol.placeVolume(remainderYEnvVolume, dd4hep::Transform3D(copyRotation, remainderYTrans));
              remainderYEnvPhys.addPhysVolID("chamber", chamberID);
              dd4hep::DetElement envDE(sideDE, chamberName, chamberID);
              envDE.setPlacement(remainderYEnvPhys);
              remainderYEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());   


              for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
                dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
                dd4hep::Box layerShape(layerDet.x(), remainderY * layerDet.y(), layerDet.z());
                std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
                dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
                dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
                dd4hep::PlacedVolume layerPlacedVolume = remainderYEnvVolume.placeVolume(layerVolume, dd4hep::Transform3D(copyRotation, transLayer));
              //  dd4hep::DetElement layerDE(envDE, "layerDE", layerIdx);
              //  layerDE.setPlacement(layerPlacedVolume);

                if (layer.hasAttr("vis")) {
                  layerVolume.setVisAttributes(lcdd, layerDet.visStr());
                }
                if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
                  dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                  sensDet.setType(sdType.typeStr());
                  layerVolume.setSensitiveDetector(sensDet);
                  layerPlacedVolume.addPhysVolID("layer", sensitiveLayerIndex);
                  sensitiveLayerIndex++;
                }
                layer.m_node = layers.next();
              }

            } 

// -----------------------------------------------

            if (copyIndexZ == (numCopiesZ - 1) && copyIndexY != (numCopiesY - 1)) {

              dd4hep::Box remainderZEnvelope(dimensions.x(), dimensions.y(),  remainderZ * dimensions.z());
              dd4hep::Volume remainderZEnvVolume(name + "remainderZ", remainderZEnvelope, lcdd.material(dimensions.materialStr()));


              double remainderZEnvYPos = (copyIndexY * 2 * dimensions.y()) - overlapYValue + dimensions.y_offset() - sideY/2.0 + dimensions.y();   
              double remainderZEnvZPos = (copyIndexZ * 2 * dimensions.z()) + (remainderZ * dimensions.z()) - overlapZValue + dimensions.z_offset() - sideZ/2.0; 
              dd4hep::Position remainderZTrans(dimensions.x_offset(), remainderZEnvYPos, remainderZEnvZPos);
              dd4hep::PlacedVolume remainderZEnvPhys = sideVol.placeVolume(remainderZEnvVolume, dd4hep::Transform3D(copyRotation, remainderZTrans));
              remainderZEnvPhys.addPhysVolID("chamber", chamberID);
              dd4hep::DetElement envDE(sideDE, chamberName, chamberID);
              envDE.setPlacement(remainderZEnvPhys);
              remainderZEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());   


              for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
                dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
                dd4hep::Box layerShape(layerDet.x(), layerDet.y(), remainderZ * layerDet.z());
                std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
                dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
                dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
                dd4hep::PlacedVolume layerPlacedVolume = remainderZEnvVolume.placeVolume(layerVolume, dd4hep::Transform3D(copyRotation, transLayer));
              //  dd4hep::DetElement layerDE(envDE, "layerDE", layerIdx);
              //  layerDE.setPlacement(layerPlacedVolume);

                if (layer.hasAttr("vis")) {
                  layerVolume.setVisAttributes(lcdd, layerDet.visStr());
                }
                if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
                  dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                  sensDet.setType(sdType.typeStr());
                  layerVolume.setSensitiveDetector(sensDet);
                  layerPlacedVolume.addPhysVolID("layer", sensitiveLayerIndex);
                  sensitiveLayerIndex++;
                }
                layer.m_node = layers.next();
              }


            }
        
// ------------------------------------------------------------

            if (copyIndexY == (numCopiesY - 1) && copyIndexZ == (numCopiesZ - 1)) {            

              dd4hep::Box remainderYZEnvelope(dimensions.x(), remainderY * dimensions.y(), remainderZ * dimensions.z());
              dd4hep::Volume remainderYZEnvVolume(name + "remainderYZ", remainderYZEnvelope, lcdd.material(dimensions.materialStr()));


              double remainderYZEnvYPos = (copyIndexY * 2 * dimensions.y()) + (remainderY * dimensions.y()) - overlapYValue + dimensions.y_offset() - sideY/2.0;   
              double remainderYZEnvZPos = (copyIndexZ * 2 * dimensions.z()) + (remainderZ * dimensions.z()) - overlapZValue + dimensions.z_offset() - sideZ/2.0; 
              dd4hep::Position remainderYZTrans(dimensions.x_offset(), remainderYZEnvYPos, remainderYZEnvZPos);
              dd4hep::PlacedVolume remainderYZEnvPhys = sideVol.placeVolume(remainderYZEnvVolume, dd4hep::Transform3D(copyRotation, remainderYZTrans));
              remainderYZEnvPhys.addPhysVolID("chamber", chamberID);
              dd4hep::DetElement envDE(sideDE, chamberName, chamberID);
              envDE.setPlacement(remainderYZEnvPhys);
              remainderYZEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());   


              for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
                dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
                dd4hep::Box layerShape(layerDet.x(), remainderY * layerDet.y(), remainderZ * layerDet.z());
                std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
                dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
                dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
                dd4hep::PlacedVolume layerPlacedVolume = remainderYZEnvVolume.placeVolume(layerVolume, dd4hep::Transform3D(copyRotation, transLayer));
              //  dd4hep::DetElement layerDE(envDE, "layerDE", layerIdx);
              //  layerDE.setPlacement(layerPlacedVolume);

                if (layer.hasAttr("vis")) {
                  layerVolume.setVisAttributes(lcdd, layerDet.visStr());
                }
                if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
                  dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                  sensDet.setType(sdType.typeStr());
                  layerVolume.setSensitiveDetector(sensDet);
                  layerPlacedVolume.addPhysVolID("layer", sensitiveLayerIndex);
                  sensitiveLayerIndex++;
                }
                layer.m_node = layers.next();
              }

            } 

        }
// -------------------------------------------------------------

        else{

        dd4hep::Box envelope(dimensions.x(), dimensions.y(), dimensions.z());
        dd4hep::Volume envVolume(name, envelope, lcdd.material(dimensions.materialStr()));




        double envYPos = (copyIndexY * 2 * dimensions.y()) - overlapYValue + dimensions.y_offset() - sideY/2.0 + dimensions.y();   
        double envZPos = (copyIndexZ * 2 * dimensions.z()) - overlapZValue + dimensions.z_offset() - sideZ/2.0 + dimensions.z(); 
        dd4hep::Position trans(dimensions.x_offset(), envYPos, envZPos);
        dd4hep::PlacedVolume envPhys = sideVol.placeVolume(envVolume, dd4hep::Transform3D(copyRotation, trans));
        envPhys.addPhysVolID("chamber", chamberID);
        dd4hep::DetElement envDE(sideDE, chamberName, chamberID);
        envDE.setPlacement(envPhys);
        envVolume.setVisAttributes(lcdd, xmlDet.visStr());


       // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() << std::endl;

      for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
        dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
        dd4hep::Box layerShape(layerDet.x(), layerDet.y(), layerDet.z());
        std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
        dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
        dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
        dd4hep::PlacedVolume layerPlacedVolume = envVolume.placeVolume(layerVolume, dd4hep::Transform3D(copyRotation, transLayer));
      //  dd4hep::DetElement layerDE(envDE, "layerDE", layerIdx);
      //  layerDE.setPlacement(layerPlacedVolume);

        if (layer.hasAttr("vis")) {
          layerVolume.setVisAttributes(lcdd, layerDet.visStr());
        }
        if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
          dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
          sensDet.setType(sdType.typeStr());
          layerVolume.setSensitiveDetector(sensDet);
          layerPlacedVolume.addPhysVolID("layer", sensitiveLayerIndex);
          sensitiveLayerIndex++;
        }
        layer.m_node = layers.next();
        }


      
     }


    }
    }

  
    }
  }
// ---------------------------------------// R A D I A T O R S //------------------------------------------------

  if ( numBarrelRadiators > 0) {


    for (int side = 0; side < numSides; ++side) {

      dd4hep::Box barrelRadiatorEnvelope(barrelRadiatorEnvX, barrelRadiatorEnvY, barrelRadiatorEnvZ);
      std::string barrelRadiatorEnvelopeName = dd4hep::xml::_toString(side, "barrelRadiatorEnvelope%d");
      dd4hep::Volume barrelRadiatorEnvVol(barrelRadiatorEnvelopeName, barrelRadiatorEnvelope, mat);

      double angle_degrees = shapeAngle * side; // Calculate the angle for each side
      double angle_radians = angle_degrees * M_PI / 180.0;

      double barrelRadiatorXOffset = barrelRadiatorLayerRadius * std::cos(angle_radians);
      double barrelRadiatorYOffset = barrelRadiatorLayerRadius * std::sin(angle_radians);

      dd4hep::RotationZ barrelRadiatorRotationZ(angle_radians);
      dd4hep::Rotation3D barrelRadiatorRotation = dd4hep::Rotation3D(barrelRadiatorRotationZ);

      double barrelRadiatorXPos = barrelRadiatorXOffset ; 
      double barrelRadiatorYPos = barrelRadiatorYOffset ;   
      double barrelRadiatorZPos = 0.0; 

      dd4hep::Position barrelRadiatorEnvelopeTrans(barrelRadiatorXPos, barrelRadiatorYPos, barrelRadiatorZPos);
      dd4hep::PlacedVolume barrelRadiatorEnvelopePhys = experimentalHall.placeVolume(barrelRadiatorEnvVol, dd4hep::Transform3D(barrelRadiatorRotation, barrelRadiatorEnvelopeTrans));
      barrelRadiatorEnvelopePhys.addPhysVolID("barrelRadiatorEnvelope", side);
      dd4hep::DetElement barrelRadiatorEnvelopeDE(detElement, barrelRadiatorEnvelopeName + "DE", side);
      barrelRadiatorEnvelopeDE.setPlacement(barrelRadiatorEnvelopePhys);
      barrelRadiatorEnvVol.setVisAttributes(lcdd, xmlDet.visStr()); 

// -----------------------------------------------------------

      dd4hep::Trapezoid barrelRadiator(barrelLength/2.0, barrelLength/2.0, barrelRadiatorSideLength/2.0, barrelRadiatorSideLength2/2.0, barrelRadiatorThickness/2.0);
      std::string barrelRadiatorName = dd4hep::xml::_toString(side, "barrelRadiator%d");
      dd4hep::Volume barrelRadiatorVol(barrelRadiatorName, barrelRadiator, barrelRadiatorMaterial);


      dd4hep::Position barrelRadiatorTrans(0, 0, 0);
      dd4hep::PlacedVolume barrelRadiatorPhys = barrelRadiatorEnvVol.placeVolume(barrelRadiatorVol, dd4hep::Transform3D(dd4hep::RotationY(90.* dd4hep::degree), barrelRadiatorTrans));
      barrelRadiatorPhys.addPhysVolID("barrelRadiator", side);
      dd4hep::DetElement barrelRadiatorDE(barrelRadiatorEnvelopeDE, barrelRadiatorName + "DE", side);
      barrelRadiatorDE.setPlacement(barrelRadiatorPhys);
      barrelRadiatorVol.setVisAttributes(lcdd.visAttributes("yoke_vis"));       


    }  
  }

// -------------------------------------------------------------------------------------------
//                                    // E N D C A P //
// -------------------------------------------------------------------------------------------
//---------------------------------------Detectors--------------------------------------------


  if ( numEndcapDetectorLayers > 0) {

    dd4hep::Box endcapDetectorEnvelope(endcapDetectorEnvX, endcapDetectorEnvY, endcapDetectorEnvZ);
    std::string endcapDetectorEnvelopeName = dd4hep::xml::_toString(numEndcapDetectorLayers, "endcapDetectorEnvelope%d");
    dd4hep::Volume endcapDetectorEnvVol(endcapDetectorEnvelopeName, endcapDetectorEnvelope, mat);



    double endcapDetectorEnvXPos = 0.0 ; 
    double endcapDetectorEnvYPos = 0.0 ;   
    double endcapDetectorEnvZPos = endcapDetectorZOffset; 

    dd4hep::Position endcapDetectorEnvelopeTrans(endcapDetectorEnvXPos, endcapDetectorEnvYPos, endcapDetectorEnvZPos);
    dd4hep::PlacedVolume endcapDetectorEnvelopePhys = experimentalHall.placeVolume(endcapDetectorEnvVol, dd4hep::Transform3D(dd4hep::RotationZ(0.), endcapDetectorEnvelopeTrans));
    endcapDetectorEnvelopePhys.addPhysVolID("endcapDetectorEnvelope", numEndcapDetectorLayers);
    dd4hep::DetElement endcapDetectorEnvelopeDE(detElement, endcapDetectorEnvelopeName + "DE", numEndcapDetectorLayers); // remember to loop over numEndcapDetectorLayers.. because now it still just one number.
    endcapDetectorEnvelopeDE.setPlacement(endcapDetectorEnvelopePhys);
    endcapDetectorEnvVol.setVisAttributes(lcdd, xmlDet.visStr()); 


// -----------------------------------------------

    for (int side = 0; side < numSides; ++side) {

      dd4hep::Trapezoid endcapDetectorSideEnvelope(endcapDetectorEnvZ/2.0, endcapDetectorEnvZ/2.0, endcapDetectorSideLength/2.0, endcapDetectorSideLength2/2.0, endcapDetectorYLength/2.0); //maybe "endcapDetectorEnvZ" shouldn't be divided over 2.
      std::string endcapDetectorSideEnvName = dd4hep::xml::_toString(side, "endcapDetectorSideEnv%d");
      dd4hep::Volume endcapDetectorSideEnvVol(endcapDetectorSideEnvName, endcapDetectorSideEnvelope, mat);


      double angle_degrees = shapeAngle * side; // Calculate the angle for each side
      double angle_radians = angle_degrees * M_PI / 180.0;

      double endcapDetectorMidRadius = endcapDetectorLayerInnerRadius + (endcapDetectorYLength /2.0);

      double endcapDetectorXOffset = endcapDetectorMidRadius * std::cos(angle_radians);
      double endcapDetectorYOffset = endcapDetectorMidRadius * std::sin(angle_radians);

      dd4hep::RotationZ endcapDetectorRotationZ(angle_radians);
      dd4hep::Rotation3D endcapDetectorRotation = dd4hep::Rotation3D(endcapDetectorRotationZ);

      double endcapDetectorSideEnvXPos = endcapDetectorXOffset ; 
      double endcapDetectorSideEnvYPos = endcapDetectorYOffset ;   
      double endcapDetectorSideEnvZPos = - endcapDetectorEnvZ/2.0 ;   
      double endcapDetectorSideEnvZPos2 = endcapDetectorEnvZ/2.0;

// ---------  here I start to divide the two z-positions //maybe by odd and even numbers

      if (side % 2 == 0) {

        dd4hep::Position endcapDetectorSideEnvTrans(endcapDetectorSideEnvXPos, endcapDetectorSideEnvYPos, endcapDetectorSideEnvZPos);
        dd4hep::PlacedVolume endcapDetectorSideEnvPhys = endcapDetectorEnvVol.placeVolume(endcapDetectorSideEnvVol, dd4hep::Transform3D(endcapDetectorRotation * dd4hep::RotationY(90.* dd4hep::degree), endcapDetectorSideEnvTrans));
        endcapDetectorSideEnvPhys.addPhysVolID("endcapDetectorSideEnvelope", side);
        dd4hep::DetElement endcapDetectorSideEnvDE(endcapDetectorEnvelopeDE, endcapDetectorSideEnvName + "DE", side);
        endcapDetectorSideEnvDE.setPlacement(endcapDetectorSideEnvPhys);
        endcapDetectorSideEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

        // ----- dividing the trapezoid envelope to smaller pieces (rectangles)

        int numChambersY = endcapDetectorYLength / (2 * dimensions.y() - overlapY); // numbers of the rectangles in each teapezoids.. depends on the number of the chambers that can be overlapped in Y direction
        
        for (int rectangle = 0; rectangle < numChambersY; ++rectangle) {

          double rectangleEnvX = endcapDetectorEnvZ/4.0;
          double rectangleEnvY = (endcapDetectorLayerOuterRadius - rectangle * (2 * dimensions.y() - overlapY)) * std::tan(shapeAngle_radians); // without multiplying by 2 .. because it is the half length // it should be dimensions.x() instead of z, but in the endcap its perpendicular to usual dimension set
          double rectangleEnvZ = dimensions.z();

          dd4hep::Box rectangleEnvelope(rectangleEnvX, rectangleEnvY, rectangleEnvZ);
          std::string rectangleEnvelopeName = dd4hep::xml::_toString(rectangle, "rectangleEnvelope%d");
          dd4hep::Volume rectangleEnvVol(rectangleEnvelopeName, rectangleEnvelope, mat);

          double rectangleEnvXPos = 0.0; 
          double rectangleEnvYPos = 0.0;   
          double rectangleEnvZPos = endcapDetectorYLength/2.0 - dimensions.y() - (rectangle * (2 * dimensions.y() - overlapY));

          double yRotation = std::atan(dimensions.x() / dimensions.y());
          dd4hep::RotationY rotationY(yRotation);
       //   dd4hep::Rotation3D overlapRotation = rotationX;

     //     if (rectangle == numChambersY) {} // To make the last rectangle overlapping over the previous one, so, it fits on the inner radius.

          dd4hep::Position rectangeEnvelopeTrans(rectangleEnvXPos, rectangleEnvYPos, rectangleEnvZPos);
          dd4hep::PlacedVolume rectangleEnvelopePhys = endcapDetectorSideEnvVol.placeVolume(rectangleEnvVol, dd4hep::Transform3D(rotationY, rectangeEnvelopeTrans));
          rectangleEnvelopePhys.addPhysVolID("rectangleEnvelope", rectangle);
          dd4hep::DetElement rectangleEnvelopeDE(endcapDetectorSideEnvDE, rectangleEnvelopeName + "DE", rectangle); // remember to loop over numEndcapDetectorLayers.. because now it still just one number.
          rectangleEnvelopeDE.setPlacement(rectangleEnvelopePhys);
          rectangleEnvVol.setVisAttributes(lcdd, xmlDet.visStr());
        
          // ------------------------ start to build the chamber envelopes -------------------


          int numChambersInRectangle = 2 * rectangleEnvY / (2 * dimensions.y() - overlapY); // number of the chambers in each rectangle

          for (int chamberIndex = 0; chamberIndex < (numChambersInRectangle + 1); chamberIndex++) {

          dd4hep::Box envelope(dimensions.x(), dimensions.y(), dimensions.z());
          dd4hep::Volume envVolume(name, envelope, lcdd.material(dimensions.materialStr())); // "name" has to be changed

          double rectangleRemainderY = std::fmod(2 * rectangleEnvY, (2 * dimensions.y() - overlapY)) / (2 * dimensions.y());
          //double rectangleRemainderYLength = rectangleRemainderY * 2 * dimensions.y();

          dd4hep::Box rectangleRemainderYEnvelope(dimensions.x(), rectangleRemainderY * dimensions.y(), dimensions.z());
          dd4hep::Volume rectangleRemainderYEnvVolume(name + "rectangleRemainderY", rectangleRemainderYEnvelope, lcdd.material(dimensions.materialStr()));          

          double envYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() - rectangleEnvY + dimensions.y();
          double rectangleRemainderREnvYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() - rectangleEnvY + rectangleRemainderY * dimensions.y();
          
          double zRotation = std::atan(dimensions.x() / dimensions.z());
          dd4hep::RotationZ chamberRotation(zRotation);

          int endcapChamberID = side * 1000 + rectangle * 100 + chamberIndex;  // later you should add layer number to distinguish between different endcap layers.

          auto layers = xmlElement.children(_Unicode(layer));
          auto numLayers = xmlElement.numChildren(_Unicode(layer), true);
          dd4hep::xml::Handle_t layer(layers.reset());
          int sensitiveLayerIndex = 0;

          std::stringstream nameStream;
          nameStream << "envDE_" << endcapChamberID;
          std::string endcapChamberName = nameStream.str();

          // --- two cases : one for full chambers ----------------------------------------------

          if (chamberIndex == numChambersInRectangle) {

          dd4hep::Position rectangleRemainderTrans(dimensions.x_offset(), rectangleRemainderREnvYPos, 0.0);
          dd4hep::PlacedVolume rectangleRemainderEnvPhys = rectangleEnvVol.placeVolume(rectangleRemainderYEnvVolume, dd4hep::Transform3D(chamberRotation, rectangleRemainderTrans));
          rectangleRemainderEnvPhys.addPhysVolID("chamber", endcapChamberID);
          dd4hep::DetElement rectangleRemainderEnvDE(rectangleEnvelopeDE, endcapChamberName, endcapChamberID);
          rectangleRemainderEnvDE.setPlacement(rectangleRemainderEnvPhys);
          rectangleRemainderYEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());


         // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() << std::endl;

        for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
          dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
          dd4hep::Box layerShape(layerDet.x(), rectangleRemainderY * layerDet.y(), layerDet.z());
          std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
          dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
          dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
          dd4hep::PlacedVolume layerPlacedVolume = rectangleRemainderYEnvVolume.placeVolume(layerVolume, dd4hep::Transform3D(chamberRotation, transLayer));
        //  dd4hep::DetElement layerDE(envDE, "layerDE", layerIdx);
        //  layerDE.setPlacement(layerPlacedVolume);

          if (layer.hasAttr("vis")) {
            layerVolume.setVisAttributes(lcdd, layerDet.visStr());
          }
          if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
            dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
            sensDet.setType(sdType.typeStr());
            layerVolume.setSensitiveDetector(sensDet);
            layerPlacedVolume.addPhysVolID("layer", sensitiveLayerIndex);
            sensitiveLayerIndex++;
          }
        layer.m_node = layers.next();
        }




          // ----------------
          } else {

          dd4hep::Position trans(dimensions.x_offset(), envYPos, 0.0);
          dd4hep::PlacedVolume envPhys = rectangleEnvVol.placeVolume(envVolume, dd4hep::Transform3D(chamberRotation, trans));
          envPhys.addPhysVolID("chamber", endcapChamberID);
          dd4hep::DetElement envDE(rectangleEnvelopeDE, endcapChamberName, endcapChamberID);
          envDE.setPlacement(envPhys);
          envVolume.setVisAttributes(lcdd, xmlDet.visStr());


         // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() << std::endl;

        for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
          dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
          dd4hep::Box layerShape(layerDet.x(), layerDet.y(), layerDet.z());
          std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
          dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
          dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
          dd4hep::PlacedVolume layerPlacedVolume = envVolume.placeVolume(layerVolume, dd4hep::Transform3D(chamberRotation, transLayer));
        //  dd4hep::DetElement layerDE(envDE, "layerDE", layerIdx);
        //  layerDE.setPlacement(layerPlacedVolume);

          if (layer.hasAttr("vis")) {
            layerVolume.setVisAttributes(lcdd, layerDet.visStr());
          }
          if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
            dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
            sensDet.setType(sdType.typeStr());
            layerVolume.setSensitiveDetector(sensDet);
            layerPlacedVolume.addPhysVolID("layer", sensitiveLayerIndex);
            sensitiveLayerIndex++;
          }
        layer.m_node = layers.next();
        }

        }
        }          

        }
// -----------------------------------------------------------------------------
      } else {

        dd4hep::Position endcapDetectorSideEnvTrans(endcapDetectorSideEnvXPos, endcapDetectorSideEnvYPos, endcapDetectorSideEnvZPos2);
        dd4hep::PlacedVolume endcapDetectorSideEnvPhys = endcapDetectorEnvVol.placeVolume(endcapDetectorSideEnvVol, dd4hep::Transform3D(endcapDetectorRotation * dd4hep::RotationY(90.* dd4hep::degree), endcapDetectorSideEnvTrans));
        endcapDetectorSideEnvPhys.addPhysVolID("endcapDetectorSideEnvelope", side);
        dd4hep::DetElement endcapDetectorSideEnvDE(endcapDetectorEnvelopeDE, endcapDetectorSideEnvName + "DE", side);
        endcapDetectorSideEnvDE.setPlacement(endcapDetectorSideEnvPhys);
        endcapDetectorSideEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

        // ----- dividing the trapezoid envelope to smaller pieces (rectangles)

        int numChambersY = endcapDetectorYLength / (2 * dimensions.y() - overlapY); // numbers of the rectangles in each teapezoids.. depends on the number of the chambers that can be overlapped in Y direction
        
        for (int rectangle = 0; rectangle < numChambersY; ++rectangle) {

          double rectangleEnvX = endcapDetectorEnvZ/4.0;
          double rectangleEnvY = (endcapDetectorLayerOuterRadius - rectangle * (2 * dimensions.z() - overlapY)) * std::tan(shapeAngle_radians); // without multiplying by 2 .. because it is the half length // it should be dimensions.x() instead of z, but in the endcap its perpendicular to usual dimension set
          double rectangleEnvZ = dimensions.y();

          dd4hep::Box rectangleEnvelope(rectangleEnvX, rectangleEnvY, rectangleEnvZ);
          std::string rectangleEnvelopeName = dd4hep::xml::_toString(rectangle, "rectangleEnvelope%d");
          dd4hep::Volume rectangleEnvVol(rectangleEnvelopeName, rectangleEnvelope, mat);

          double rectangleEnvXPos = 0.0 ; 
          double rectangleEnvYPos = 0.0;   
          double rectangleEnvZPos = endcapDetectorYLength/2.0 - dimensions.y() - (rectangle * (2 * dimensions.y() - overlapY));

          double yRotation = std::atan(dimensions.x() / dimensions.y());
          dd4hep::RotationY rotationY(yRotation);
       //   dd4hep::Rotation3D overlapRotation = rotationX;

          dd4hep::Position rectangeEnvelopeTrans(rectangleEnvXPos, rectangleEnvYPos, rectangleEnvZPos);
          dd4hep::PlacedVolume rectangleEnvelopePhys = endcapDetectorSideEnvVol.placeVolume(rectangleEnvVol, dd4hep::Transform3D(rotationY, rectangeEnvelopeTrans));
          rectangleEnvelopePhys.addPhysVolID("rectangleEnvelope", rectangle);
          dd4hep::DetElement rectangleEnvelopeDE(endcapDetectorSideEnvDE, rectangleEnvelopeName + "DE", rectangle); // remember to loop over numEndcapDetectorLayers.. because now it still just one number.
          rectangleEnvelopeDE.setPlacement(rectangleEnvelopePhys);
          rectangleEnvVol.setVisAttributes(lcdd, xmlDet.visStr());
        
          // ------------------------ start to build the chamber envelopes -------------------

          int numChambersInRectangle = 2 * rectangleEnvY / (2 * dimensions.y() - overlapY); // number of the chambers in each rectangle

          for (int chamberIndex = 0; chamberIndex < (numChambersInRectangle + 1); chamberIndex++) {

          dd4hep::Box envelope(dimensions.x(), dimensions.y(), dimensions.z());
          dd4hep::Volume envVolume(name, envelope, lcdd.material(dimensions.materialStr())); // "name" has to be changed

          double rectangleRemainderY = std::fmod(2 * rectangleEnvY, (2 * dimensions.y() - overlapY)) / (2 * dimensions.y());
          //double rectangleRemainderYLength = rectangleRemainderY * 2 * dimensions.y();

          dd4hep::Box rectangleRemainderYEnvelope(dimensions.x(), rectangleRemainderY * dimensions.y(), dimensions.z());
          dd4hep::Volume rectangleRemainderYEnvVolume(name + "rectangleRemainderY", rectangleRemainderYEnvelope, lcdd.material(dimensions.materialStr()));          

          double envYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() - rectangleEnvY + dimensions.y();
          double rectangleRemainderREnvYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() - rectangleEnvY + rectangleRemainderY * dimensions.y();
          
          double zRotation = std::atan(dimensions.x() / dimensions.z());
          dd4hep::RotationZ chamberRotation(zRotation);

          int endcapChamberID = side * 1000 + rectangle * 100 + chamberIndex;  // later you should add layer number to distinguish between different endcap layers.

          auto layers = xmlElement.children(_Unicode(layer));
          auto numLayers = xmlElement.numChildren(_Unicode(layer), true);
          dd4hep::xml::Handle_t layer(layers.reset());
          int sensitiveLayerIndex = 0;

          std::stringstream nameStream;
          nameStream << "envDE_" << endcapChamberID;
          std::string endcapChamberName = nameStream.str();

          // --- two cases : one for full chambers ----------------------------------------------

          if (chamberIndex == numChambersInRectangle) {

          dd4hep::Position rectangleRemainderTrans(dimensions.x_offset(), rectangleRemainderREnvYPos, 0.0);
          dd4hep::PlacedVolume rectangleRemainderEnvPhys = rectangleEnvVol.placeVolume(rectangleRemainderYEnvVolume, dd4hep::Transform3D(chamberRotation, rectangleRemainderTrans));
          rectangleRemainderEnvPhys.addPhysVolID("chamber", endcapChamberID);
          dd4hep::DetElement rectangleRemainderEnvDE(rectangleEnvelopeDE, endcapChamberName, endcapChamberID);
          rectangleRemainderEnvDE.setPlacement(rectangleRemainderEnvPhys);
          rectangleRemainderYEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());


         // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() << std::endl;

        for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
          dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
          dd4hep::Box layerShape(layerDet.x(), rectangleRemainderY * layerDet.y(), layerDet.z());
          std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
          dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
          dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
          dd4hep::PlacedVolume layerPlacedVolume = rectangleRemainderYEnvVolume.placeVolume(layerVolume, dd4hep::Transform3D(chamberRotation, transLayer));
        //  dd4hep::DetElement layerDE(envDE, "layerDE", layerIdx);
        //  layerDE.setPlacement(layerPlacedVolume);

          if (layer.hasAttr("vis")) {
            layerVolume.setVisAttributes(lcdd, layerDet.visStr());
          }
          if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
            dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
            sensDet.setType(sdType.typeStr());
            layerVolume.setSensitiveDetector(sensDet);
            layerPlacedVolume.addPhysVolID("layer", sensitiveLayerIndex);
            sensitiveLayerIndex++;
          }
        layer.m_node = layers.next();
        }




          // ----------------
          } else {

          dd4hep::Position trans(dimensions.x_offset(), envYPos, 0.0);
          dd4hep::PlacedVolume envPhys = rectangleEnvVol.placeVolume(envVolume, dd4hep::Transform3D(chamberRotation, trans));
          envPhys.addPhysVolID("chamber", endcapChamberID);
          dd4hep::DetElement envDE(rectangleEnvelopeDE, endcapChamberName, endcapChamberID);
          envDE.setPlacement(envPhys);
          envVolume.setVisAttributes(lcdd, xmlDet.visStr());


         // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() << std::endl;

        for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
          dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
          dd4hep::Box layerShape(layerDet.x(), layerDet.y(), layerDet.z());
          std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
          dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
          dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
          dd4hep::PlacedVolume layerPlacedVolume = envVolume.placeVolume(layerVolume, dd4hep::Transform3D(chamberRotation, transLayer));
        //  dd4hep::DetElement layerDE(envDE, "layerDE", layerIdx);
        //  layerDE.setPlacement(layerPlacedVolume);

          if (layer.hasAttr("vis")) {
            layerVolume.setVisAttributes(lcdd, layerDet.visStr());
          }
          if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
            dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
            sensDet.setType(sdType.typeStr());
            layerVolume.setSensitiveDetector(sensDet);
            layerPlacedVolume.addPhysVolID("layer", sensitiveLayerIndex);
            sensitiveLayerIndex++;
          }
        layer.m_node = layers.next();
        }

        }
        }          

        }


      }      

    }
  }




// --------------------------------------Radiators--------------------------------------------
  if ( numEndcapRadiators > 0) {

    dd4hep::Box endcapRadiatorEnvelope(endcapRadiatorEnvX, endcapRadiatorEnvY, endcapRadiatorEnvZ);
    std::string endcapRadiatorEnvelopeName = dd4hep::xml::_toString(numEndcapRadiators, "endcapRadiatorEnvelope%d");
    dd4hep::Volume endcapRadiatorEnvVol(endcapRadiatorEnvelopeName, endcapRadiatorEnvelope, mat);



    double endcapRadiatorEnvXPos = 0.0 ; 
    double endcapRadiatorEnvYPos = 0.0 ;   
    double endcapRadiatorEnvZPos = endcapRadiatorZOffset; 

    dd4hep::Position endcapRadiatorEnvelopeTrans(endcapRadiatorEnvXPos, endcapRadiatorEnvYPos, endcapRadiatorEnvZPos);
    dd4hep::PlacedVolume endcapRadiatorEnvelopePhys = experimentalHall.placeVolume(endcapRadiatorEnvVol, dd4hep::Transform3D(dd4hep::RotationZ(0.), endcapRadiatorEnvelopeTrans));
    endcapRadiatorEnvelopePhys.addPhysVolID("endcapRadiatorEnvelope", numEndcapRadiators);
    dd4hep::DetElement endcapRadiatorEnvelopeDE(detElement, endcapRadiatorEnvelopeName + "DE", numEndcapRadiators);
    endcapRadiatorEnvelopeDE.setPlacement(endcapRadiatorEnvelopePhys);
    endcapRadiatorEnvVol.setVisAttributes(lcdd, xmlDet.visStr()); 


// ---------------------------------------------

    for (int side = 0; side < numSides; ++side) {

      dd4hep::Trapezoid endcapRadiator(endcapRadiatorThickness/2.0, endcapRadiatorThickness/2.0, endcapRadiatorSideLength/2.0, endcapRadiatorSideLength2/2.0, endcapYLength/2.0);
      std::string endcapRadiatorName = dd4hep::xml::_toString(side, "endcapRadiator%d");
      dd4hep::Volume endcapRadiatorVol(endcapRadiatorName, endcapRadiator, endcapRadiatorMaterial);


      double angle_degrees = shapeAngle * side; // Calculate the angle for each side
      double angle_radians = angle_degrees * M_PI / 180.0;

      double endcapRadiatorMidRadius = endcapRadiatorLayerInnerRadius + (endcapYLength /2.0);

      double endcapRadiatorXOffset = endcapRadiatorMidRadius * std::cos(angle_radians);
      double endcapRadiatorYOffset = endcapRadiatorMidRadius * std::sin(angle_radians);

      dd4hep::RotationZ endcapRadiatorRotationZ(angle_radians);
      dd4hep::Rotation3D endcapRadiatorRotation = dd4hep::Rotation3D(endcapRadiatorRotationZ);

      double endcapRadiatorXPos = endcapRadiatorXOffset ; 
      double endcapRadiatorYPos = endcapRadiatorYOffset ;   
      double endcapRadiatorZPos = 0.0; 


      dd4hep::Position endcapRadiatorTrans(endcapRadiatorXPos, endcapRadiatorYPos, endcapRadiatorZPos);
      dd4hep::PlacedVolume endcapRadiatorPhys = endcapRadiatorEnvVol.placeVolume(endcapRadiatorVol, dd4hep::Transform3D(endcapRadiatorRotation * dd4hep::RotationY(90.* dd4hep::degree), endcapRadiatorTrans));
      endcapRadiatorPhys.addPhysVolID("endcapRadiator", side);
      dd4hep::DetElement endcapRadiatorDE(endcapRadiatorEnvelopeDE, endcapRadiatorName + "DE", side);
      endcapRadiatorDE.setPlacement(endcapRadiatorPhys);
      endcapRadiatorVol.setVisAttributes(lcdd.visAttributes("yoke_vis"));       

    
    }
  }



// ------------------------------------------------------------------------------------------- 
  dd4hep::PlacedVolume barrelPhys = experimentalHall.placeVolume(barrelVolume);
  barrelPhys.addPhysVolID("system", xmlDet.id());
  detElement.setPlacement(barrelPhys);
  barrelVolume.setVisAttributes(lcdd.visAttributes("no_vis")); 
  return detElement;
}
}

DECLARE_DETELEMENT(muonSystem_o1_v00, det::createmuonSystem_o1_v00)
