/*
 @author: Mahmoud Ali
mahmoud.ali@cern.ch

Factory for IDEA muon system
Expected xml structure (the 'sensitive' keyword is optional and defaults to false):
<detector type="muonSystemMuRWELL_o1_v01" ...>
  <envelope rmin="..." rmax="..." z="..."  z_offset="..." material="..."/>  <!-- dimensions of the whole detector envelope-->
  <dimensions x="..." y="..." z="..." z_offset="..." x_offset="..." y_offset="...">  <!--  dimension of the local chamber envelope. x: the half length of the thickness of the chamber, y&z: the half length of the 2D plane dimensions of the chambers-->
  <sensitive type="tracker"/>

  <!-- Specify the detector parameters and the overlap // if you want exclude any component, e.g: endcap, just put endcapDetectorParameters =0 // radius is put in the middle, so its not the inner neither the outer  -->
  <generalParameters numSides="..." overlapY="..." overlapZ="..." clearance="..."/>
  <barrelDetectorParameters numBarrelDetectorLayers ="... detectorVolumeThickness="... radius="... barrelLength="...  />
  <barrelRadiatorParameters numBarrelRadiators="... barrelRadiatorThickness="... barrelRadiatorLayerRadius="... material=".../> 
  <endcapRadiatorParameters numEndcapRadiators="... endcapRadiatorThickness="...   endcapRadiatorLayerInnerRadius="... endcapRadiatorLayerOuterRadius="... endcapRadiatorZOffset="... material=".../>
  <endcapDetectorParameters numEndcapDetectorLayers="... endcapDetectorLayerInnerRadius="... endcapDetectorLayerOuterRadius="... endcapDetectorZOffset="... />


  <layer x="..." y="..." z="..." z_offset="..." x_offset="..." y_offset="..." material="...">
  . . . .
  <layer x="..." y="..." z="..." z_offset="..." x_offset="..." y_offset="..." material="..." sensitive="true">
</detector>

If used with sensitive layers, the readout must contain a "gasLayer" field
*/


#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/XMLElements.h"
#include <sstream> 
#include <iostream>
#include <cmath>

//namespace det {
using namespace std;
using namespace dd4hep;

static dd4hep::Ref_t createmuonSystemMuRWELL_o1_v01(dd4hep::Detector& lcdd,
                                                           dd4hep::xml::Handle_t xmlElement,
                                                           dd4hep::SensitiveDetector sensDet) {
  dd4hep::xml::DetElement xmlDet = static_cast<dd4hep::xml::DetElement>(xmlElement);
  std::string name = xmlDet.nameStr();
  dd4hep::DetElement detElement(name, xmlDet.id());
  dd4hep::Volume experimentalHall = lcdd.pickMotherVolume(detElement);
  

  xml_comp_t dimensions(xmlDet.dimensions());
  // xml_comp_t envelopes(xmlDet.envelopes());
  dd4hep::xml::Dimension envelopeDimensions(xmlDet.envelope());

  //                         --- General parameters ---

  auto generalParameters = xmlElement.child(_Unicode(generalParameters));
  int numSides = generalParameters.attr<int>("numSides");
  double overlapY = generalParameters.attr<double>("overlapY");
  double overlapZ = generalParameters.attr<double>("overlapZ");
  double clearance = generalParameters.attr<double>("clearance"); // it's a small distance to be used to avoid overlapping between the different volumes ~ 1 mm

  //                            --------

  dd4hep::PolyhedraRegular  detectorEnvelope(numSides, envelopeDimensions.rmin(), envelopeDimensions.rmax(), envelopeDimensions.z());
  // dd4hep::Box detectorEnvelope(envelopeDimensions.x() , envelopeDimensions.y() , envelopeDimensions.z());
  dd4hep::Material mat = lcdd.material("Air");
  dd4hep::Volume detectorVolume(name+"Envelope", detectorEnvelope, mat);

  // ----------------------------------------------------------------------------------------------------

  //                         --- Barrel parameters ---

  auto barrelDetectorParameters = xmlElement.child(_Unicode(barrelDetectorParameters));
  int numBarrelDetectorLayers = barrelDetectorParameters.attr<double>("numBarrelDetectorLayers");
  double detectorVolumeThickness = barrelDetectorParameters.attr<double>("detectorVolumeThickness");
  double radius = barrelDetectorParameters.attr<double>("radius");  
  double barrelLength = barrelDetectorParameters.attr<double>("barrelLength");


  auto barrelRadiatorParameters = xmlElement.child(_Unicode(barrelRadiatorParameters));
  int numBarrelRadiators = barrelRadiatorParameters.attr<double>("numBarrelRadiators");
  double barrelRadiatorThickness = barrelRadiatorParameters.attr<double>("barrelRadiatorThickness");
  double barrelRadiatorLayerRadius = barrelRadiatorParameters.attr<double>("barrelRadiatorLayerRadius");
  dd4hep::Material barrelRadiatorMaterial = lcdd.material(barrelRadiatorParameters.attr<std::string>("material"));
  // double barrelRadiatorLayer2Radius = barrelDetectorParameters.attr<double>("barrelRadiatorLayer2Radius");

  //                    --- Endcap parameters ---

  auto endcapDetectorParameters = xmlElement.child(_Unicode(endcapDetectorParameters));
  int numEndcapDetectorLayers = endcapDetectorParameters.attr<double>("numEndcapDetectorLayers");
  double endcapDetectorVolumeThickness = endcapDetectorParameters.attr<double>("endcapDetectorVolumeThickness");
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

  double shapeAngle = 360.0 / numSides; // it's the full angle
  double shapeAngle_radians = M_PI / numSides;  // Isn't it a half angle!!
  double angle_clearance = 0.0; // 0.02 works good but needs the detector volume thick to be more than 60 mm. // it's less than 1 degree, we use the clearnce to avoid overlapping

  double sideLength = 2 * radius * std::tan(shapeAngle_radians); 
  double sideLength2 = 2 * (radius + 1.5 * detectorVolumeThickness) * std::tan(shapeAngle_radians);

  //  double sideTRDLength1 = 2 * (radius - detectorVolumeThickness/2.0) * std::tan(shapeAngle_radians);
  //  double sideTRDLength2 = 2 * (radius + detectorVolumeThickness/2.0) * std::tan(shapeAngle_radians);

  double barrelRadiatorSideLength = 2 * (barrelRadiatorLayerRadius - barrelRadiatorThickness/2.0) * std::tan(shapeAngle_radians);
  double barrelRadiatorSideLength2 = 2 * (barrelRadiatorLayerRadius + barrelRadiatorThickness/2.0) * std::tan(shapeAngle_radians);

  //  double endcapDetectorSideVirtualLength = 2 * (endcapDetectorLayerInnerRadius) * std::tan(shapeAngle_radians); 
  // double endcapDetectorSideVirtualLength2 = 2 * (endcapDetectorLayerOuterRadius) * std::tan(shapeAngle_radians);

  double endcapDetectorSideLength = (2 * (endcapDetectorLayerInnerRadius) * std::tan(shapeAngle_radians)) + 42.0; // shouldn't be hardcoded, but // the offset of 42.0 cm to avoid the overalp with the rectangles
  //double endcapDetectorSideLength2 = (2 * (endcapDetectorLayerOuterRadius) * std::tan(shapeAngle_radians)) + 42.0;  // the offset of 42.0 cm to avoid the overalp with the rectangles

  double endcapDetectorSideTrapLength = (2 * (endcapDetectorLayerOuterRadius - 2 * dimensions.z()) * std::tan(shapeAngle_radians)) + 42.0;  // the offset of 42.0 cm to avoid the overalp with the rectangles
  double endcapDetectorSideTrapYLength = endcapDetectorLayerOuterRadius - 2 * dimensions.z() - endcapDetectorLayerInnerRadius;
  double endcapDetectorSideBoxLength = 2 * (endcapDetectorLayerOuterRadius) * std::tan(shapeAngle_radians);
  double endcapDetectorSideBoxYLength = 2 * dimensions.z();

  double endcapRadiatorSideLength = 2 * (endcapRadiatorLayerInnerRadius) * std::tan(shapeAngle_radians); //it is also the same for the endcap detector layers.
  double endcapRadiatorSideLength2 = 2 * (endcapRadiatorLayerOuterRadius) * std::tan(shapeAngle_radians); //it is also the same for the endcap detector layers.

  double endcapDetectorYLength = endcapDetectorLayerOuterRadius - endcapDetectorLayerInnerRadius;
  double endcapYLength = endcapRadiatorLayerOuterRadius - endcapRadiatorLayerInnerRadius; // It is the distance betwwen the inner and the outer radius of the endcap, it can be in both Y and X dimensions //it is also the same for the endcap detector layers.

  //  int numCopiesY = sideLength / (2 * dimensions.y() - overlapY) + 1;
  //  int numCopiesZ = barrelLength / (2 * dimensions.z() - overlapZ) + 1;

  //  double remainderY = std::fmod((sideLength - 2 * clearance), (2 * dimensions.y() - overlapY)) / (2 * dimensions.y());
  double remainderZ = std::fmod((barrelLength - 2 * clearance), (2 * dimensions.z() - overlapZ)) / (2 * dimensions.z()) - (2 * clearance / dimensions.z());
  double endcapRemainderZ = std::fmod((endcapDetectorYLength - 2 * clearance), (2 * dimensions.z() - overlapZ)) / (2 * dimensions.z()) - (2 * clearance / dimensions.z());
  //  double remainderY2 = std::fmod((sideLength2 - 2 * clearance), (2 * dimensions.y() - overlapY)) / (2 * dimensions.y());

  //  double remainderYLength = remainderY * 2 * dimensions.y();
  //  double remainderZLength = remainderZ * 2 * dimensions.z();

  //  double totalOverlapValueY = (numCopiesY - 1) * overlapY;
  //  double totalOverlapValueZ = (numCopiesZ - 1) * overlapZ;

  double sideEnvX = detectorVolumeThickness/2.0; // 
  double sideEnvY = (sideLength / 2.0); // + clearance * 15;
  double sideEnvZ = (barrelLength / 2.0); // + clearance/2.0;
  double sideEnvY2 = (sideLength2 / 2.0);

  //  double barrelRadiatorEnvX = barrelRadiatorThickness/2.0;
  //  double barrelRadiatorEnvY = barrelRadiatorSideLength2 / 2.0;
  //  double barrelRadiatorEnvZ = barrelLength / 2.0;

  //  double endcapRadiatorEnvX = endcapRadiatorLayerOuterRadius;  // it depends on outer radius
  //  double endcapRadiatorEnvY = endcapRadiatorLayerOuterRadius;  // outer radius too
  //  double endcapRadiatorEnvZ = endcapRadiatorThickness / 2.0;  // layer thickness

  //  double endcapDetectorEnvX = endcapDetectorLayerOuterRadius;  // it depends on outer radius
  //  double endcapDetectorEnvY = endcapDetectorLayerOuterRadius;  // outer radius too
  double endcapDetectorEnvZ = detectorVolumeThickness;  // layer thickness

  //  double sideY = (numCopiesY -1) * 2 * dimensions.y() + remainderYLength - totalOverlapValueY;
  //  double sideZ = (numCopiesZ -1) * 2 * dimensions.z() + remainderZLength - totalOverlapValueZ;

  //  int nameCounter = 0;
  int barrelIdCounter = 1;
  int endcapIdCounter = 1;

// ----------------------------------------------------------------------------------------------------
// ------------------------------// B A R R E L // ----------------------------------------------------


  if (numBarrelDetectorLayers > 0){

  int numBarrelLayer = 1;
  double barrelLayerRMin = radius - detectorVolumeThickness/2.0;   // have to be changed depends on xml impementation
  double barrelLayerRMax = radius + 1.5 * detectorVolumeThickness;

  dd4hep::PolyhedraRegular  BarrelDetectorLayer(numSides, barrelLayerRMin, barrelLayerRMax, barrelLength);
  dd4hep::Volume BarrelDetectorLayerVolume(" barrelDetectorLayer" + numBarrelLayer , BarrelDetectorLayer, mat);


  for (int side = 0; side < numSides; ++side) {

      //dd4hep::Trapezoid sideEnvelope(barrelLength/2.0, barrelLength/2.0, sideTRDLength1/2.0, sideTRDLength2/2.0, detectorVolumeThickness/2.0);      
      dd4hep::Box sideEnvelope(sideEnvX, sideEnvY, sideEnvZ);
      dd4hep::Box sideEnvelope2(sideEnvX, sideEnvY2, sideEnvZ);      
      std::string sideName = dd4hep::xml::_toString(side, "side%d");
      dd4hep::Volume sideVol(sideName, sideEnvelope, mat);
      dd4hep::Volume sideVol2(sideName, sideEnvelope2, mat);      

      double angle_degrees = shapeAngle * side; // Calculate the angle for each side
      double angle_radians = angle_degrees * M_PI / 180.0;

      double sideXOffset = radius * std::cos(angle_radians+shapeAngle_radians);
      double sideYOffset = radius * std::sin(angle_radians+shapeAngle_radians);

      double sideXOffset2 = (radius + detectorVolumeThickness) * std::cos(angle_radians+shapeAngle_radians);
      double sideYOffset2 = (radius + detectorVolumeThickness) * std::sin(angle_radians+shapeAngle_radians);      

      dd4hep::RotationZ sideRotationZ(angle_radians+shapeAngle_radians+angle_clearance);
      //dd4hep::RotationY sideRotationY(90.* dd4hep::degree);     
      dd4hep::Rotation3D sideRotation = dd4hep::Rotation3D(sideRotationZ);      

      double sideXPos = sideXOffset ; 
      double sideYPos = sideYOffset ;   
      double sideZPos = 0.0; 
      double sideXPos2 = sideXOffset2 ; 
      double sideYPos2 = sideYOffset2 ;       

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

      sideDE = dd4hep::DetElement(detElement, sideName + "DE", side);
      sideDE.setPlacement(sidePhys);
      sideVol.setVisAttributes(lcdd, xmlDet.visStr());
      sideVol2.setVisAttributes(lcdd, xmlDet.visStr());

      // -------------------------------------------------------------------------------------------------------
      //  Dividing every side to small rectangles  //
      // -------------------------------------------------------------------------------------------------------

      int numRectangles = barrelLength / (2 * dimensions.z() - overlapZ); // numbers of the rectangles in each side.. depends on the number of the chambers that can be overlapped in Z direction
        
      for (int rectangle = 0; rectangle < (numRectangles + 1); ++rectangle) {

          double rectangleEnvX = detectorVolumeThickness/4.5;  // dividing by 4.5 gets the best thickness for the recatngle to avoid any overlap ~ in our case the uRWELL the best rectangle thickness is 26.667 mm which is 120/4.5.
          double rectangleEnvY;
          if (side % 2 == 0) {
            rectangleEnvY = sideLength / 2.0; // + clearance;
          } else {
            rectangleEnvY = sideLength2 / 2.0; // + clearance;
          } 
          double rectangleEnvZ = dimensions.z() + clearance/2.0;
          double rectangleRemainderEnvZ = (remainderZ * dimensions.z()) + clearance/4.0; // this is the last rectangle in the side, and it is smaller than the usual one, with larger rotation angle. so it need to be shorter in length to avoid overlap.          

          dd4hep::Box rectangleEnvelope(rectangleEnvX, rectangleEnvY, rectangleEnvZ);
          std::string rectangleEnvelopeName = dd4hep::xml::_toString(rectangle, "rectangleEnvelope%d");
          dd4hep::Volume rectangleEnvVol(rectangleEnvelopeName, rectangleEnvelope, mat);

          dd4hep::Box rectangleRemainderEnvelope(rectangleEnvX, rectangleEnvY, rectangleRemainderEnvZ);
          std::string rectangleRemainderEnvelopeName = dd4hep::xml::_toString(rectangle, "rectangleRemainderEnvelope%d");
          dd4hep::Volume rectangleRemainderEnvVol(rectangleRemainderEnvelopeName, rectangleRemainderEnvelope, mat);          

          double rectangleEnvXPos = 0.0; 
          double rectangleEnvYPos = 0.0;   
          double rectangleEnvZPos = barrelLength/2.0 - dimensions.z() - (rectangle * (2 * dimensions.z() - overlapZ)) - clearance;
          double rectangleRemainderEnvXPos = sideEnvX/2.0;
          double rectangleRemainderEnvZPos = barrelLength/2.0 - (remainderZ * dimensions.z()) - (rectangle * (2 * dimensions.z() - overlapZ)) - clearance;

          double yRotation = std::atan(rectangleEnvX / (rectangleEnvZ - (2 * overlapZ)));
          double yRemainderRotation = std::atan(rectangleEnvX / (rectangleRemainderEnvZ - (2 * overlapZ)));
          dd4hep::RotationY rotationY(yRotation);
          dd4hep::RotationY remainderRotationY;
          if (rectangleRemainderEnvZ <= rectangleEnvZ/10.0) {
            remainderRotationY = dd4hep::RotationY(yRotation);
          } else {
            remainderRotationY = dd4hep::RotationY(yRemainderRotation);
          }  

          // _________________________________________________________________________________________________
          if (rectangle == numRectangles) { 
    
            dd4hep::Position rectangeEnvelopeTrans;
            if (rectangleRemainderEnvZ <= rectangleEnvZ/10.0) {
              rectangeEnvelopeTrans = dd4hep::Position(rectangleRemainderEnvXPos, rectangleEnvYPos, rectangleRemainderEnvZPos);
            } else {
              rectangeEnvelopeTrans = dd4hep::Position(rectangleEnvXPos, rectangleEnvYPos, rectangleRemainderEnvZPos);
            }
            dd4hep::PlacedVolume rectangleEnvelopePhys;
            if (side % 2 == 0) {
              rectangleEnvelopePhys = sideVol.placeVolume(rectangleRemainderEnvVol, dd4hep::Transform3D(remainderRotationY, rectangeEnvelopeTrans));
            } else {
             rectangleEnvelopePhys = sideVol2.placeVolume(rectangleRemainderEnvVol, dd4hep::Transform3D(remainderRotationY, rectangeEnvelopeTrans)); 
            }
            dd4hep::DetElement rectangleEnvelopeDE(sideDE, rectangleRemainderEnvelopeName + "DE", rectangle); 
            rectangleEnvelopeDE.setPlacement(rectangleEnvelopePhys);
            rectangleEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

            // ------------------------ start to build the chamber envelopes -------------------

            int numChambersInRectangle = 2 * rectangleEnvY / (2 * dimensions.y() - overlapY); // number of the chambers in each rectangle

            for (int chamberIndex = 0; chamberIndex < (numChambersInRectangle + 1); chamberIndex++) {

            //  int barrelChamberID = side * 1000 + rectangle * 100 + chamberIndex;  // later you should add layer number to distinguish between different barrel layers.
            
            std::stringstream barrelNameStream;
            //   barrelNameStream << "MuRWELL_Barrel_" << barrelChamberID++;
            barrelNameStream << "MuRWELL_Barrel_" << barrelIdCounter++;            
            std::string BarrelChamberName = barrelNameStream.str();

            dd4hep::Box envelope(dimensions.x(), dimensions.y(), remainderZ *  dimensions.z());
            dd4hep::Volume envVolume(BarrelChamberName, envelope, lcdd.material(dimensions.materialStr())); 

            double rectangleRemainderY = std::fmod(2 * (rectangleEnvY - clearance), (2 * dimensions.y() - overlapY)) / (2 * dimensions.y());

            dd4hep::Box rectangleRemainderYEnvelope(dimensions.x(), rectangleRemainderY * dimensions.y(), remainderZ *  dimensions.z());
            dd4hep::Volume rectangleRemainderYEnvVolume(BarrelChamberName + "rectangleRemainderY", rectangleRemainderYEnvelope, lcdd.material(dimensions.materialStr()));          

            double envYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() - rectangleEnvY + dimensions.y() + clearance/20.0 ; // found that the positioning of the chambers inside the rectangle had an overlap with the mother volume ~ 45 um.
            double rectangleRemainderREnvYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() - rectangleEnvY + rectangleRemainderY * dimensions.y() + clearance/20.0;
          
            double zRotation = std::atan(dimensions.x() / (dimensions.y() - (2 * overlapY)));
            dd4hep::RotationZ chamberRotation(zRotation);

            double rectangleRemainderZRotation = std::atan(dimensions.x() / (rectangleRemainderY * dimensions.y() - (2 * overlapY))); // Y and Z are reversed in local remainder
            dd4hep::RotationZ rectangleRemainderRotationZ(rectangleRemainderZRotation);

            auto layers = xmlElement.children(_Unicode(layer));
            auto numLayers = xmlElement.numChildren(_Unicode(layer), true);
            dd4hep::xml::Handle_t layer(layers.reset());
            int sensitiveLayerIndex = 0;

            //std::stringstream nameStream;
            //nameStream << "envDE_" << barrelIdCounter;
            //std::string barrelChamberName = nameStream.str();

            // --- two cases : one for remainder chambers ----------------------------------------------

            if (chamberIndex == numChambersInRectangle) {

              dd4hep::Position rectangleRemainderTrans(dimensions.x_offset(), rectangleRemainderREnvYPos, 0.0);
              dd4hep::PlacedVolume rectangleRemainderEnvPhys = rectangleRemainderEnvVol.placeVolume(rectangleRemainderYEnvVolume, dd4hep::Transform3D(rectangleRemainderRotationZ, rectangleRemainderTrans));
              rectangleRemainderEnvPhys.addPhysVolID("chamber", barrelIdCounter);
              dd4hep::DetElement rectangleRemainderEnvDE(rectangleEnvelopeDE, BarrelChamberName, barrelIdCounter);
              rectangleRemainderEnvDE.setPlacement(rectangleRemainderEnvPhys);
              rectangleRemainderYEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());

              // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() << std::endl;

              for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
              dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
              dd4hep::Box layerShape(layerDet.x(), rectangleRemainderY * layerDet.y(), remainderZ * layerDet.z());
              std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
              dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
              dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
              dd4hep::PlacedVolume layerPlacedVolume = rectangleRemainderYEnvVolume.placeVolume(layerVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), transLayer));

              if (layer.hasAttr("vis")) {
                layerVolume.setVisAttributes(lcdd, layerDet.visStr());
              }
              if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
                dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                sensDet.setType(sdType.typeStr());
                layerVolume.setSensitiveDetector(sensDet);
                layerPlacedVolume.addPhysVolID("gasLayer", sensitiveLayerIndex);
              //  dd4hep::printout(dd4hep::INFO,"Sensitive layer has been created at", name, BarrelChamberName);
                sensitiveLayerIndex++;
              }
              layer.m_node = layers.next();
              }

              // ---------------- Second case: for the full chambers
            } else {

              dd4hep::Position trans(dimensions.x_offset(), envYPos, 0.0);
              dd4hep::PlacedVolume envPhys = rectangleRemainderEnvVol.placeVolume(envVolume, dd4hep::Transform3D(chamberRotation, trans));
              envPhys.addPhysVolID("chamber", barrelIdCounter);
              dd4hep::DetElement envDE(rectangleEnvelopeDE, BarrelChamberName, barrelIdCounter);
              envDE.setPlacement(envPhys);
              envVolume.setVisAttributes(lcdd, xmlDet.visStr());

              // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() << std::endl;

              for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
                  dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
                  dd4hep::Box layerShape(layerDet.x(), layerDet.y(), remainderZ * layerDet.z());
                  std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
                  dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
                  dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
                  dd4hep::PlacedVolume layerPlacedVolume = envVolume.placeVolume(layerVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), transLayer));


                  if (layer.hasAttr("vis")) {
                    layerVolume.setVisAttributes(lcdd, layerDet.visStr());
                  }
                  if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
                    dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                    sensDet.setType(sdType.typeStr());
                    layerVolume.setSensitiveDetector(sensDet);
                    layerPlacedVolume.addPhysVolID("gasLayer", sensitiveLayerIndex);
                    // dd4hep::printout(dd4hep::INFO,"Sensitive layer has been created at", name, BarrelChamberName);
                    sensitiveLayerIndex++;
                  }
                  layer.m_node = layers.next();
              }
            }
          }

          } else {
          // ......_____________________________________________________________________________________________

          dd4hep::Position rectangeEnvelopeTrans(rectangleEnvXPos, rectangleEnvYPos, rectangleEnvZPos);
          dd4hep::PlacedVolume rectangleEnvelopePhys;
          if (side % 2 == 0) {
            rectangleEnvelopePhys = sideVol.placeVolume(rectangleEnvVol, dd4hep::Transform3D(rotationY, rectangeEnvelopeTrans));
          } else {
            rectangleEnvelopePhys = sideVol2.placeVolume(rectangleEnvVol, dd4hep::Transform3D(rotationY, rectangeEnvelopeTrans));
          }
          dd4hep::DetElement rectangleEnvelopeDE(sideDE, rectangleEnvelopeName + "DE", rectangle); 
          rectangleEnvelopeDE.setPlacement(rectangleEnvelopePhys);
          rectangleEnvVol.setVisAttributes(lcdd, xmlDet.visStr());


        //  ------------------------ start to build the chamber envelopes -------------------

          int numChambersInRectangle = 2 * rectangleEnvY / (2 * dimensions.y() - overlapY); // number of the chambers in each rectangle

          for (int chamberIndex = 0; chamberIndex < (numChambersInRectangle + 1); chamberIndex++) {

            // int barrelChamberID = side * 1000 + rectangle * 100 + chamberIndex;  // later you should add layer number to distinguish between different barrel layers.
            
            std::stringstream barrelNameStream;
            //barrelNameStream << "MuRWELL_Barrel_" << barrelChamberID++;
            barrelNameStream << "MuRWELL_Barrel_" << barrelIdCounter++;            
            std::string BarrelChamberName = barrelNameStream.str();

            dd4hep::Box envelope(dimensions.x(), dimensions.y(), dimensions.z());
            dd4hep::Volume envVolume(BarrelChamberName, envelope, lcdd.material(dimensions.materialStr())); 

            double rectangleRemainderY = std::fmod(2 * (rectangleEnvY - clearance), (2 * dimensions.y() - overlapY)) / (2 * dimensions.y());

            dd4hep::Box rectangleRemainderYEnvelope(dimensions.x(), rectangleRemainderY * dimensions.y(), dimensions.z());
            dd4hep::Volume rectangleRemainderYEnvVolume(BarrelChamberName + "rectangleRemainderY", rectangleRemainderYEnvelope, lcdd.material(dimensions.materialStr()));          

            double envYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() - rectangleEnvY + dimensions.y() + clearance/20.0 ; // found that the positioning of the chambers inside the rectangle had an overlap with the mother volume ~ 45 um.
            double rectangleRemainderREnvYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() - rectangleEnvY + rectangleRemainderY * dimensions.y() + clearance/20.0;
          
            double zRotation = std::atan(dimensions.x() / (dimensions.y() - (2 * overlapY)));
            dd4hep::RotationZ chamberRotation(zRotation);

            double rectangleRemainderZRotation = std::atan(dimensions.x() / (rectangleRemainderY * dimensions.y() - (2 * overlapY))); // Y and Z are reversed in local remainder
            dd4hep::RotationZ rectangleRemainderRotationZ(rectangleRemainderZRotation);

            auto layers = xmlElement.children(_Unicode(layer));
            auto numLayers = xmlElement.numChildren(_Unicode(layer), true);
            dd4hep::xml::Handle_t layer(layers.reset());
            int sensitiveLayerIndex = 0;

            // std::stringstream nameStream;
            // nameStream << "envDE_" << barrelIdCounter;
            // std::string barrelChamberName = nameStream.str();

          // --- two cases : one for remainder chambers ----------------------------------------------

            if (chamberIndex == numChambersInRectangle) {

              dd4hep::Position rectangleRemainderTrans(dimensions.x_offset(), rectangleRemainderREnvYPos, 0.0);
              dd4hep::PlacedVolume rectangleRemainderEnvPhys = rectangleEnvVol.placeVolume(rectangleRemainderYEnvVolume, dd4hep::Transform3D(rectangleRemainderRotationZ, rectangleRemainderTrans));
              rectangleRemainderEnvPhys.addPhysVolID("chamber", barrelIdCounter);
              dd4hep::DetElement rectangleRemainderEnvDE(rectangleEnvelopeDE, BarrelChamberName, barrelIdCounter);
              rectangleRemainderEnvDE.setPlacement(rectangleRemainderEnvPhys);
              rectangleRemainderYEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());

              // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() << std::endl;

              for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
              dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
              dd4hep::Box layerShape(layerDet.x(), rectangleRemainderY * layerDet.y(), layerDet.z());
              std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
              dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
              dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
              dd4hep::PlacedVolume layerPlacedVolume = rectangleRemainderYEnvVolume.placeVolume(layerVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), transLayer));

              if (layer.hasAttr("vis")) {
                layerVolume.setVisAttributes(lcdd, layerDet.visStr());
              }
              if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
                dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                sensDet.setType(sdType.typeStr());
                layerVolume.setSensitiveDetector(sensDet);
                layerPlacedVolume.addPhysVolID("gasLayer", sensitiveLayerIndex);
              //  dd4hep::printout(dd4hep::INFO,"Sensitive layer has been created at", name, BarrelChamberName);
                sensitiveLayerIndex++;
              }
              layer.m_node = layers.next();
              }

          // ---------------- Second case: for the full chambers
            } else {

              dd4hep::Position trans(dimensions.x_offset(), envYPos, 0.0);
              dd4hep::PlacedVolume envPhys = rectangleEnvVol.placeVolume(envVolume, dd4hep::Transform3D(chamberRotation, trans));
              envPhys.addPhysVolID("chamber", barrelIdCounter);
              dd4hep::DetElement envDE(rectangleEnvelopeDE, BarrelChamberName, barrelIdCounter);
              envDE.setPlacement(envPhys);
              envVolume.setVisAttributes(lcdd, xmlDet.visStr());

              // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() << std::endl;

              for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
                  dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
                  dd4hep::Box layerShape(layerDet.x(), layerDet.y(), layerDet.z());
                  std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
                  dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
                  dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
                  dd4hep::PlacedVolume layerPlacedVolume = envVolume.placeVolume(layerVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), transLayer));


                  if (layer.hasAttr("vis")) {
                    layerVolume.setVisAttributes(lcdd, layerDet.visStr());
                  }
                  if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
                    dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
                    sensDet.setType(sdType.typeStr());
                    layerVolume.setSensitiveDetector(sensDet);
                    layerPlacedVolume.addPhysVolID("gasLayer", sensitiveLayerIndex);
                    // dd4hep::printout(dd4hep::INFO,"Sensitive layer has been created at", name, BarrelChamberName);
                    sensitiveLayerIndex++;
                  }
                  layer.m_node = layers.next();
              }
            }
          }
        }          
    }
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------

    }

  dd4hep::Position detectorLayerTrans(0., 0., 0.);
  dd4hep::PlacedVolume detectorLayerPhys = detectorVolume.placeVolume(BarrelDetectorLayerVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), detectorLayerTrans));
  detectorLayerPhys.addPhysVolID("layer", numBarrelLayer); 
  dd4hep::DetElement BarrelDetectorLayerDE(detElement, " barrelDetectorLayerDE", numBarrelLayer);
  BarrelDetectorLayerDE.setPlacement(detectorLayerPhys);
  BarrelDetectorLayerVolume.setVisAttributes(lcdd.visAttributes("no_vis")); 
  numBarrelLayer++;

  }
// ---------------------------------------// B A R R E L // R A D I A T O R S //------------------------------------------------

  if ( numBarrelRadiators > 0) {

    int numBarrelRadiatorLayer = 1;
    double barrelRadiatorLayerRMin = barrelRadiatorLayerRadius - barrelRadiatorThickness/2.0;  // have to be changed depends on xml impementation
    double barrelRadiatorLayerRMax = barrelRadiatorLayerRadius + barrelRadiatorThickness/2.0;    

    dd4hep::PolyhedraRegular  BarrelRadiatorLayer(numSides, barrelRadiatorLayerRMin, barrelRadiatorLayerRMax, barrelLength);
    dd4hep::Volume BarrelRadiatorLayerVolume(" barrelRadiatorLayer" + numBarrelRadiatorLayer , BarrelRadiatorLayer, mat);

    for (int side = 0; side < numSides; ++side) {

      dd4hep::Trapezoid barrelRadiatorEnvelope(barrelLength/2.0, barrelLength/2.0, barrelRadiatorSideLength/2.0, barrelRadiatorSideLength2/2.0, barrelRadiatorThickness/2.0);
      //  dd4hep::Box barrelRadiatorEnvelope(barrelRadiatorEnvX, barrelRadiatorEnvY, barrelRadiatorEnvZ);
      std::string barrelRadiatorEnvelopeName = dd4hep::xml::_toString(side, "barrelRadiatorEnvelope%d");
      dd4hep::Volume barrelRadiatorEnvVol(barrelRadiatorEnvelopeName, barrelRadiatorEnvelope, mat);

      double angle_degrees = shapeAngle * side; // Calculate the angle for each side
      double angle_radians = angle_degrees * M_PI / 180.0;

      double barrelRadiatorXOffset = barrelRadiatorLayerRadius * std::cos(angle_radians+shapeAngle_radians);
      double barrelRadiatorYOffset = barrelRadiatorLayerRadius * std::sin(angle_radians+shapeAngle_radians);

      dd4hep::RotationZ barrelRadiatorRotationZ(angle_radians+shapeAngle_radians);
      dd4hep::RotationY barrelRadiatorRotationY(90.* dd4hep::degree);

      dd4hep::Rotation3D barrelRadiatorRotation = barrelRadiatorRotationZ * barrelRadiatorRotationY;

      double barrelRadiatorXPos = barrelRadiatorXOffset ; 
      double barrelRadiatorYPos = barrelRadiatorYOffset ;   
      double barrelRadiatorZPos = 0.0; 

      dd4hep::Position barrelRadiatorEnvelopeTrans(barrelRadiatorXPos, barrelRadiatorYPos, barrelRadiatorZPos);
      dd4hep::PlacedVolume barrelRadiatorEnvelopePhys = BarrelRadiatorLayerVolume.placeVolume(barrelRadiatorEnvVol, dd4hep::Transform3D(barrelRadiatorRotation, barrelRadiatorEnvelopeTrans));
      //  barrelRadiatorEnvelopePhys.addPhysVolID("barrelRadiatorEnvelope", side);
      dd4hep::DetElement barrelRadiatorEnvelopeDE(detElement, barrelRadiatorEnvelopeName + "DE", side);
      barrelRadiatorEnvelopeDE.setPlacement(barrelRadiatorEnvelopePhys);
      barrelRadiatorEnvVol.setVisAttributes(lcdd, xmlDet.visStr()); 

    // ----------------------------------------- Building the yokes ------------------

      dd4hep::Trapezoid barrelRadiator(barrelLength/2.0, barrelLength/2.0, barrelRadiatorSideLength/2.0, barrelRadiatorSideLength2/2.0, barrelRadiatorThickness/2.0);
      std::string barrelRadiatorName = dd4hep::xml::_toString(side, "barrelRadiator%d");
      dd4hep::Volume barrelRadiatorVol(barrelRadiatorName, barrelRadiator, barrelRadiatorMaterial);


      dd4hep::Position barrelRadiatorTrans(0, 0, 0);
      dd4hep::PlacedVolume barrelRadiatorPhys = barrelRadiatorEnvVol.placeVolume(barrelRadiatorVol, dd4hep::Transform3D(dd4hep::RotationY(0.), barrelRadiatorTrans));
      //  barrelRadiatorPhys.addPhysVolID("barrelRadiator", side);
      dd4hep::DetElement barrelRadiatorDE(barrelRadiatorEnvelopeDE, barrelRadiatorName + "DE", side);
      barrelRadiatorDE.setPlacement(barrelRadiatorPhys);
      barrelRadiatorVol.setVisAttributes(lcdd.visAttributes("yoke_vis"));       

    }  

    dd4hep::Position radiatorLayerTrans(0., 0., 0.);
    dd4hep::PlacedVolume radiatorLayerPhys = detectorVolume.placeVolume(BarrelRadiatorLayerVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), radiatorLayerTrans));
    dd4hep::DetElement BarrelRadiatorLayerDE(detElement, " Barrel_RadiatorLayerDE", numBarrelRadiatorLayer);
    BarrelRadiatorLayerDE.setPlacement(radiatorLayerPhys);
    BarrelRadiatorLayerVolume.setVisAttributes(lcdd.visAttributes("no_vis")); 
    numBarrelRadiatorLayer++;

  }

// -------------------------------------------------------------------------------------------
//  ----------------------------------// E N D C A P //---------------------------------------
//--------------------------------------- Endcap Detectors------------------------------------

  if ( numEndcapDetectorLayers > 0) {

    dd4hep::PolyhedraRegular  endcapDetectorEnvelope(numSides, endcapDetectorLayerInnerRadius, endcapDetectorLayerOuterRadius, endcapDetectorVolumeThickness);
    //dd4hep::Box endcapDetectorEnvelope(endcapDetectorEnvX, endcapDetectorEnvY, endcapDetectorEnvZ);
    std::string endcapDetectorEnvelopeName = dd4hep::xml::_toString(numEndcapDetectorLayers, "endcapDetectorEnvelope%d");
    dd4hep::Volume endcapDetectorEnvVol(endcapDetectorEnvelopeName, endcapDetectorEnvelope, mat);

    double endcapDetectorEnvXPos = 0.0 ; 
    double endcapDetectorEnvYPos = 0.0 ;   
    double endcapDetectorEnvZPos = endcapDetectorZOffset; 

    dd4hep::Position endcapDetectorEnvelopeTrans(endcapDetectorEnvXPos, endcapDetectorEnvYPos, endcapDetectorEnvZPos);
    dd4hep::PlacedVolume endcapDetectorEnvelopePhys = detectorVolume.placeVolume(endcapDetectorEnvVol, dd4hep::Transform3D(dd4hep::RotationZ(0.), endcapDetectorEnvelopeTrans));
    dd4hep::DetElement endcapDetectorEnvelopeDE(detElement, endcapDetectorEnvelopeName + "DE", numEndcapDetectorLayers); // remember to loop over numEndcapDetectorLayers.. because now it still just one number.
    endcapDetectorEnvelopeDE.setPlacement(endcapDetectorEnvelopePhys);
    endcapDetectorEnvVol.setVisAttributes(lcdd, xmlDet.visStr()); 

    // --------------------------Building detector sides ---------------------

    for (int side = 0; side < numSides; ++side) {

      dd4hep::Trapezoid endcapDetectorSideTrap(endcapDetectorEnvZ/2.0, endcapDetectorEnvZ/2.0, endcapDetectorSideLength/2.0, endcapDetectorSideTrapLength/2.0, endcapDetectorSideTrapYLength/2.0);
      dd4hep::Box endcapDetectorSideBox(endcapDetectorEnvZ/2.0, endcapDetectorSideBoxLength/2.0, endcapDetectorSideBoxYLength/2.0);

      double boxOffsetZ = endcapDetectorYLength/2.0;

      dd4hep::Position boxPos(0.0 , 0.0, boxOffsetZ);
      dd4hep::Rotation3D boxRot(dd4hep::RotationY(0.0 * dd4hep::degree));
      
      dd4hep::Transform3D boxTransform(boxRot, boxPos);

      //Combining two shapes by UnionSolid: the first shape is centralized and the second transform around the first..
      dd4hep::UnionSolid endcapDetectorSideEnvelope(endcapDetectorSideTrap, endcapDetectorSideBox, boxTransform);
      std::string endcapDetectorSideEnvName = dd4hep::xml::_toString(side, "endcapDetectorSideEnv%d");
      dd4hep::Volume endcapDetectorSideEnvVol(endcapDetectorSideEnvName, endcapDetectorSideEnvelope, mat);

      double angle_degrees = shapeAngle * side; // Calculate the angle for each side
      double angle_radians = angle_degrees * M_PI / 180.0;

      //double endcapDetectorMidRadius = endcapDetectorLayerInnerRadius + (endcapDetectorYLength /2.0);
      double endcapDetectorTrapCenterRadius = endcapDetectorLayerInnerRadius + (endcapDetectorSideTrapYLength/2.0);      

      double endcapDetectorXOffset = endcapDetectorTrapCenterRadius * std::cos(angle_radians+shapeAngle_radians);
      double endcapDetectorYOffset = endcapDetectorTrapCenterRadius * std::sin(angle_radians+shapeAngle_radians);

      dd4hep::RotationZ endcapDetectorRotationZ(angle_radians+shapeAngle_radians);
      dd4hep::Rotation3D endcapDetectorRotation = dd4hep::Rotation3D(endcapDetectorRotationZ);

      double endcapDetectorSideEnvXPos = endcapDetectorXOffset ; 
      double endcapDetectorSideEnvYPos = endcapDetectorYOffset ;   
      double endcapDetectorSideEnvZPos = - endcapDetectorEnvZ/2.0 ;   
      double endcapDetectorSideEnvZPos2 = endcapDetectorEnvZ/2.0;

     // ---------  here I start to divide the two z-positions // by odd and even numbers

      dd4hep::Position endcapDetectorSideEnvTrans;
      dd4hep::PlacedVolume endcapDetectorSideEnvPhys;
      dd4hep::DetElement endcapDetectorSideEnvDE;

      if (side % 2 == 0) {

        endcapDetectorSideEnvTrans = dd4hep::Position(endcapDetectorSideEnvXPos, endcapDetectorSideEnvYPos, endcapDetectorSideEnvZPos);
        endcapDetectorSideEnvPhys = endcapDetectorEnvVol.placeVolume(endcapDetectorSideEnvVol, dd4hep::Transform3D(endcapDetectorRotation * dd4hep::RotationY(90.0 * dd4hep::degree) , endcapDetectorSideEnvTrans));
        endcapDetectorSideEnvDE = dd4hep::DetElement(endcapDetectorEnvelopeDE, endcapDetectorSideEnvName + "DE", side);
        endcapDetectorSideEnvDE.setPlacement(endcapDetectorSideEnvPhys);
        endcapDetectorSideEnvVol.setVisAttributes(lcdd, xmlDet.visStr());

      } else {

        endcapDetectorSideEnvTrans = dd4hep::Position(endcapDetectorSideEnvXPos, endcapDetectorSideEnvYPos, endcapDetectorSideEnvZPos2);
        endcapDetectorSideEnvPhys = endcapDetectorEnvVol.placeVolume(endcapDetectorSideEnvVol, dd4hep::Transform3D(endcapDetectorRotation * dd4hep::RotationY(90.0 * dd4hep::degree) , endcapDetectorSideEnvTrans));
        endcapDetectorSideEnvDE = dd4hep::DetElement(endcapDetectorEnvelopeDE, endcapDetectorSideEnvName + "DE", side);
        endcapDetectorSideEnvDE.setPlacement(endcapDetectorSideEnvPhys);
        endcapDetectorSideEnvVol.setVisAttributes(lcdd, xmlDet.visStr());
      }

     // ----- dividing the trapezoid envelope to smaller pieces (rectangles)

      int numRectangles = endcapDetectorYLength / (2 * dimensions.z() - overlapZ); // numbers of the rectangles in each teapezoids.. depends on the number of the chambers that can be overlapped in Y direction
        
      for (int rectangle = 0; rectangle < (numRectangles + 1); ++rectangle) {

          double rectangleEnvX = endcapDetectorEnvZ/4.5;  // dividing by 4.5 gets the best thickness for the recatngle to avoid any overlap ~ in our case the uRWELL the thickness is 26.667 mm which is 120/4.5.
          double rectangleEnvY = (endcapDetectorLayerOuterRadius - rectangle * (2 * dimensions.y() - overlapY)) * std::tan(shapeAngle_radians) ; // without multiplying by 2 .. because it is the half length // it should be dimensions.x() instead of z, but in the endcap its perpendicular to usual dimension set
          double rectangleEnvZ;
          if (rectangle == numRectangles) {
            rectangleEnvZ =  (endcapRemainderZ * dimensions.z()) + clearance/4.0;
          } else {
            rectangleEnvZ = dimensions.z() + clearance/2.0;
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
            rectangleEnvZPos = endcapDetectorYLength/2.0 + ((1-endcapRemainderZ) * dimensions.z()) - (rectangle * (2 * dimensions.z() - overlapZ)) - clearance;
          } else {
            rectangleEnvZPos = endcapDetectorYLength/2.0 - (rectangle * (2 * dimensions.y() - overlapY)) - clearance;
          }  
           
          double yRotation = std::atan(rectangleEnvX / (rectangleEnvZ - (2 * overlapY)));
          dd4hep::RotationY rotationY(yRotation);

          dd4hep::Position rectangeEnvelopeTrans(rectangleEnvXPos, rectangleEnvYPos, rectangleEnvZPos);
          dd4hep::PlacedVolume rectangleEnvelopePhys = endcapDetectorSideEnvVol.placeVolume(rectangleEnvVol, dd4hep::Transform3D(rotationY, rectangeEnvelopeTrans));
          dd4hep::DetElement rectangleEnvelopeDE(endcapDetectorSideEnvDE, rectangleEnvelopeName + "DE", rectangle); // remember to loop over numEndcapDetectorLayers.. because now it still just one number.
          rectangleEnvelopeDE.setPlacement(rectangleEnvelopePhys);
          rectangleEnvVol.setVisAttributes(lcdd, xmlDet.visStr());
        
          // ------------------------ start to build the chamber envelopes -------------------

          int numChambersInRectangle = 2 * rectangleEnvY / (2 * dimensions.y() - overlapY); // number of the chambers in each rectangle

          for (int chamberIndex = 0; chamberIndex < (numChambersInRectangle + 1); chamberIndex++) {

          //  int endcapChamberID = side * 1000 + rectangle * 100 + chamberIndex;  // later you should add layer number to distinguish between different endcap layers.
            
          std::stringstream endcapNameStream;
          endcapNameStream << "MuRWELL_Endcap_" << endcapIdCounter++;
          std::string EndcapChamberName = endcapNameStream.str();

          dd4hep::Box envelope;
          if (rectangle == numRectangles) {
            envelope = dd4hep::Box(dimensions.x(), dimensions.y(),endcapRemainderZ * dimensions.z());
          } else {
            envelope = dd4hep::Box(dimensions.x(), dimensions.y(), dimensions.z());
          } 

          dd4hep::Volume envVolume(EndcapChamberName, envelope, lcdd.material(dimensions.materialStr())); 

          double rectangleRemainderY = std::fmod(2 * (rectangleEnvY - clearance), (2 * dimensions.y() - overlapY)) / (2 * dimensions.y());
          //double rectangleRemainderYLength = rectangleRemainderY * 2 * dimensions.y();

          dd4hep::Box rectangleRemainderYEnvelope;
          if (rectangle == numRectangles) {
            rectangleRemainderYEnvelope = dd4hep::Box(dimensions.x(), rectangleRemainderY * dimensions.y(), endcapRemainderZ * dimensions.z());
          } else {
            rectangleRemainderYEnvelope = dd4hep::Box(dimensions.x(), rectangleRemainderY * dimensions.y(), dimensions.z());
          } 

          dd4hep::Volume rectangleRemainderYEnvVolume(EndcapChamberName + "rectangleRemainderY", rectangleRemainderYEnvelope, lcdd.material(dimensions.materialStr()));          

          double envYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() - rectangleEnvY + dimensions.y() + 0.005 ; // found that the positioning of the chambers inside the rectangle had an overlap with the mother volume ~ 45 um.
          double rectangleRemainderREnvYPos = (chamberIndex * 2 * dimensions.y()) - (overlapY * chamberIndex) + dimensions.y_offset() - rectangleEnvY + rectangleRemainderY * dimensions.y() + 0.005;
          
          double zRotation = std::atan(dimensions.x() / (dimensions.z() - (2 * overlapZ)));
          dd4hep::RotationZ chamberRotation(zRotation);

          double rectangleRemainderZRotation = std::atan(dimensions.x() / (rectangleRemainderY * dimensions.z() - (2 * overlapZ))); // Y and Z are reversed in local remainder
          dd4hep::RotationZ rectangleRemainderRotationZ(rectangleRemainderZRotation);

          auto layers = xmlElement.children(_Unicode(layer));
          auto numLayers = xmlElement.numChildren(_Unicode(layer), true);
          dd4hep::xml::Handle_t layer(layers.reset());
          int sensitiveLayerIndex = 0;

          //  std::stringstream nameStream;
          //  nameStream << "envDE_" << endcapIdCounter;
          //  std::string endcapChamberName = nameStream.str();

          // --- two cases : one for full chambers ----------------------------------------------

          if (chamberIndex == numChambersInRectangle) {

          dd4hep::Position rectangleRemainderTrans(dimensions.x_offset(), rectangleRemainderREnvYPos, 0.0);
          dd4hep::PlacedVolume rectangleRemainderEnvPhys = rectangleEnvVol.placeVolume(rectangleRemainderYEnvVolume, dd4hep::Transform3D(rectangleRemainderRotationZ , rectangleRemainderTrans));
          rectangleRemainderEnvPhys.addPhysVolID("chamber", endcapIdCounter);
          dd4hep::DetElement rectangleRemainderEnvDE(rectangleEnvelopeDE, EndcapChamberName, endcapIdCounter);
          rectangleRemainderEnvDE.setPlacement(rectangleRemainderEnvPhys);
          rectangleRemainderYEnvVolume.setVisAttributes(lcdd, xmlDet.visStr());

          // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() << std::endl;

         for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
          dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
          dd4hep::Box layerShape;
          if (rectangle == numRectangles) {
            layerShape = dd4hep::Box(layerDet.x(), rectangleRemainderY * layerDet.y(), endcapRemainderZ * layerDet.z());
          } else {
            layerShape = dd4hep::Box(layerDet.x(), rectangleRemainderY * layerDet.y(), layerDet.z());
          }  
          std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
          dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
          dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
          dd4hep::PlacedVolume layerPlacedVolume = rectangleRemainderYEnvVolume.placeVolume(layerVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), transLayer));
          //  dd4hep::DetElement layerDE(envDE, "layerDE", layerIdx);
          //  layerDE.setPlacement(layerPlacedVolume);

          if (layer.hasAttr("vis")) {
            layerVolume.setVisAttributes(lcdd, layerDet.visStr());
          }
          if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
            dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
            sensDet.setType(sdType.typeStr());
            layerVolume.setSensitiveDetector(sensDet);
            layerPlacedVolume.addPhysVolID("gasLayer", sensitiveLayerIndex);
          //  dd4hep::printout(dd4hep::INFO,"Sensitive layer has been created at", name,EndcapChamberName);
            sensitiveLayerIndex++;
          }
          layer.m_node = layers.next();
         }

          // ----------------
          } else {

          dd4hep::Position trans(dimensions.x_offset(), envYPos, 0.0);
          dd4hep::PlacedVolume envPhys = rectangleEnvVol.placeVolume(envVolume, dd4hep::Transform3D(chamberRotation , trans));
          envPhys.addPhysVolID("chamber", endcapIdCounter);
          dd4hep::DetElement envDE(rectangleEnvelopeDE, EndcapChamberName, endcapIdCounter);
          envDE.setPlacement(envPhys);
          envVolume.setVisAttributes(lcdd, xmlDet.visStr());

          // std::cout << "Adding detector element: " << detElement.name() << " to path: " << detElement.path() << std::endl;

        for (unsigned layerIdx = 0; layerIdx < numLayers; ++layerIdx) {
          dd4hep::xml::DetElement layerDet = static_cast<dd4hep::xml::DetElement>(layer);
          dd4hep::Box layerShape;
          if (rectangle == numRectangles) {
            layerShape = dd4hep::Box(layerDet.x(), layerDet.y(), endcapRemainderZ * layerDet.z());
          } else {
            layerShape = dd4hep::Box(layerDet.x(), layerDet.y(), layerDet.z());
          }            
          std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
          dd4hep::Volume layerVolume(layerName, layerShape, lcdd.material(layer.attr<std::string>("material")));
          dd4hep::Position transLayer(layerDet.x_offset(), layerDet.y_offset(), layerDet.z_offset());
          dd4hep::PlacedVolume layerPlacedVolume = envVolume.placeVolume(layerVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), transLayer));
          //  dd4hep::DetElement layerDE(envDE, "layerDE", layerIdx);
          //  layerDE.setPlacement(layerPlacedVolume);

          if (layer.hasAttr("vis")) {
            layerVolume.setVisAttributes(lcdd, layerDet.visStr());
          }
          if (layer.hasAttr("sensitive") && layerDet.isSensitive()) {
            dd4hep::xml::Dimension sdType(xmlElement.child(_U(sensitive)));
            sensDet.setType(sdType.typeStr());
            layerVolume.setSensitiveDetector(sensDet);
            layerPlacedVolume.addPhysVolID("gasLayer", sensitiveLayerIndex);
           // dd4hep::printout(dd4hep::INFO,"Sensitive layer has been created at", name, EndcapChamberName);
            sensitiveLayerIndex++;
          }
        layer.m_node = layers.next();
        }
        }
        }          
      }
    }
  }          
// --------------------------------------Radiators--------------------------------------------
  if ( numEndcapRadiators > 0) {

    dd4hep::PolyhedraRegular  endcapRadiatorEnvelope(numSides, endcapRadiatorLayerInnerRadius, endcapRadiatorLayerOuterRadius, endcapRadiatorThickness);
    //dd4hep::Box endcapRadiatorEnvelope(endcapRadiatorEnvX, endcapRadiatorEnvY, endcapRadiatorEnvZ);
    std::string endcapRadiatorEnvelopeName = dd4hep::xml::_toString(numEndcapRadiators, "endcapRadiatorEnvelope%d");
    dd4hep::Volume endcapRadiatorEnvVol(endcapRadiatorEnvelopeName, endcapRadiatorEnvelope, mat);

    double endcapRadiatorEnvXPos = 0.0 ; 
    double endcapRadiatorEnvYPos = 0.0 ;   
    double endcapRadiatorEnvZPos = endcapRadiatorZOffset; 

    dd4hep::Position endcapRadiatorEnvelopeTrans(endcapRadiatorEnvXPos, endcapRadiatorEnvYPos, endcapRadiatorEnvZPos);
    dd4hep::PlacedVolume endcapRadiatorEnvelopePhys = detectorVolume.placeVolume(endcapRadiatorEnvVol, dd4hep::Transform3D(dd4hep::RotationZ(0.), endcapRadiatorEnvelopeTrans));
    dd4hep::DetElement endcapRadiatorEnvelopeDE(detElement, endcapRadiatorEnvelopeName + "DE", numEndcapRadiators);
    endcapRadiatorEnvelopeDE.setPlacement(endcapRadiatorEnvelopePhys);
    endcapRadiatorEnvVol.setVisAttributes(lcdd, xmlDet.visStr()); 

  // ----------------Building radiator sides------------------

    for (int side = 0; side < numSides; ++side) {

      dd4hep::Trapezoid endcapRadiator(endcapRadiatorThickness/2.0, endcapRadiatorThickness/2.0, endcapRadiatorSideLength/2.0, endcapRadiatorSideLength2/2.0, endcapYLength/2.0);
      std::string endcapRadiatorName = dd4hep::xml::_toString(side, "endcapRadiator%d");
      dd4hep::Volume endcapRadiatorVol(endcapRadiatorName, endcapRadiator, endcapRadiatorMaterial);

      double angle_degrees = shapeAngle * side; // Calculate the angle for each side
      double angle_radians = angle_degrees * M_PI / 180.0;  //radian angle for each side (rotating angle)

      double endcapRadiatorMidRadius = endcapRadiatorLayerInnerRadius + (endcapYLength /2.0);

      double endcapRadiatorXOffset = endcapRadiatorMidRadius * std::cos(angle_radians+shapeAngle_radians);
      double endcapRadiatorYOffset = endcapRadiatorMidRadius * std::sin(angle_radians+shapeAngle_radians);

      dd4hep::RotationZ endcapRadiatorRotationZ(angle_radians+shapeAngle_radians);
      dd4hep::Rotation3D endcapRadiatorRotation = dd4hep::Rotation3D(endcapRadiatorRotationZ);

      double endcapRadiatorXPos = endcapRadiatorXOffset ; 
      double endcapRadiatorYPos = endcapRadiatorYOffset ;   
      double endcapRadiatorZPos = 0.0; 

      dd4hep::Position endcapRadiatorTrans(endcapRadiatorXPos, endcapRadiatorYPos, endcapRadiatorZPos);
      dd4hep::PlacedVolume endcapRadiatorPhys = endcapRadiatorEnvVol.placeVolume(endcapRadiatorVol, dd4hep::Transform3D(endcapRadiatorRotation * dd4hep::RotationY(90.* dd4hep::degree), endcapRadiatorTrans));
      dd4hep::DetElement endcapRadiatorDE(endcapRadiatorEnvelopeDE, endcapRadiatorName + "DE", side);
      endcapRadiatorDE.setPlacement(endcapRadiatorPhys);
      endcapRadiatorVol.setVisAttributes(lcdd.visAttributes("yoke_vis"));       
    
    }
  }

  // ------------------------------------------------------------------------------------------- 
  dd4hep::Position detectorTrans(0., 0., envelopeDimensions.z_offset());
  dd4hep::PlacedVolume detectorPhys = experimentalHall.placeVolume(detectorVolume, dd4hep::Transform3D(dd4hep::RotationZ(shapeAngle_radians), detectorTrans));
  detectorPhys.addPhysVolID("system", xmlDet.id());
  detElement.setPlacement(detectorPhys);
  detectorVolume.setVisAttributes(lcdd.visAttributes("no_vis")); 
  return detElement;
}

//}

DECLARE_DETELEMENT(muonSystemMuRWELL_o1_v01, createmuonSystemMuRWELL_o1_v01)
