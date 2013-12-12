//#include "DDBeamCal.hh"


#include <DD4hep/DetFactoryHelper.h>
#include <XML/Layering.h>

#include <string>


static DD4hep::Geometry::Ref_t create_detector(DD4hep::Geometry::LCDD& lcdd,
					       DD4hep::XML::Handle_t xmlHandle,
					       DD4hep::Geometry::SensitiveDetector sens) {

  std::cout << __PRETTY_FUNCTION__  << std::endl;
  std::cout << "Here is my BeamCal"  << std::endl;
  //Materials
  DD4hep::Geometry::Material air = lcdd.air();

  //Access to the XML File
  DD4hep::XML::DetElement xmlBeamCal = xmlHandle;
  const std::string detName = xmlBeamCal.nameStr();

  DD4hep::XML::Dimension dimensions =  xmlBeamCal.dimensions();

  //BeamCal Dimensions
  const double bcalInnerR = dimensions.inner_r();
  const double bcalOuterR = dimensions.outer_r();
  const double bcalInnerZ = dimensions.inner_z();
  const double bcalThickness = DD4hep::Layering(xmlBeamCal).totalThickness();
  const double bcalCentreZ = bcalInnerZ+bcalThickness*0.5;

  // counter for the current layer to be placed
  int thisLayerId = 1;

  //GlobalParameters we have to know about
  const double fullCrossingAngle = 0.020/*rad*/;
  const double incomingBeamPipeRadius = 3.5/*mm*/;


  ////////////////////////////////////////////////////////////////////////////////
  //Calculations for the position of the incoming beampipe
  ////////////////////////////////////////////////////////////////////////////////

  // this needs the full crossing angle, because the beamCal is centred on the
  // outgoing beampipe, and the incoming beampipe is a full crossing agle away
  const DD4hep::Geometry::Position    incomingBeamPipePosition( std::tan(fullCrossingAngle) * bcalCentreZ, 0.0, 0.0);
  //The extra parenthesis are paramount!
  const DD4hep::Geometry::Rotation3D  incomingBeamPipeRotation( ( DD4hep::Geometry::RotationY( fullCrossingAngle ) ) );
  const DD4hep::Geometry::Transform3D incomingBPTransform( incomingBeamPipeRotation, incomingBeamPipePosition );

  //Envelope to place the layers in
  DD4hep::Geometry::Tube envelopeTube (bcalInnerR,bcalOuterR,bcalThickness/2,0,2*M_PI);
  DD4hep::Geometry::Tube incomingBeamPipe (0.0, incomingBeamPipeRadius, bcalThickness);//we want this to be longer than the BeamCal
  DD4hep::Geometry::SubtractionSolid envelope (envelopeTube, incomingBeamPipe, incomingBPTransform);
  DD4hep::Geometry::Volume     envelopeVol(detName+"_envelope",envelope,air);



  //Create the section cutout for the sensor and readout volumes
#warning "These need to be picked up from XML"
  const double bcalCutOutStart =-M_PI *  20/180;
  const double bcalCutOutEnd   = M_PI *  20/180;

  //This should be calculated or at least cross-checked
  DD4hep::Geometry::Position incomingBeamPipeAtEndOfBeamCalPosition(-incomingBeamPipeRadius, incomingBeamPipeRadius, bcalInnerZ+bcalThickness);
  incomingBeamPipeAtEndOfBeamCalPosition = DD4hep::Geometry::RotationY(-fullCrossingAngle) * incomingBeamPipeAtEndOfBeamCalPosition;
  const double cutOutRadius = (incomingBeamPipeAtEndOfBeamCalPosition.Rho())+0.1;
  std::cout << "cutOutRadius: " << cutOutRadius << " mm " << std::endl;

  //const double cutOutRadius    = 80;

  DD4hep::Geometry::Tube cutOutTube (0.0, cutOutRadius, bcalThickness, bcalCutOutStart, bcalCutOutEnd);

  ////////////////////////////////////////////////////////////////////////////////
  // Create all the layers
  ////////////////////////////////////////////////////////////////////////////////

  //Loop over all the layer (repeat=NN) sections
  //_U is define to create unicodestring, right?
  for(DD4hep::XML::Collection_t coll(xmlBeamCal,_U(layer)); coll; ++coll)  {
    DD4hep::XML::Component xmlLayer(coll); //we know this thing is a layer


    //This just calculates the total size of a single layer
    //Why no convenience function for this?
    double layerThickness = 0;
    for(DD4hep::XML::Collection_t l(xmlLayer,_U(slice)); l; ++l)
      layerThickness += xml_comp_t(l).thickness();

    std::cout << "Total Length " << bcalThickness  << std::endl;
    std::cout << "Layer Thickness " << layerThickness  << std::endl;

    //Loop for repeat=NN
    for(int i=0, repeat=xmlLayer.repeat(); i<repeat; ++i)  {

      std::string layer_name = detName + DD4hep::XML::_toString(thisLayerId,"_layer%d");
      DD4hep::Geometry::Volume layer_vol(layer_name,DD4hep::Geometry::Tube(bcalInnerR,bcalOuterR,layerThickness/2),air);

      int sliceID=1;
      double inThisLayerPosition = -layerThickness*0.5;
      for(DD4hep::XML::Collection_t collSlice(xmlLayer,_U(slice)); collSlice; ++collSlice)  {
	DD4hep::XML::Component compSlice = collSlice;
	const double      sliceThickness = compSlice.thickness();
	const std::string sliceName = layer_name + DD4hep::XML::_toString(sliceID,"slice%d");
	DD4hep::Geometry::Material   slice_mat  = lcdd.material(compSlice.materialStr());

	bool isAbsorberStructure(false);
	try {
	  std::string sliceType = compSlice.attr< std::string >(_U(type));
	  if ( sliceType.compare("absorber") == 0 ){
	    isAbsorberStructure=true;
	  } // else {
	  //   throw std::runtime_error("Unknown type of slice in BeamCal, use \"absorber\" or nothing");
	  // }//Do we want this to fail or not?
	} catch (std::runtime_error &e) {
	  //std::cout << "Catching " << e.what()  << std::endl;
	  //std::cout << e.what()  << std::endl;
	}

	DD4hep::Geometry::Tube sliceBase(bcalInnerR,bcalOuterR,sliceThickness/2);
	DD4hep::Geometry::SubtractionSolid slice_subtracted;


	if(isAbsorberStructure) {
	  //If we have the absorber structure then we create the slice with a hole at the position of the outgoing beam pipe
	  ///In This case we have to know the global position of the slice, because the cutout depends on the outgoing beam pipe position
	  const double thisPositionZ = bcalCentreZ - bcalThickness*0.5 + (thisLayerId-1)*layerThickness+ sliceThickness*0.5;
	  const DD4hep::Geometry::Position thisBPPosition( std::tan(fullCrossingAngle) * thisPositionZ, 0.0, 0.0);
	  //The extra parenthesis are paramount!
	  const DD4hep::Geometry::Transform3D thisBPTransform( incomingBeamPipeRotation, thisBPPosition );
	  slice_subtracted = DD4hep::Geometry::SubtractionSolid(sliceBase, incomingBeamPipe, thisBPTransform);
	} else {
	  //If we do not have the absorber structure then we create the slice with a wedge cutout, i.e, keyhole shape
	  /// Is it better to join two pieces or subtract two pieces?
	  slice_subtracted = DD4hep::Geometry::SubtractionSolid(sliceBase, cutOutTube, DD4hep::Geometry::Transform3D() );
	}

	DD4hep::Geometry::Volume slice_vol (sliceName,slice_subtracted,slice_mat);

	if ( compSlice.isSensitive() )  {
	  sens.setType("calorimeter");
	  slice_vol.setSensitiveDetector(sens);
	}

	slice_vol.setAttributes(lcdd,compSlice.regionStr(),compSlice.limitsStr(),compSlice.visStr());
	DD4hep::Geometry::PlacedVolume pv = layer_vol.placeVolume(slice_vol,
								  DD4hep::Geometry::Position(0,0,inThisLayerPosition+sliceThickness*0.5));
	pv.addPhysVolID("slice",sliceID);
	inThisLayerPosition += sliceThickness;
	++sliceID;
      }//For all slices in this layer

      //Why are we doing this for each layer, this just needs to be done once and then placed multiple times
      //Do we need unique IDs for each piece?
      layer_vol.setVisAttributes(lcdd,xmlLayer.visStr());

      DD4hep::Geometry::Position layer_pos(0,0,-bcalThickness*0.5+(thisLayerId-0.5)*layerThickness);

      DD4hep::Geometry::PlacedVolume pv = envelopeVol.placeVolume(layer_vol,layer_pos);
      pv.addPhysVolID("layer",thisLayerId);

      ++thisLayerId;

    }//for all layers

  }// for all layer collections

  const DD4hep::Geometry::Position bcForwardPos (std::tan(-0.5*fullCrossingAngle)*bcalCentreZ,0.0, bcalCentreZ);
  const DD4hep::Geometry::Position bcBackwardPos(std::tan(-0.5*fullCrossingAngle)*bcalCentreZ,0.0,-bcalCentreZ);
  const DD4hep::Geometry::Rotation3D bcForwardRot ( DD4hep::Geometry::RotationY(-fullCrossingAngle*0.5 ) );//Minus??
  const DD4hep::Geometry::Rotation3D bcBackwardRot( DD4hep::Geometry::RotationZYX ( (M_PI), (M_PI-fullCrossingAngle*0.5), (0.0)));

  // const DD4hep::Geometry::Position bcForwardPos (std::tan(-fullCrossingAngle)*bcalCentreZ,0.0, bcalCentreZ);
  // const DD4hep::Geometry::Position bcBackwardPos(std::tan(-fullCrossingAngle)*bcalCentreZ,0.0,-bcalCentreZ);
  // //const DD4hep::Geometry::Rotation3D bcForwardRot ( (DD4hep::Geometry::RotationY(fullCrossingAngle)) );
  // const DD4hep::Geometry::Rotation3D   bcBackwardRot( DD4hep::Geometry::RotationZYX ( (M_PI), (M_PI-fullCrossingAngle), (0.0)));

  DD4hep::Geometry::DetElement sdet ( detName, xmlBeamCal.id() );
  DD4hep::Geometry::DetElement rdet ( detName, xmlBeamCal.id() );
  DD4hep::Geometry::Volume motherVol = lcdd.pickMotherVolume(sdet);

  DD4hep::Geometry::PlacedVolume pv =
    motherVol.placeVolume(envelopeVol, DD4hep::Geometry::Transform3D( bcForwardRot, bcForwardPos ) );
  pv.addPhysVolID("system",xmlBeamCal.id());
  pv.addPhysVolID("barrel", 1);
  sdet.setPlacement(pv);

  DD4hep::Geometry::PlacedVolume pv2 =
    motherVol.placeVolume(envelopeVol, DD4hep::Geometry::Transform3D( bcBackwardRot, bcBackwardPos ) );
  pv2.addPhysVolID("system",xmlBeamCal.id());
  pv2.addPhysVolID("barrel", 2);
  rdet.setPlacement(pv2);

  return sdet;
}

DECLARE_DETELEMENT(BeamCal,create_detector)
