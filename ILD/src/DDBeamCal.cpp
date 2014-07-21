#include <DD4hep/DetFactoryHelper.h>
#include <XML/Layering.h>

#include <string>


static DD4hep::Geometry::Ref_t create_detector(DD4hep::Geometry::LCDD& lcdd,
					       DD4hep::XML::Handle_t xmlHandle,
					       DD4hep::Geometry::SensitiveDetector sens) {

  std::cout << __PRETTY_FUNCTION__  << std::endl;
  std::cout << "Here is my BeamCal"  << std::endl;
  std::cout << " and this is the sensitive detector: " << &sens  << std::endl;
  sens.setType("calorimeter");
  //Materials
  DD4hep::Geometry::Material air = lcdd.air();

  //Access to the XML File
  DD4hep::XML::DetElement xmlBeamCal = xmlHandle;
  const std::string detName = xmlBeamCal.nameStr();

 //--------------------------------
  DD4hep::Geometry::Assembly assembly( detName + "_assembly"  ) ;
  //--------------------------------

  DD4hep::XML::Dimension dimensions =  xmlBeamCal.dimensions();

  //BeamCal Dimensions
  const double bcalInnerR = dimensions.inner_r();
  const double bcalOuterR = dimensions.outer_r();
  const double bcalInnerZ = dimensions.inner_z();
  const double bcalThickness = DD4hep::Layering(xmlBeamCal).totalThickness();
  const double bcalCentreZ = bcalInnerZ+bcalThickness*0.5;

  // counter for the current layer to be placed
  int thisLayerId = 1;

  //Parameters we have to know about
  DD4hep::XML::Component xmlParameter = xmlBeamCal.child(_Unicode(parameter));
  const double mradFullCrossingAngle  = xmlParameter.attr< double >(_Unicode(crossingangle));
  std::cout << " The crossing angle is: " << mradFullCrossingAngle << " radian"  << std::endl;

  //Create the section cutout for the sensor and readout volumes
  const double bcalCutOutSpan  = xmlParameter.attr< double >(_Unicode(cutoutspanningangle));
  const double bcalCutOutStart = -bcalCutOutSpan*0.5;
  const double bcalCutOutEnd   = bcalCutOutStart + bcalCutOutSpan;
  const double incomingBeamPipeRadius = xmlParameter.attr< double >( _Unicode(incomingbeampiperadius) );

  std::cout << "bcalCutOutSpan  "<< bcalCutOutSpan  << " Radian"<< std::endl;
  std::cout << "bcalCutOutSpan  "<< bcalCutOutSpan*180/M_PI  << " DEGREE"<< std::endl;
  std::cout << "bcalCutOutStart "<< bcalCutOutStart << " Radian"<< std::endl;
  std::cout << "bcalCutOutEnd   "<< bcalCutOutEnd   << " Radian"<< std::endl;
  std::cout << "incommingBeamPipeRadius: "<< incomingBeamPipeRadius  << std::endl;
  ////////////////////////////////////////////////////////////////////////////////
  //Calculations for the position of the incoming beampipe
  ////////////////////////////////////////////////////////////////////////////////

  // this needs the full crossing angle, because the beamCal is centred on the
  // outgoing beampipe, and the incoming beampipe is a full crossing agle away
  const DD4hep::Geometry::Position    incomingBeamPipePosition( std::tan(mradFullCrossingAngle) * bcalCentreZ, 0.0, 0.0);
  //The extra parenthesis are paramount!
  const DD4hep::Geometry::Rotation3D  incomingBeamPipeRotation( ( DD4hep::Geometry::RotationY( mradFullCrossingAngle ) ) );
  const DD4hep::Geometry::Transform3D incomingBPTransform( incomingBeamPipeRotation, incomingBeamPipePosition );

  //Envelope to place the layers in
  DD4hep::Geometry::Tube envelopeTube (bcalInnerR, bcalOuterR, bcalThickness*0.5 );
  DD4hep::Geometry::Tube incomingBeamPipe (0.0, incomingBeamPipeRadius, bcalThickness);//we want this to be longer than the BeamCal
  DD4hep::Geometry::SubtractionSolid envelope (envelopeTube, incomingBeamPipe, incomingBPTransform);
  DD4hep::Geometry::Volume     envelopeVol(detName+"_envelope",envelope,air);

  //This should be calculated or at least cross-checked
#pragma message("What about the special case when the incoming beampipe would be completely inside the outgoing beam pipe?")
  DD4hep::Geometry::Position incomingBeamPipeAtEndOfBeamCalPosition(-incomingBeamPipeRadius, incomingBeamPipeRadius, bcalInnerZ+bcalThickness);
  //This Rotation needs to be the fullCrossing angle, because the incoming beampipe has that much to the outgoing beam pipe
  //And the BeamCal is centred on the outgoing beampipe
  incomingBeamPipeAtEndOfBeamCalPosition = DD4hep::Geometry::RotationY(-mradFullCrossingAngle) * incomingBeamPipeAtEndOfBeamCalPosition;
  const double cutOutRadius = ( incomingBeamPipeAtEndOfBeamCalPosition.Rho() ) + 0.01/*cm*/;
  std::cout << "cutOutRadius: " << cutOutRadius << " cm " << std::endl;

  //we use bcalThickness on purpose to make the subtraction work properly
  DD4hep::Geometry::Tube cutOutTube (0.0, cutOutRadius, bcalThickness, bcalCutOutStart, bcalCutOutEnd);

  ////////////////////////////////////////////////////////////////////////////////
  // Create all the layers
  ////////////////////////////////////////////////////////////////////////////////

  //Loop over all the layer (repeat=NN) sections
  //This is the starting point to place all layers, we need this when we have more than one layer block
  double referencePosition = -bcalThickness*0.5;
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
	  const std::string& sliceType = compSlice.attr< std::string >(_Unicode(layerType));
	  if ( sliceType.compare("holeForIncomingBeampipe") == 0 ){
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
	  //If we have the absorber structure then we create the slice with a
	  //hole at the position of the outgoing beam pipe. In This case we have
	  //to know the global position of the slice, because the cutout depends
	  //on the outgoing beam pipe position
	  const double thisPositionZ = bcalCentreZ - bcalThickness*0.5 + (thisLayerId-1)*layerThickness+ sliceThickness*0.5;
	  const DD4hep::Geometry::Position thisBPPosition( std::tan(mradFullCrossingAngle) * thisPositionZ, 0.0, 0.0);
	  //The extra parenthesis are paramount! But.. there are none
	  const DD4hep::Geometry::Transform3D thisBPTransform( incomingBeamPipeRotation, thisBPPosition );
	  slice_subtracted = DD4hep::Geometry::SubtractionSolid(sliceBase, incomingBeamPipe, thisBPTransform);
	} else {
	  //If we do not have the absorber structure then we create the slice with a wedge cutout, i.e, keyhole shape
	  /// Is it better to join two pieces or subtract two pieces?
	  slice_subtracted = DD4hep::Geometry::SubtractionSolid(sliceBase, cutOutTube, DD4hep::Geometry::Transform3D() );
	}

	DD4hep::Geometry::Volume slice_vol (sliceName,slice_subtracted,slice_mat);

	if ( compSlice.isSensitive() )  {
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

      DD4hep::Geometry::Position layer_pos(0,0,referencePosition+0.5*layerThickness);
      referencePosition += layerThickness;
      DD4hep::Geometry::PlacedVolume pv = envelopeVol.placeVolume(layer_vol,layer_pos);
      pv.addPhysVolID("layer",thisLayerId);

      ++thisLayerId;

    }//for all layers

  }// for all layer collections

  const DD4hep::Geometry::Position bcForwardPos (std::tan(-0.5*mradFullCrossingAngle)*bcalCentreZ,0.0, bcalCentreZ);
  const DD4hep::Geometry::Position bcBackwardPos(std::tan(-0.5*mradFullCrossingAngle)*bcalCentreZ,0.0,-bcalCentreZ);
  const DD4hep::Geometry::Rotation3D bcForwardRot ( DD4hep::Geometry::RotationY(-mradFullCrossingAngle*0.5 ) );
  const DD4hep::Geometry::Rotation3D bcBackwardRot( DD4hep::Geometry::RotationZYX ( (M_PI), (M_PI-mradFullCrossingAngle*0.5), (0.0)));

  DD4hep::Geometry::DetElement beamcals ( detName, xmlBeamCal.id() );
  DD4hep::Geometry::DetElement sdet ( "BeamCal01", xmlBeamCal.id() );
  DD4hep::Geometry::DetElement rdet ( "BeamCal02", xmlBeamCal.id() );

  DD4hep::Geometry::Volume motherVol = lcdd.pickMotherVolume(sdet);

  DD4hep::Geometry::PlacedVolume pv =
    assembly.placeVolume(envelopeVol, DD4hep::Geometry::Transform3D( bcForwardRot, bcForwardPos ) );
  pv.addPhysVolID("system",xmlBeamCal.id());
  pv.addPhysVolID("barrel", 1);
  sdet.setPlacement(pv);

  DD4hep::Geometry::PlacedVolume pv2 =
    assembly.placeVolume(envelopeVol, DD4hep::Geometry::Transform3D( bcBackwardRot, bcBackwardPos ) );
  pv2.addPhysVolID("system",xmlBeamCal.id());
  pv2.addPhysVolID("barrel", 2);
  rdet.setPlacement(pv2);
 
  beamcals.add(sdet);
  beamcals.add(rdet);

  pv = motherVol.placeVolume( assembly ) ;
  pv.addPhysVolID("system",xmlBeamCal.id());
  beamcals.setPlacement( pv ) ;

  return beamcals;
}

DECLARE_DETELEMENT(BeamCal,create_detector)
