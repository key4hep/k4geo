#include <DD4hep/DetFactoryHelper.h>
#include "DD4hep/DetType.h"
#include <XML/Layering.h>
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

#include <string>

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

static DD4hep::Geometry::Ref_t create_detector(DD4hep::Geometry::LCDD& lcdd,
					       DD4hep::XML::Handle_t element,
					       DD4hep::Geometry::SensitiveDetector sens) {

  std::cout << "This is the BeamCal"  << std::endl;
  sens.setType("calorimeter");
  //Materials
  DD4hep::Geometry::Material air = lcdd.air();

  //Access to the XML File
  DD4hep::XML::DetElement xmlBeamCal  = element;
  std::string   detName     = xmlBeamCal.nameStr();

  DD4hep::Geometry::DetElement sdet ( detName, xmlBeamCal.id() );

  // --- create an envelope volume and position it into the world ---------------------
  
  DD4hep::Geometry::Volume envelope = DD4hep::XML::createPlacedEnvelope( lcdd,  element , sdet ) ;
  
  DD4hep::XML::setDetectorTypeFlag( element, sdet ) ;
  
  if( lcdd.buildType() == DD4hep::BUILD_ENVELOPE ) return sdet ;
  
  //-----------------------------------------------------------------------------------




  DD4hep::XML::Dimension dimensions =  xmlBeamCal.dimensions();

  //BeamCal Dimensions
  const double bcalInnerR = dimensions.inner_r();
  const double bcalOuterR = dimensions.outer_r();
  const double bcalInnerZ = dimensions.inner_z();
  const double bcalThickness = DD4hep::Layering(xmlBeamCal).totalThickness();
  const double bcalCentreZ = bcalInnerZ+bcalThickness*0.5;

  double BeamCal_cell_size      = lcdd.constant<double>("BeamCal_cell_size");
  //========== fill data for reconstruction ============================
  DD4hep::DDRec::LayeredCalorimeterData* caloData = new DD4hep::DDRec::LayeredCalorimeterData ;
  caloData->layoutType = DD4hep::DDRec::LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = 0  ; // hardcoded tube
  caloData->outer_symmetry = 0  ;
  caloData->phi0 = 0 ;

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = bcalInnerR ;
  caloData->extent[1] = bcalOuterR ;
  caloData->extent[2] = bcalInnerZ ;
  caloData->extent[3] = bcalInnerZ + bcalThickness ;

  // counter for the current layer to be placed
  int thisLayerId = 1;

  //Parameters we have to know about
  DD4hep::XML::Component xmlParameter = xmlBeamCal.child(_Unicode(parameter));
  const double fullCrossingAngle  = xmlParameter.attr< double >(_Unicode(crossingangle));
  std::cout << " The crossing angle is: " << fullCrossingAngle << " radian"  << std::endl;

  //Create the section cutout for the sensor and readout volumes
  // The cutout sits on the negative x-axis, i.e. phi~=180degrees
  const double bcalCutOutSpan  = xmlParameter.attr< double >(_Unicode(cutoutspanningangle));
  const double bcalCutOutStart = 180.0*dd4hep::degree-bcalCutOutSpan*0.5;
  const double bcalCutOutEnd   = bcalCutOutStart + bcalCutOutSpan;
  const double incomingBeamPipeRadius = xmlParameter.attr< double >( _Unicode(incomingbeampiperadius) );

  std::cout << "bcalCutOutSpan  "<< bcalCutOutSpan/dd4hep::mrad  << " Radian"<< std::endl;
  std::cout << "bcalCutOutSpan  "<< bcalCutOutSpan/dd4hep::degree  << " DEGREE"<< std::endl;
  std::cout << "bcalCutOutStart "<< bcalCutOutStart << " Radian"<< std::endl;
  std::cout << "bcalCutOutEnd   "<< bcalCutOutEnd   << " Radian"<< std::endl;
  std::cout << "incommingBeamPipeRadius: "<< incomingBeamPipeRadius/dd4hep::cm << " cm"  << std::endl;
  ////////////////////////////////////////////////////////////////////////////////
  //Calculations for the position of the incoming beampipe
  ////////////////////////////////////////////////////////////////////////////////

  // this needs the full crossing angle, because the beamCal is centred on the
  // outgoing beampipe, and the incoming beampipe is a full crossing angle away
  const DD4hep::Geometry::Position    incomingBeamPipePosition( std::tan(-fullCrossingAngle) * bcalCentreZ, 0.0, 0.0);
  //The extra parenthesis are paramount!
  const DD4hep::Geometry::Rotation3D  incomingBeamPipeRotation( ( DD4hep::Geometry::RotationY( -fullCrossingAngle ) ) );
  const DD4hep::Geometry::Transform3D incomingBPTransform( incomingBeamPipeRotation, incomingBeamPipePosition );

  //Envelope to place the layers in
  
  DD4hep::Geometry::Tube envelopeTube (bcalInnerR, bcalOuterR, bcalThickness*0.5 );
  DD4hep::Geometry::Tube incomingBeamPipe (0.0, incomingBeamPipeRadius, bcalThickness);//we want this to be longer than the BeamCal
  DD4hep::Geometry::SubtractionSolid BeamCalModule (envelopeTube, incomingBeamPipe, incomingBPTransform);
  DD4hep::Geometry::Volume     envelopeVol(detName+"_module",BeamCalModule,air);
  envelopeVol.setVisAttributes(lcdd,xmlBeamCal.visStr());
  
  DD4hep::Geometry::Position incomingBeamPipeAtEndOfBeamCalPosition(-incomingBeamPipeRadius, incomingBeamPipeRadius, bcalInnerZ+bcalThickness);
  //This Rotation needs to be the fullCrossing angle, because the incoming beampipe has that much to the outgoing beam pipe
  //And the BeamCal is centred on the outgoing beampipe
  incomingBeamPipeAtEndOfBeamCalPosition = DD4hep::Geometry::RotationY(-fullCrossingAngle) * incomingBeamPipeAtEndOfBeamCalPosition;
  const double cutOutRadius = ( incomingBeamPipeAtEndOfBeamCalPosition.Rho() ) + 0.01*dd4hep::cm;
  std::cout << "cutOutRadius: " << cutOutRadius/dd4hep::cm << " cm " << std::endl;

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

    std::cout << "Total Length "    << bcalThickness/dd4hep::cm  << " cm" << std::endl;
    std::cout << "Layer Thickness " << layerThickness/dd4hep::cm << " cm" << std::endl;
    DD4hep::DDRec::LayeredCalorimeterData::Layer caloLayer ;
    double nRadiationLengths=0.;
    double nInteractionLengths=0.;
    double thickness_sum=0;
            
    //Loop for repeat=NN
    for(int i=0, repeat=xmlLayer.repeat(); i<repeat; ++i)  {

      std::string layer_name = detName + DD4hep::XML::_toString(thisLayerId,"_layer%d");
      DD4hep::Geometry::Tube layer_base(bcalInnerR,bcalOuterR,layerThickness*0.5);
      DD4hep::Geometry::SubtractionSolid layer_subtracted;
      { // put this in extra block to limit scope
	const double thisPositionZ = bcalCentreZ + referencePosition + layerThickness*0.5;
	const DD4hep::Geometry::Position thisBPPosition( std::tan(-fullCrossingAngle) * thisPositionZ, 0.0, 0.0);
	const DD4hep::Geometry::Transform3D thisBPTransform( incomingBeamPipeRotation, thisBPPosition );
	layer_subtracted = DD4hep::Geometry::SubtractionSolid(layer_base, incomingBeamPipe, thisBPTransform);
      }

      DD4hep::Geometry::Volume layer_vol(layer_name,layer_subtracted,air);


      int sliceID=1;
      double inThisLayerPosition = -layerThickness*0.5;

      double radiator_thickness = 0.0;

      for(DD4hep::XML::Collection_t collSlice(xmlLayer,_U(slice)); collSlice; ++collSlice)  {
	DD4hep::XML::Component compSlice = collSlice;
	const double      slice_thickness = compSlice.thickness();
	const std::string sliceName = layer_name + DD4hep::XML::_toString(sliceID,"slice%d");
	DD4hep::Geometry::Material   slice_material  = lcdd.material(compSlice.materialStr());

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

	DD4hep::Geometry::Tube sliceBase(bcalInnerR,bcalOuterR,slice_thickness/2);
	DD4hep::Geometry::SubtractionSolid slice_subtracted;


	if(isAbsorberStructure) {
	  //If we have the absorber structure then we create the slice with a
	  //hole at the position of the outgoing beam pipe. In This case we have
	  //to know the global position of the slice, because the cutout depends
	  //on the outgoing beam pipe position
	  const double thisPositionZ = bcalCentreZ + referencePosition + 0.5*layerThickness + inThisLayerPosition + slice_thickness*0.5;
	  const DD4hep::Geometry::Position thisBPPosition( std::tan(-fullCrossingAngle) * thisPositionZ, 0.0, 0.0);
	  //The extra parenthesis are paramount! But.. there are none
	  const DD4hep::Geometry::Transform3D thisBPTransform( incomingBeamPipeRotation, thisBPPosition );
	  slice_subtracted = DD4hep::Geometry::SubtractionSolid(sliceBase, incomingBeamPipe, thisBPTransform);
	  radiator_thickness = slice_thickness;
	} else {
	  //If we do not have the absorber structure then we create the slice with a wedge cutout, i.e, keyhole shape
	  /// Is it better to join two pieces or subtract two pieces?
	  slice_subtracted = DD4hep::Geometry::SubtractionSolid(sliceBase, cutOutTube, DD4hep::Geometry::Transform3D() );
	}

	DD4hep::Geometry::Volume slice_vol (sliceName,slice_subtracted,slice_material);

        nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
        nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
        thickness_sum += slice_thickness/2;
                
	if ( compSlice.isSensitive() )  {
	  slice_vol.setSensitiveDetector(sens);
          
#if DD4HEP_VERSION_GE( 0, 15 )
          //Store "inner" quantities
          caloLayer.inner_nRadiationLengths = nRadiationLengths;
          caloLayer.inner_nInteractionLengths = nInteractionLengths;
          caloLayer.inner_thickness = thickness_sum;
          //Store scintillator thickness
          caloLayer.sensitive_thickness = slice_thickness;
#endif      
          //Reset counters to measure "outside" quantitites
          nRadiationLengths=0.;
          nInteractionLengths=0.;
          thickness_sum = 0.;
	}
	
        nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
        nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
        thickness_sum += slice_thickness/2;
                	

	slice_vol.setAttributes(lcdd,compSlice.regionStr(),compSlice.limitsStr(),compSlice.visStr());
	DD4hep::Geometry::PlacedVolume pv = layer_vol.placeVolume(slice_vol,
								  DD4hep::Geometry::Position(0,0,inThisLayerPosition+slice_thickness*0.5));
	pv.addPhysVolID("slice",sliceID);
	inThisLayerPosition += slice_thickness;
	++sliceID;
      }//For all slices in this layer

      //-----------------------------------------------------------------------------------------
      
      caloLayer.distance = bcalCentreZ + referencePosition;

#if DD4HEP_VERSION_GE( 0, 15 )
      caloLayer.outer_nRadiationLengths = nRadiationLengths;
      caloLayer.outer_nInteractionLengths = nInteractionLengths;
      caloLayer.outer_thickness = thickness_sum;
#endif            
      caloLayer.cellSize0 = BeamCal_cell_size ;
      caloLayer.cellSize1 = BeamCal_cell_size ;
      
      caloData->layers.push_back( caloLayer ) ;
      //-----------------------------------------------------------------------------------------

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

  const DD4hep::Geometry::Position bcForwardPos (std::tan(0.5*fullCrossingAngle)*bcalCentreZ,0.0, bcalCentreZ);
  const DD4hep::Geometry::Position bcBackwardPos(std::tan(0.5*fullCrossingAngle)*bcalCentreZ,0.0,-bcalCentreZ);
  const DD4hep::Geometry::Rotation3D bcForwardRot ( DD4hep::Geometry::RotationY(+fullCrossingAngle*0.5 ) );
  const DD4hep::Geometry::Rotation3D bcBackwardRot( DD4hep::Geometry::RotationZYX ( (180*dd4hep::degree), (180*dd4hep::degree-fullCrossingAngle*0.5), (0.0)));
  
  DD4hep::Geometry::PlacedVolume pv =
    envelope.placeVolume(envelopeVol, DD4hep::Geometry::Transform3D( bcForwardRot, bcForwardPos ) );
  pv.addPhysVolID("barrel", 1);


  DD4hep::Geometry::PlacedVolume pv2 =
    envelope.placeVolume(envelopeVol, DD4hep::Geometry::Transform3D( bcBackwardRot, bcBackwardPos ) );
  pv2.addPhysVolID("barrel", 2);

  sdet.addExtension< DD4hep::DDRec::LayeredCalorimeterData >( caloData ) ;

  return sdet;
}

DECLARE_DETELEMENT(BeamCal_o1_v01,create_detector)
