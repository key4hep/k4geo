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

static dd4hep::Ref_t create_detector(dd4hep::Detector& theDetector,
                                     dd4hep::xml::Handle_t element,
                                     dd4hep::SensitiveDetector sens) {

  std::cout << "This is the BeamCal"  << std::endl;
  sens.setType("calorimeter");

  //Access to the XML File
  dd4hep::xml::DetElement xmlBeamCal  = element;
  std::string   detName     = xmlBeamCal.nameStr();

  dd4hep::DetElement sdet ( detName, xmlBeamCal.id() );

  // --- create an envelope volume and position it into the world ---------------------
  
  dd4hep::Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;
  
  dd4hep::xml::setDetectorTypeFlag( element, sdet ) ;
  
  if( theDetector.buildType() == dd4hep::BUILD_ENVELOPE ) { std::cout << "Building envelope.\n"; return sdet ; }
  else { std::cout << "Building all.\n"; }
  
  //-----------------------------------------------------------------------------------




  dd4hep::xml::Dimension dimensions =  xmlBeamCal.dimensions();

  //BeamCal Dimensions
  const double bcalInnerR = dimensions.inner_r();
  const double bcalOuterR = dimensions.outer_r();
  const double bcalInnerZ = dimensions.inner_z();
  const double bcalThickness = dd4hep::Layering(xmlBeamCal).totalThickness();
  const double bcalCentreZ = bcalInnerZ+bcalThickness*0.5;

  double BeamCal_cell_size      = theDetector.constant<double>("BeamCal_cell_size");
  //========== fill data for reconstruction ============================
  dd4hep::rec::LayeredCalorimeterData* caloData = new dd4hep::rec::LayeredCalorimeterData ;
  caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::EndcapLayout ;
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
  dd4hep::xml::Component xmlParameter = xmlBeamCal.child(_Unicode(parameter));
  const double fullCrossingAngle  = xmlParameter.attr< double >(_Unicode(crossingangle));
  std::cout << " The crossing angle is: " << fullCrossingAngle << " radian"  << std::endl;


  //Create the section cutout for the sensor and readout volumes
  // The cutout sits on the negative x-axis, i.e. phi~=180degrees
  const double bcalCutOutSpan  = xmlParameter.attr< double >(_Unicode(cutoutspanningangle));
  const double bcalCutOutStart = 180.0*dd4hep::degree-bcalCutOutSpan*0.5;
  const double bcalCutOutEnd   = bcalCutOutStart + bcalCutOutSpan;
  const double incomingBeamPipeRadius = xmlParameter.attr< double >( _Unicode(incomingbeampiperadius) );

  std::cout << "bcalCutOutSpan  "<< bcalCutOutSpan/dd4hep::mrad  << " mrad"<< std::endl;
  std::cout << "bcalCutOutSpan  "<< bcalCutOutSpan/dd4hep::degree  << " degree"<< std::endl;
  std::cout << "bcalCutOutStart "<< bcalCutOutStart << " Radian"<< std::endl;
  std::cout << "bcalCutOutEnd   "<< bcalCutOutEnd   << " Radian"<< std::endl;
  std::cout << "incommingBeamPipeRadius: "<< incomingBeamPipeRadius/dd4hep::cm << " cm"  << std::endl;
  ////////////////////////////////////////////////////////////////////////////////
  //Calculations for the position of the incoming beampipe
  ////////////////////////////////////////////////////////////////////////////////

  // this needs the full crossing angle, because the beamCal is centred on the
  // outgoing beampipe, and the incoming beampipe is a full crossing angle away
  const dd4hep::Position    incomingBeamPipePosition( std::tan(-fullCrossingAngle) * bcalCentreZ, 0.0, 0.0);
  //The extra parenthesis are paramount!
  const dd4hep::Rotation3D  incomingBeamPipeRotation( ( dd4hep::RotationY( -fullCrossingAngle ) ) );
  const dd4hep::Transform3D incomingBPTransform( incomingBeamPipeRotation, incomingBeamPipePosition );

  //Envelope to place the layers in
  
  // dd4hep::Tube envelopeTube (bcalInnerR, bcalOuterR, bcalThickness*0.5 );
  dd4hep::Tube incomingBeamPipe (0.0, incomingBeamPipeRadius, bcalThickness);//we want this to be longer than the BeamCal

  dd4hep::Assembly     envelopeVol(detName+"_module");
  dd4hep::DetElement beamCalDE_1(sdet,"Calorimeter1",1);
  dd4hep::DetElement beamCalDE_2(sdet,"Calorimeter2",2);

  envelopeVol.setVisAttributes(theDetector,xmlBeamCal.visStr());
  
  dd4hep::Position incomingBeamPipeAtEndOfBeamCalPosition(-incomingBeamPipeRadius, incomingBeamPipeRadius, bcalInnerZ+bcalThickness);
  //This Rotation needs to be the fullCrossing angle, because the incoming beampipe has that much to the outgoing beam pipe
  //And the BeamCal is centred on the outgoing beampipe
  incomingBeamPipeAtEndOfBeamCalPosition = dd4hep::RotationY(-fullCrossingAngle) * incomingBeamPipeAtEndOfBeamCalPosition;
  const double cutOutRadius = ( incomingBeamPipeAtEndOfBeamCalPosition.Rho() ) + 0.01*dd4hep::cm;
  std::cout << "cutOutRadius: " << cutOutRadius/dd4hep::cm << " cm " << std::endl;

  //we use bcalThickness on purpose to make the subtraction work properly
  dd4hep::Tube cutOutTube (0.0, cutOutRadius, bcalThickness, bcalCutOutStart, bcalCutOutEnd);

  ////////////////////////////////////////////////////////////////////////////////
  // Create all the layers
  ////////////////////////////////////////////////////////////////////////////////

  //Loop over all the layer (repeat=NN) sections
  //This is the starting point to place all layers, we need this when we have more than one layer block
  double referencePosition = -bcalThickness*0.5;
  for(dd4hep::xml::Collection_t coll(xmlBeamCal,_U(layer)); coll; ++coll)  {
    dd4hep::xml::Component xmlLayer(coll); //we know this thing is a layer


    //This just calculates the total size of a single layer
    //Why no convenience function for this?
    double layerThickness = 0;
    for(dd4hep::xml::Collection_t l(xmlLayer,_U(slice)); l; ++l)
      layerThickness += xml_comp_t(l).thickness();

    std::cout << "Total Length "    << bcalThickness/dd4hep::cm  << " cm" << std::endl;
    std::cout << "Layer Thickness " << layerThickness/dd4hep::cm << " cm" << std::endl;
    dd4hep::rec::LayeredCalorimeterData::Layer caloLayer ;
            
    //Loop for repeat=NN
    for(int i=0, repeat=xmlLayer.repeat(); i<repeat; ++i)  {
      double nRadiationLengths=0.;
      double nInteractionLengths=0.;
      double thickness_sum=0;

      std::string layer_name = detName + dd4hep::xml::_toString(thisLayerId,"_layer%d");

      dd4hep::Assembly layer_vol(layer_name) ;

      int sliceID=1;
      double inThisLayerPosition = -layerThickness*0.5;


      for(dd4hep::xml::Collection_t collSlice(xmlLayer,_U(slice)); collSlice; ++collSlice)  {
	dd4hep::xml::Component compSlice = collSlice;
	const double      slice_thickness = compSlice.thickness();
	const std::string sliceName = layer_name + dd4hep::xml::_toString(sliceID,"slice%d");
	dd4hep::Material   slice_material  = theDetector.material(compSlice.materialStr());

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

	// Check if a separate outer_radius is declared.
	double outerR = bcalOuterR;
	try {
	  outerR = compSlice.outer_radius();
	} catch (std::runtime_error &e) {
	  // Nothing to catch. Everything is fine.
	}

	dd4hep::Tube sliceBase(bcalInnerR, outerR, slice_thickness/2);
	dd4hep::SubtractionSolid slice_subtracted;

	if(isAbsorberStructure) {
	  //If we have the absorber structure then we create the slice with a
	  //hole at the position of the outgoing beam pipe. In This case we have
	  //to know the global position of the slice, because the cutout depends
	  //on the outgoing beam pipe position
	  const double thisPositionZ = bcalCentreZ + referencePosition + 0.5*layerThickness + inThisLayerPosition + slice_thickness*0.5;
	  const dd4hep::Position thisBPPosition( std::tan(-fullCrossingAngle) * thisPositionZ, 0.0, 0.0);
	  //The extra parenthesis are paramount! But.. there are none
	  const dd4hep::Transform3D thisBPTransform( incomingBeamPipeRotation, thisBPPosition );
	  slice_subtracted = dd4hep::SubtractionSolid(sliceBase, incomingBeamPipe, thisBPTransform);
	} else {
	  //If we do not have the absorber structure then we create the slice with a wedge cutout, i.e, keyhole shape
	  /// Is it better to join two pieces or subtract two pieces?
	  slice_subtracted = dd4hep::SubtractionSolid(sliceBase, cutOutTube, dd4hep::Transform3D() );
	}

	dd4hep::Volume slice_vol (sliceName,slice_subtracted,slice_material);

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
                	

	slice_vol.setAttributes(theDetector,compSlice.regionStr(),compSlice.limitsStr(),compSlice.visStr());
	dd4hep::PlacedVolume pv = layer_vol.placeVolume(slice_vol,
								  dd4hep::Position(0,0,inThisLayerPosition+slice_thickness*0.5));
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
      layer_vol.setVisAttributes(theDetector,xmlLayer.visStr());

      dd4hep::Position layer_pos(0,0,referencePosition+0.5*layerThickness);
      referencePosition += layerThickness;
      dd4hep::PlacedVolume pv = envelopeVol.placeVolume(layer_vol,layer_pos);
      pv.addPhysVolID("layer",thisLayerId);

      ++thisLayerId;

    }//for all layers

  }// for all layer collections

  const dd4hep::Position bcForwardPos (std::tan(0.5*fullCrossingAngle)*bcalCentreZ,0.0, bcalCentreZ);
  const dd4hep::Position bcBackwardPos(std::tan(0.5*fullCrossingAngle)*bcalCentreZ,0.0,-bcalCentreZ);
  const dd4hep::Rotation3D bcForwardRot ( dd4hep::RotationY(+fullCrossingAngle*0.5 ) );
  const dd4hep::Rotation3D bcBackwardRot( dd4hep::RotationZYX ( (180*dd4hep::degree), (180*dd4hep::degree-fullCrossingAngle*0.5), (0.0)));
  
  dd4hep::PlacedVolume pv =
    envelope.placeVolume(envelopeVol, dd4hep::Transform3D( bcForwardRot, bcForwardPos ) );
  pv.addPhysVolID("barrel", 1);
  beamCalDE_1.setPlacement(pv);
  dd4hep::PlacedVolume pv2 =
    envelope.placeVolume(envelopeVol, dd4hep::Transform3D( bcBackwardRot, bcBackwardPos ) );
  pv2.addPhysVolID("barrel", 2);
  beamCalDE_2.setPlacement(pv2);

  sdet.addExtension< dd4hep::rec::LayeredCalorimeterData >( caloData ) ;

  return sdet;
}

DECLARE_DETELEMENT(BeamCal_o1_v02,create_detector)
