#include <DD4hep/DetFactoryHelper.h>
#include <XML/Layering.h>

#include <string>

using dd4hep::Assembly;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::SensitiveDetector;
using dd4hep::Tube;
using dd4hep::Volume;

static Ref_t create_detector(Detector& theDetector,
                             xml_h xmlHandle,
                             SensitiveDetector sens) {

  std::cout << __PRETTY_FUNCTION__  << std::endl;
  std::cout << "Here is my TestBeamSetup"  << std::endl;
  std::cout << " and this is the sensitive detector: " << &sens  << std::endl;
  sens.setType("calorimeter");
  //Materials
  Material air = theDetector.air();

  //Access to the XML File
  xml_det_t xmlTestBeamSetup = xmlHandle;
  const std::string detName = xmlTestBeamSetup.nameStr();

 //--------------------------------
  Assembly assembly( detName + "_assembly"  ) ;
  //--------------------------------

  //Parameters we have to know about
  dd4hep::xml::Component xmlParameter = xmlTestBeamSetup.child(_Unicode(parameter));
  const double sensorCentreOffset = xmlParameter.attr< double >(_Unicode(sensorCentreOffset));
  const double squareSideLength   = xmlParameter.attr< double >(_Unicode(squareSideLength));

  dd4hep::xml::Dimension dimensions =  xmlTestBeamSetup.dimensions();

  //TestBeamSetup Dimensions
  const double lcalInnerR = dimensions.inner_r();
  const double lcalOuterR = dimensions.outer_r();
  const double lcalInnerZ = dimensions.inner_z();
  const double lcalThickness = Layering( xmlTestBeamSetup ).totalThickness();
  const double lcalCentreZ = lcalInnerZ+lcalThickness*0.5;
  const double startAngle = -15*dd4hep::deg+90*dd4hep::degree;
  const double endAngle   =  15*dd4hep::deg+90*dd4hep::degree;

  const double offsetArc    = (lcalOuterR+lcalInnerR)*0.5;
  const double centreSquare = (lcalOuterR-sensorCentreOffset)-(lcalOuterR+lcalInnerR)*0.5;
  // counter for the current layer to be placed
  int thisLayerId = 1;



  //Envelope to place the layers in
  Box envelopeBox (squareSideLength*0.5, squareSideLength*0.5, lcalThickness*0.5);
  Volume envelopeVol(detName+"_envelope",envelopeBox,air);
  envelopeVol.setVisAttributes(theDetector,xmlTestBeamSetup.visStr());

  ////////////////////////////////////////////////////////////////////////////////
  // Create all the layers
  ////////////////////////////////////////////////////////////////////////////////

  //Loop over all the layer (repeat=NN) sections
  //This is the starting point to place all layers, we need this when we have more than one layer block
  double referencePosition = -lcalThickness*0.5;
  for(dd4hep::xml::Collection_t coll(xmlTestBeamSetup,_U(layer)); coll; ++coll)  {
    dd4hep::xml::Component xmlLayer(coll); //we know this thing is a layer


    //This just calculates the total size of a single layer
    //Why no convenience function for this?
    double layerThickness = 0;
    for(dd4hep::xml::Collection_t l(xmlLayer,_U(slice)); l; ++l)
      layerThickness += xml_comp_t(l).thickness();

    std::cout << "Total Length "    << lcalThickness/dd4hep::cm  << " cm" << std::endl;
    std::cout << "Layer Thickness " << layerThickness/dd4hep::cm << " cm" << std::endl;

    //Loop for repeat=NN
    for(int i=0, repeat=xmlLayer.repeat(); i<repeat; ++i)  {

      std::string layer_name = detName + dd4hep::xml::_toString(thisLayerId,"_layer%d");
      Box layer_base(squareSideLength*0.5, squareSideLength*0.5, layerThickness*0.5);

      Volume layer_vol(layer_name,layer_base,air);

      int sliceID=1;
      double inThisLayerPosition = -layerThickness*0.5;
      for(dd4hep::xml::Collection_t collSlice(xmlLayer,_U(slice)); collSlice; ++collSlice)  {
	dd4hep::xml::Component compSlice = collSlice;
	const double      sliceThickness = compSlice.thickness();
	const std::string sliceName = layer_name + dd4hep::xml::_toString(sliceID,"slice%d");
	Material   slice_mat  = theDetector.material(compSlice.materialStr());


	const std::string& sliceType = compSlice.attr< std::string >(_Unicode(type));

	if ( sliceType.compare("Square") == 0 ){

	  Box sliceBase(squareSideLength*0.5,squareSideLength*0.5,sliceThickness/2);
	  Volume slice_vol (sliceName,sliceBase,slice_mat);

	  if ( compSlice.isSensitive() )  {
	    slice_vol.setSensitiveDetector(sens);
	  }

	  slice_vol.setAttributes(theDetector,compSlice.regionStr(),compSlice.limitsStr(),compSlice.visStr());
	  PlacedVolume pv =
	    layer_vol.placeVolume(slice_vol,
				  Position(0,0.0,inThisLayerPosition+sliceThickness*0.5));
	  pv.addPhysVolID("slice",sliceID);

	} else {


	  Tube sliceBase(lcalInnerR,lcalOuterR,sliceThickness/2, startAngle, endAngle);
	  Volume slice_vol (sliceName,sliceBase,slice_mat);

	  if ( compSlice.isSensitive() )  {
	    slice_vol.setSensitiveDetector(sens);
	  }

	  slice_vol.setAttributes(theDetector,compSlice.regionStr(),compSlice.limitsStr(),compSlice.visStr());
	  PlacedVolume pv = layer_vol.placeVolume(slice_vol,
								    Position(0,-offsetArc,inThisLayerPosition+sliceThickness*0.5));
	  pv.addPhysVolID("slice",sliceID);

      }//not a square





	inThisLayerPosition += sliceThickness;
	++sliceID;
      }//For all slices in this layer

      //Why are we doing this for each layer, this just needs to be done once and then placed multiple times
      //Do we need unique IDs for each piece?
      layer_vol.setVisAttributes(theDetector,xmlLayer.visStr());

      Position layer_pos(0,0,referencePosition+0.5*layerThickness);
      referencePosition += layerThickness;
      PlacedVolume pv = envelopeVol.placeVolume(layer_vol,layer_pos);
      pv.addPhysVolID("layer",thisLayerId);

      ++thisLayerId;

    }//for all layers

  }// for all layer collections

  const Position bcForwardPos ( 0.0 , -centreSquare, lcalCentreZ);

  DetElement TestBeamSetups ( detName, xmlTestBeamSetup.id() );
  DetElement sdet ( "TestBeamSetup01", xmlTestBeamSetup.id() );

  Volume motherVol = theDetector.pickMotherVolume(sdet);

  PlacedVolume pv =
    assembly.placeVolume(envelopeVol, bcForwardPos );
  pv.addPhysVolID("system",xmlTestBeamSetup.id());
  pv.addPhysVolID("barrel", 1);
  sdet.setPlacement(pv);

  TestBeamSetups.add(sdet);

  pv = motherVol.placeVolume( assembly ) ;
  pv.addPhysVolID("system",xmlTestBeamSetup.id());
  TestBeamSetups.setPlacement( pv ) ;

  return TestBeamSetups;
}

DECLARE_DETELEMENT(FCalTB_2014,create_detector)
