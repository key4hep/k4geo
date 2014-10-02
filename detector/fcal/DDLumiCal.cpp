#include <DD4hep/DetFactoryHelper.h>
#include <XML/Layering.h>

#include <string>


static DD4hep::Geometry::Ref_t create_detector(DD4hep::Geometry::LCDD& lcdd,
					       DD4hep::XML::Handle_t xmlHandle,
					       DD4hep::Geometry::SensitiveDetector sens) {

  std::cout << __PRETTY_FUNCTION__  << std::endl;
  std::cout << "Here is my LumiCal"  << std::endl;
  std::cout << " and this is the sensitive detector: " << &sens  << std::endl;
  sens.setType("calorimeter");
  //Materials
  DD4hep::Geometry::Material air = lcdd.air();

  //Access to the XML File
  DD4hep::XML::DetElement xmlLumiCal = xmlHandle;
  const std::string detName = xmlLumiCal.nameStr();

 //--------------------------------
  DD4hep::Geometry::Assembly assembly( detName + "_assembly"  ) ;
  //--------------------------------

  DD4hep::XML::Dimension dimensions =  xmlLumiCal.dimensions();

  //LumiCal Dimensions
  const double lcalInnerR = dimensions.inner_r();
  const double lcalOuterR = dimensions.outer_r();
  const double lcalInnerZ = dimensions.inner_z();
  const double lcalThickness = DD4hep::Layering(xmlLumiCal).totalThickness();
  const double lcalCentreZ = lcalInnerZ+lcalThickness*0.5;

  // counter for the current layer to be placed
  int thisLayerId = 1;

  //Parameters we have to know about
  DD4hep::XML::Component xmlParameter = xmlLumiCal.child(_Unicode(parameter));
  const double mradFullCrossingAngle  = xmlParameter.attr< double >(_Unicode(crossingangle));
  std::cout << " The crossing angle is: " << mradFullCrossingAngle << " radian"  << std::endl;


  //Envelope to place the layers in
  DD4hep::Geometry::Tube envelopeTube (lcalInnerR, lcalOuterR, lcalThickness*0.5 );
  DD4hep::Geometry::Volume     envelopeVol(detName+"_envelope",envelopeTube,air);
  envelopeVol.setVisAttributes(lcdd,xmlLumiCal.visStr());

  ////////////////////////////////////////////////////////////////////////////////
  // Create all the layers
  ////////////////////////////////////////////////////////////////////////////////

  //Loop over all the layer (repeat=NN) sections
  //This is the starting point to place all layers, we need this when we have more than one layer block
  double referencePosition = -lcalThickness*0.5;
  for(DD4hep::XML::Collection_t coll(xmlLumiCal,_U(layer)); coll; ++coll)  {
    DD4hep::XML::Component xmlLayer(coll); //we know this thing is a layer


    //This just calculates the total size of a single layer
    //Why no convenience function for this?
    double layerThickness = 0;
    for(DD4hep::XML::Collection_t l(xmlLayer,_U(slice)); l; ++l)
      layerThickness += xml_comp_t(l).thickness();

    std::cout << "Total Length "    << lcalThickness/dd4hep::cm  << " cm" << std::endl;
    std::cout << "Layer Thickness " << layerThickness/dd4hep::cm << " cm" << std::endl;

    //Loop for repeat=NN
    for(int i=0, repeat=xmlLayer.repeat(); i<repeat; ++i)  {

      std::string layer_name = detName + DD4hep::XML::_toString(thisLayerId,"_layer%d");
      DD4hep::Geometry::Tube layer_base(lcalInnerR,lcalOuterR,layerThickness*0.5);
      
      DD4hep::Geometry::Volume layer_vol(layer_name,layer_base,air);


      int sliceID=1;
      double inThisLayerPosition = -layerThickness*0.5;
      for(DD4hep::XML::Collection_t collSlice(xmlLayer,_U(slice)); collSlice; ++collSlice)  {
	DD4hep::XML::Component compSlice = collSlice;
	const double      sliceThickness = compSlice.thickness();
	const std::string sliceName = layer_name + DD4hep::XML::_toString(sliceID,"slice%d");
	DD4hep::Geometry::Material   slice_mat  = lcdd.material(compSlice.materialStr());

	DD4hep::Geometry::Tube sliceBase(lcalInnerR,lcalOuterR,sliceThickness/2);

	DD4hep::Geometry::Volume slice_vol (sliceName,sliceBase,slice_mat);

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

  const DD4hep::Geometry::Position bcForwardPos (std::tan(-0.5*mradFullCrossingAngle)*lcalCentreZ,0.0, lcalCentreZ);
  const DD4hep::Geometry::Position bcBackwardPos(std::tan(-0.5*mradFullCrossingAngle)*lcalCentreZ,0.0,-lcalCentreZ);
  const DD4hep::Geometry::Rotation3D bcForwardRot ( DD4hep::Geometry::RotationY(-mradFullCrossingAngle*0.5 ) );
  const DD4hep::Geometry::Rotation3D bcBackwardRot( DD4hep::Geometry::RotationZYX ( (M_PI), (M_PI-mradFullCrossingAngle*0.5), (0.0)));

  DD4hep::Geometry::DetElement LumiCals ( detName, xmlLumiCal.id() );
  DD4hep::Geometry::DetElement sdet ( "LumiCal01", xmlLumiCal.id() );
  DD4hep::Geometry::DetElement rdet ( "LumiCal02", xmlLumiCal.id() );

  DD4hep::Geometry::Volume motherVol = lcdd.pickMotherVolume(sdet);

  DD4hep::Geometry::PlacedVolume pv =
    assembly.placeVolume(envelopeVol, DD4hep::Geometry::Transform3D( bcForwardRot, bcForwardPos ) );
  pv.addPhysVolID("system",xmlLumiCal.id());
  pv.addPhysVolID("barrel", 1);
  sdet.setPlacement(pv);

  DD4hep::Geometry::PlacedVolume pv2 =
    assembly.placeVolume(envelopeVol, DD4hep::Geometry::Transform3D( bcBackwardRot, bcBackwardPos ) );
  pv2.addPhysVolID("system",xmlLumiCal.id());
  pv2.addPhysVolID("barrel", 2);
  rdet.setPlacement(pv2);
 
  LumiCals.add(sdet);
  LumiCals.add(rdet);

  pv = motherVol.placeVolume( assembly ) ;
  pv.addPhysVolID("system",xmlLumiCal.id());
  LumiCals.setPlacement( pv ) ;

  return LumiCals;
}

DECLARE_DETELEMENT(LumiCal,create_detector)
