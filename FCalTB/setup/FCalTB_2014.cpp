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
  const double BoltR = 4.*dd4hep::mm ; 
  const double BoltEarThickenss = 2.*dd4hep::mm ; 
  const int NSectors = 12;
  const int NBolts = 3;
  double ear_gap = 10. *dd4hep::mm;
  const double lcalExtra = 20 * dd4hep::cm;

  // counter for the current layer to be placed
  int thisLayerId = 1;

  //Parameters we have to know about
  DD4hep::XML::Component xmlParameter = xmlLumiCal.child(_Unicode(parameter));
  const double mradFullCrossingAngle  = xmlParameter.attr< double >(_Unicode(crossingangle));
  std::cout << " The crossing angle is: " << mradFullCrossingAngle << " radian"  << std::endl;


  //Envelope to place the layers in
//  DD4hep::Geometry::Tube envelopeTube (lcalInnerR, 2*lcalOuterR, lcalThickness*0.5 );
DD4hep::Geometry::Box envelopeTube (5*dd4hep::m, 5*dd4hep::m, 5*dd4hep::m );
  DD4hep::Geometry::Volume     envelopeVol(detName+"_envelope",envelopeTube,air);
  envelopeVol.setVisAttributes(lcdd,xmlLumiCal.visStr());


  //Start building the support structures 
  DD4hep::Geometry::Tube BoltTube(0, BoltR, lcalThickness*0.5 );
  DD4hep::Geometry::Volume BoltVol (detName+"_Bolt",BoltTube,air);
  BoltVol.setVisAttributes(lcdd,"BCLayerVis1");

  DD4hep::Geometry::Tube BoltEarTube1(0, 2*(BoltR+BoltEarThickenss), lcalThickness*0.5 );
  DD4hep::Geometry::Tube BoltEarTube2(0, lcalInnerR, lcalThickness*0.5*1.1,-1.0*M_PI/(double)NSectors, M_PI/(double)NSectors );

  DD4hep::Geometry::Position BoltTransform (-1.0*lcalInnerR, 0.0, 0.0);
  DD4hep::Geometry::SubtractionSolid BoltEarShape (BoltEarTube1, BoltEarTube2, BoltTransform);

  DD4hep::Geometry::Volume BoltEarVol (detName+"_BoltEar",BoltEarShape,air);
  BoltEarVol.setVisAttributes(lcdd,"BCLayerVis3");

  DD4hep::Geometry::PlacedVolume pvt = BoltEarVol.placeVolume(BoltVol,DD4hep::Geometry::Position(BoltR+BoltEarThickenss, 0.0, 0.0));



  double dPhi  = 2*M_PI/(double)NSectors;
  double suppAng = dPhi/2.0;
   
	for ( int ie=0 ; ie < NSectors; ie++ ){
		
  		const DD4hep::Geometry::Position EarPos (lcalOuterR * cos (suppAng),-1.0*lcalOuterR * sin (suppAng),0);
  		const DD4hep::Geometry::RotationZ EarRot( -1.0*suppAng ) ;
		if(ie!=0 && ie!=5 && ie!=6 && ie!=11)	// this line is terrible it has to be corrected	
		DD4hep::Geometry::PlacedVolume pvear5 = envelopeVol.placeVolume(BoltEarVol, DD4hep::Geometry::Transform3D( EarRot, EarPos ) );

		suppAng += dPhi;

		
      }




  // FE mother board and cooling
      double FE_hth = 0.5 *dd4hep::mm;
      double FE_phi0 = atan( ear_gap/lcalOuterR );
      double FE_dphi = dPhi - 2.*FE_phi0;
   // FE chips
      double FEChip_hx   = 0.4*lcalExtra;
      int    nFE_Sectors  = 64/(2*NBolts);       
      double FE_Sec_dphi  =  FE_dphi/double( nFE_Sectors );
      double FEChip_hy   = (BoltEarThickenss+BoltR)*sin( FE_Sec_dphi/2. )-2.;
      double FEChip_hdz  = 0.25 *dd4hep::mm;
      double FECool_hdz  = 1.00 *dd4hep::mm;

      DD4hep::Geometry::Tube FEMothSolid (lcalOuterR, lcalOuterR+lcalExtra, lcalThickness, FE_phi0, FE_dphi);
      DD4hep::Geometry::Tube FEBoardSolid (lcalOuterR, lcalOuterR+lcalExtra, FE_hth, FE_phi0, FE_dphi);
      DD4hep::Geometry::Tube FECoolSolid (lcalOuterR, lcalOuterR+lcalExtra, FECool_hdz, FE_phi0, FE_dphi);
      DD4hep::Geometry::Tube FESectorSolid (lcalOuterR, lcalOuterR+lcalExtra, FEChip_hdz, FE_phi0, FE_Sec_dphi);
      DD4hep::Geometry::Tube ChipMothSolid (lcalOuterR, lcalOuterR+lcalExtra, FEChip_hdz, FE_phi0, FE_dphi);
    // FE chips
      DD4hep::Geometry::Box ChipSolid(FEChip_hx, FEChip_hy, FEChip_hdz);

std::cout<<lcalOuterR<<"\t"<<lcalOuterR+lcalExtra<<"\t"<<FEChip_hdz<<"\t"<<FE_phi0<<"\t"<<FE_Sec_dphi<<std::endl;

 // FE chips
      DD4hep::Geometry::Volume FEMothLog  ("FEMotherLog", FEMothSolid, air ); 
      DD4hep::Geometry::Volume FECoolLog  ("FECoolLog", FECoolSolid, air ); //Alu);
      DD4hep::Geometry::Volume FEBoardLog ("FEBoardLog",  FEBoardSolid, air ); //fanele2);
      DD4hep::Geometry::Volume ChipMothLog ("ChipMothLog", ChipMothSolid, air ); //silicon,);
      DD4hep::Geometry::Volume FESectorLog ("FESectorLog",  FESectorSolid, air ); 
      DD4hep::Geometry::Volume FEChipLog ("FEChipLog", ChipSolid, air ); //silicon );

      double xCh = lcalOuterR+lcalExtra/2.;

  	const DD4hep::Geometry::Position FE1Pos (xCh * cos (FE_phi0 + FE_Sec_dphi/2),-1.0*xCh * sin (FE_phi0 + FE_Sec_dphi/2),0);
  	const DD4hep::Geometry::RotationZ FE1Rot( FE_phi0 + FE_Sec_dphi/2 ) ;


  FEChipLog.setVisAttributes(lcdd,"BCLayerVis1");
  FESectorLog.setVisAttributes(lcdd,"BCLayerVis3");

    DD4hep::Geometry::PlacedVolume pvear6 = FESectorLog.placeVolume(FEChipLog, DD4hep::Geometry::Transform3D( FE1Rot, FE1Pos ) );
    DD4hep::Geometry::PlacedVolume pvear7 = envelopeVol.placeVolume(FESectorLog);


        double zpos = lcalThickness - FECool_hdz;
        DD4hep::Geometry::PlacedVolume pvear8 = FEMothLog.placeVolume(FECoolLog, DD4hep::Geometry::Position( 0,0,zpos) );
  //      new G4PVPlacement ( 0, G4ThreeVector(0.,0.,zpos), FECoolLog, "LCalFECooling", FEMothLog, false, 1);
	// PCB
	zpos -= (FECool_hdz + FE_hth);
        DD4hep::Geometry::PlacedVolume pvear9 = FEMothLog.placeVolume(FEBoardLog, DD4hep::Geometry::Position( 0,0,zpos) );
        //new G4PVPlacement ( 0, G4ThreeVector(0.,0.,zpos), FEBoardLog, "LCalFEBoard", FEMothLog, false, 1);
	// chips
	zpos -= (FE_hth + FEChip_hdz);
        //new G4PVPlacement ( 0, G4ThreeVector(0.,0.,zpos), ChipMothLog, "LCalFEchip", FEMothLog, false, 1);
        DD4hep::Geometry::PlacedVolume pvear10 = FEMothLog.placeVolume(ChipMothLog, DD4hep::Geometry::Position( 0,0,zpos) );

	// everything in support layer
        zpos = -0.1*dd4hep::mm/2.;
        double phi_rot = 0.; 

        for ( int k=0; k< NBolts ; k++ ){
  		const DD4hep::Geometry::Position SubPos1 (zpos * cos (phi_rot),-1.0*zpos * sin (phi_rot),0);
  		const DD4hep::Geometry::RotationZ SubRot1( phi_rot ) ;

  		const DD4hep::Geometry::Position SubPos2 (zpos * cos (phi_rot+M_PI),-1.0*zpos * sin (phi_rot+M_PI),0);
  		const DD4hep::Geometry::RotationZ SubRot2( phi_rot+M_PI ) ;

        DD4hep::Geometry::PlacedVolume pvear11 = envelopeVol.placeVolume(FEMothLog, DD4hep::Geometry::Transform3D( SubRot1, SubPos1 ) );
        DD4hep::Geometry::PlacedVolume pvear12 = envelopeVol.placeVolume(FEMothLog, DD4hep::Geometry::Transform3D( SubRot2, SubPos2 ) );
        	phi_rot += dPhi;
	}


//      new G4PVPlacement ( transFE1, FEChipLog, "LcalFEchip", FESectorLog, false, 0); 
     //


  //DD4hep::Geometry::PlacedVolume pv = assembly.placeVolume(envelopeVol, DD4hep::Geometry::Transform3D( bcForwardRot2, bcForwardPos2 ) );

/*

  DD4hep::Geometry::Volume BoltEarVol2 (detName+"_BoltEar2",BoltEarTube1,air);
  BoltEarVol2.setVisAttributes(lcdd,"BCLayerVis1");



  ////////////////////////////////////////////////////////////////////////////////
  // Create all the layers
  ////////////////////////////////////////////////////////////////////////////////
/*
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
*/

  const DD4hep::Geometry::Position bcForwardPos (0,0,0);
  const DD4hep::Geometry::Position bcBackwardPos(0,0,0);
  const DD4hep::Geometry::Rotation3D bcForwardRot ( DD4hep::Geometry::RotationY(0 ) );
  const DD4hep::Geometry::Rotation3D bcBackwardRot ( DD4hep::Geometry::RotationY(0 ) );

  DD4hep::Geometry::DetElement LumiCals ( detName, xmlLumiCal.id() );
  DD4hep::Geometry::DetElement sdet ( "LumiCal01", xmlLumiCal.id() );
  DD4hep::Geometry::DetElement rdet ( "LumiCal02", xmlLumiCal.id() );

  DD4hep::Geometry::Volume motherVol = lcdd.pickMotherVolume(sdet);

  DD4hep::Geometry::PlacedVolume pv = assembly.placeVolume(envelopeVol, DD4hep::Geometry::Transform3D( bcForwardRot, bcForwardPos ) );
  pv.addPhysVolID("system",xmlLumiCal.id());
  pv.addPhysVolID("barrel", 1);
  sdet.setPlacement(pv);
/*
  DD4hep::Geometry::PlacedVolume pv2 =
    assembly.placeVolume(BoltEarVol, DD4hep::Geometry::Transform3D( bcBackwardRot, bcBackwardPos ) );
  pv2.addPhysVolID("system",xmlLumiCal.id());
  pv2.addPhysVolID("barrel", 2);
  rdet.setPlacement(pv2);
 */
  LumiCals.add(sdet);
//  LumiCals.add(rdet);

  pv = motherVol.placeVolume( assembly ) ;
  pv.addPhysVolID("system",xmlLumiCal.id());
  LumiCals.setPlacement( pv ) ;

  return LumiCals;
}

DECLARE_DETELEMENT(FCalTB_2014,create_detector)
