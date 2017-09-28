#include <DD4hep/DetFactoryHelper.h>
#include <XML/Layering.h>

#include <string>

using dd4hep::Assembly;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::Cone;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::IntersectionSolid;
using dd4hep::Layering;
using dd4hep::Layer;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::PolyhedraRegular;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::Rotation3D;
using dd4hep::RotationY;
using dd4hep::RotationZ;
using dd4hep::RotationZYX;
using dd4hep::SensitiveDetector;
using dd4hep::SubtractionSolid;
using dd4hep::Transform3D;
using dd4hep::Translation3D;
using dd4hep::Trapezoid;
using dd4hep::Tube;
using dd4hep::Volume;
using dd4hep::_toString;
using dd4hep::UnionSolid;
using dd4hep::IntersectionSolid;
using dd4hep::Segmentation;


static Ref_t create_detector(Detector& theDetector,
                             xml_h element,
                             SensitiveDetector sens) {

  std::cout << __PRETTY_FUNCTION__  << std::endl;
  std::cout << "Here is my LumiCal"  << std::endl;
  std::cout << " and this is the sensitive detector: " << &sens  << std::endl;
  sens.setType("calorimeter");
  //Materials
  Material air = theDetector.air();

  //Access to the XML File
  xml_det_t     xmlLumiCal    = element;
  const std::string detName = xmlLumiCal.nameStr();

 //--------------------------------
  Assembly assembly( detName + "_assembly"  ) ;
  //--------------------------------

  dd4hep::xml::Dimension dimensions =  xmlLumiCal.dimensions();

  //LumiCal Dimensions
  const double lcalInnerR = dimensions.inner_r();
  const double lcalOuterR = dimensions.outer_r();
  // const double lcalInnerZ = dimensions.inner_z();
  const double lcalThickness = Layering(xmlLumiCal).totalThickness();
  // const double lcalCentreZ = lcalInnerZ+lcalThickness*0.5;
  const double BoltR = 4.*dd4hep::mm ; 
  const double BoltEarThickenss = 2.*dd4hep::mm ; 
  const int NSectors = 12;
  const int NBolts = 3;
  double ear_gap = 10. *dd4hep::mm;
  const double lcalExtra = 20 * dd4hep::cm;

  // counter for the current layer to be placed
  // int thisLayerId = 1;

  //Parameters we have to know about
  dd4hep::xml::Component xmlParameter = xmlLumiCal.child(_Unicode(parameter));
  const double mradFullCrossingAngle  = xmlParameter.attr< double >(_Unicode(crossingangle));
  std::cout << " The crossing angle is: " << mradFullCrossingAngle << " radian"  << std::endl;


  //Envelope to place the layers in
//  Tube envelopeTube (lcalInnerR, 2*lcalOuterR, lcalThickness*0.5 );
Box envelopeTube (5*dd4hep::m, 5*dd4hep::m, 5*dd4hep::m );
  Volume     envelopeVol(detName+"_envelope",envelopeTube,air);
  envelopeVol.setVisAttributes(theDetector,xmlLumiCal.visStr());


  //Start building the support structures 
  Tube BoltTube(0, BoltR, lcalThickness*0.5 );
  Volume BoltVol (detName+"_Bolt",BoltTube,air);
  BoltVol.setVisAttributes(theDetector,"BCLayerVis1");

  Tube BoltEarTube1(0, 2*(BoltR+BoltEarThickenss), lcalThickness*0.5 );
  Tube BoltEarTube2(0, lcalInnerR, lcalThickness*0.5*1.1,-1.0*M_PI/(double)NSectors, M_PI/(double)NSectors );

  Position BoltTransform (-1.0*lcalInnerR, 0.0, 0.0);
  SubtractionSolid BoltEarShape (BoltEarTube1, BoltEarTube2, BoltTransform);

  Volume BoltEarVol (detName+"_BoltEar",BoltEarShape,air);
  BoltEarVol.setVisAttributes(theDetector,"BCLayerVis3");

  // PlacedVolume pvt =
    BoltEarVol.placeVolume(BoltVol,Position(BoltR+BoltEarThickenss, 0.0, 0.0));



  double dPhi  = 2*M_PI/(double)NSectors;
  double suppAng = dPhi/2.0;
   
	for ( int ie=0 ; ie < NSectors; ie++ ){
		
  		const Position EarPos (lcalOuterR * cos (suppAng),-1.0*lcalOuterR * sin (suppAng),0);
  		const RotationZ EarRot( -1.0*suppAng ) ;
		if(ie!=0 && ie!=5 && ie!=6 && ie!=11)	// this line is terrible it has to be corrected	
		// PlacedVolume pvear5 =
		  envelopeVol.placeVolume(BoltEarVol, Transform3D( EarRot, EarPos ) );

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

      Tube FEMothSolid (lcalOuterR, lcalOuterR+lcalExtra, lcalThickness, FE_phi0, FE_dphi);
      Tube FEBoardSolid (lcalOuterR, lcalOuterR+lcalExtra, FE_hth, FE_phi0, FE_dphi);
      Tube FECoolSolid (lcalOuterR, lcalOuterR+lcalExtra, FECool_hdz, FE_phi0, FE_dphi);
      Tube FESectorSolid (lcalOuterR, lcalOuterR+lcalExtra, FEChip_hdz, FE_phi0, FE_Sec_dphi);
      Tube ChipMothSolid (lcalOuterR, lcalOuterR+lcalExtra, FEChip_hdz, FE_phi0, FE_dphi);
    // FE chips
      Box ChipSolid(FEChip_hx, FEChip_hy, FEChip_hdz);

std::cout<<lcalOuterR<<"\t"<<lcalOuterR+lcalExtra<<"\t"<<FEChip_hdz<<"\t"<<FE_phi0<<"\t"<<FE_Sec_dphi<<std::endl;

 // FE chips
      Volume FEMothLog  ("FEMotherLog", FEMothSolid, air ); 
      Volume FECoolLog  ("FECoolLog", FECoolSolid, air ); //Alu);
      Volume FEBoardLog ("FEBoardLog",  FEBoardSolid, air ); //fanele2);
      Volume ChipMothLog ("ChipMothLog", ChipMothSolid, air ); //silicon,);
      Volume FESectorLog ("FESectorLog",  FESectorSolid, air ); 
      Volume FEChipLog ("FEChipLog", ChipSolid, air ); //silicon );

      double xCh = lcalOuterR+lcalExtra/2.;

  	const Position FE1Pos (xCh * cos (FE_phi0 + FE_Sec_dphi/2),-1.0*xCh * sin (FE_phi0 + FE_Sec_dphi/2),0);
  	const RotationZ FE1Rot( FE_phi0 + FE_Sec_dphi/2 ) ;


  FEChipLog.setVisAttributes(theDetector,"BCLayerVis1");
  FESectorLog.setVisAttributes(theDetector,"BCLayerVis3");

    // PlacedVolume pvear6 =
      FESectorLog.placeVolume(FEChipLog, Transform3D( FE1Rot, FE1Pos ) );
    // PlacedVolume pvear7 =
      envelopeVol.placeVolume(FESectorLog);


        double zpos = lcalThickness - FECool_hdz;
        // PlacedVolume pvear8 =
	  FEMothLog.placeVolume(FECoolLog, Position( 0,0,zpos) );
  //      new G4PVPlacement ( 0, G4ThreeVector(0.,0.,zpos), FECoolLog, "LCalFECooling", FEMothLog, false, 1);
	// PCB
	zpos -= (FECool_hdz + FE_hth);
        // PlacedVolume pvear9 =
	  FEMothLog.placeVolume(FEBoardLog, Position( 0,0,zpos) );
        //new G4PVPlacement ( 0, G4ThreeVector(0.,0.,zpos), FEBoardLog, "LCalFEBoard", FEMothLog, false, 1);
	// chips
	zpos -= (FE_hth + FEChip_hdz);
        //new G4PVPlacement ( 0, G4ThreeVector(0.,0.,zpos), ChipMothLog, "LCalFEchip", FEMothLog, false, 1);
        // PlacedVolume pvear10 =
	  FEMothLog.placeVolume(ChipMothLog, Position( 0,0,zpos) );

	// everything in support layer
        zpos = -0.1*dd4hep::mm/2.;
        double phi_rot = 0.; 

        for ( int k=0; k< NBolts ; k++ ){
  		const Position SubPos1 (zpos * cos (phi_rot),-1.0*zpos * sin (phi_rot),0);
  		const RotationZ SubRot1( phi_rot ) ;

  		const Position SubPos2 (zpos * cos (phi_rot+M_PI),-1.0*zpos * sin (phi_rot+M_PI),0);
  		const RotationZ SubRot2( phi_rot+M_PI ) ;

        // PlacedVolume pvear11 =
	  envelopeVol.placeVolume(FEMothLog, Transform3D( SubRot1, SubPos1 ) );
        // PlacedVolume pvear12 =
	  envelopeVol.placeVolume(FEMothLog, Transform3D( SubRot2, SubPos2 ) );
        	phi_rot += dPhi;
	}


//      new G4PVPlacement ( transFE1, FEChipLog, "LcalFEchip", FESectorLog, false, 0); 
     //


  //PlacedVolume pv = assembly.placeVolume(envelopeVol, Transform3D( bcForwardRot2, bcForwardPos2 ) );

/*

  Volume BoltEarVol2 (detName+"_BoltEar2",BoltEarTube1,air);
  BoltEarVol2.setVisAttributes(theDetector,"BCLayerVis1");



  ////////////////////////////////////////////////////////////////////////////////
  // Create all the layers
  ////////////////////////////////////////////////////////////////////////////////

  //Loop over all the layer (repeat=NN) sections
  //This is the starting point to place all layers, we need this when we have more than one layer block
  double referencePosition = -lcalThickness*0.5;
  for(dd4hep::xml::Collection_t coll(xmlLumiCal,_U(layer)); coll; ++coll)  {
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
      Tube layer_base(lcalInnerR,lcalOuterR,layerThickness*0.5);
      
      Volume layer_vol(layer_name,layer_base,air);


      int sliceID=1;
      double inThisLayerPosition = -layerThickness*0.5;
      for(dd4hep::xml::Collection_t collSlice(xmlLayer,_U(slice)); collSlice; ++collSlice)  {
	dd4hep::xml::Component compSlice = collSlice;
	const double      sliceThickness = compSlice.thickness();
	const std::string sliceName = layer_name + dd4hep::xml::_toString(sliceID,"slice%d");
	Material   slice_mat  = theDetector.material(compSlice.materialStr());

	Tube sliceBase(lcalInnerR,lcalOuterR,sliceThickness/2);

	Volume slice_vol (sliceName,sliceBase,slice_mat);

	if ( compSlice.isSensitive() )  {
	  slice_vol.setSensitiveDetector(sens);
	}

	slice_vol.setAttributes(theDetector,compSlice.regionStr(),compSlice.limitsStr(),compSlice.visStr());
	PlacedVolume pv = layer_vol.placeVolume(slice_vol,
								  Position(0,0,inThisLayerPosition+sliceThickness*0.5));
	pv.addPhysVolID("slice",sliceID);
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
*/

  const Position bcForwardPos (0,0,0);
  const Position bcBackwardPos(0,0,0);
  const Rotation3D bcForwardRot ( RotationY(0 ) );
  const Rotation3D bcBackwardRot ( RotationY(0 ) );

  DetElement LumiCals ( detName, xmlLumiCal.id() );
  DetElement sdet ( "LumiCal01", xmlLumiCal.id() );
  DetElement rdet ( "LumiCal02", xmlLumiCal.id() );

  Volume motherVol = theDetector.pickMotherVolume(sdet);

  PlacedVolume pv = assembly.placeVolume(envelopeVol, Transform3D( bcForwardRot, bcForwardPos ) );
  pv.addPhysVolID("system",xmlLumiCal.id());
  pv.addPhysVolID("barrel", 1);
  sdet.setPlacement(pv);
/*
  PlacedVolume pv2 =
    assembly.placeVolume(BoltEarVol, Transform3D( bcBackwardRot, bcBackwardPos ) );
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

DECLARE_DETELEMENT(LumiCal_o1_v02,create_detector)
