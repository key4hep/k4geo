// $Id$
#include <DD4hep/DetFactoryHelper.h>
#include "DD4hep/DetType.h"
#include <XML/Layering.h>
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

#include <string>

using dd4hep::Assembly;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::Cone;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::DetType;
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
using dd4hep::EllipticalTube;

using dd4hep::rec::LayeredCalorimeterData;

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

static Ref_t create_detector(Detector& theDetector,
                             xml_h element,
                             SensitiveDetector sens) {
    
  std::cout << __PRETTY_FUNCTION__  << std::endl << std::endl;
  std::cout << " Here is my LumiCal_o1_v03 :"  << std::endl;
  std::cout << "----------------------------"  << std::endl;
  sens.setType("calorimeter");
  //Materials
  Material air = theDetector.air();
    
    //Access to the XML File
    xml_det_t     xmlLumiCal    = element;
    const std::string detName   = xmlLumiCal.nameStr();
    
    DetElement sdet ( detName, xmlLumiCal.id() );
    
    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector, element , sdet ) ;
    DetElement lumiCalDE_1(sdet,"Calorimeter1",1);
    DetElement lumiCalDE_2(sdet,"Calorimeter2",2);

    sdet.setTypeFlag( DetType::CALORIMETER |  DetType::ENDCAP  | DetType::ELECTROMAGNETIC |  DetType::FORWARD ) ;

    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
    
    //-----------------------------------------------------------------------------------
    
    //Parameters we have to know about
    dd4hep::xml::Component xmlParameter = xmlLumiCal.child(_Unicode(parameter));
    const double fullCrossingAngle  = xmlParameter.attr< double >(_Unicode(crossingangle));

    
    //LumiCal Dimensions
    dd4hep::xml::Dimension dimensions =  xmlLumiCal.dimensions();
    const double lcalInnerR  = dimensions.inner_r();
    const double lcalOuterR  = dimensions.outer_r();
    const double lcalInnerZ  = dimensions.inner_z();
    const double layerStagger= dimensions.attr< double >("layerstagger");
    const double serviceRmax = dimensions.attr< double >("servicesRmax");
    const double lcalExtraS  = serviceRmax - lcalOuterR;
    const double lcalSectors = dimensions.attr<double>("phi_sectors");
    const double lcalRSects  = dimensions.attr<double>("r_sectors");
    const double lcalModules = dimensions.attr<double>("modules");
    const double lcalRGap    = dimensions.attr<double>("r_gap");
    const double lcalPhiGap  = dimensions.attr<double>("phi_gap");
    const double lcalPhi0    = dimensions.attr<double>("outerPhi0");

    const double lcalTileDphi  = (360./lcalModules)*dd4hep::deg;
    const double cellPhiSize   = lcalTileDphi/lcalSectors;
    // FIXME: staggerPhi must be somehow passed to DDrec
    // temporary it goes in place of obsolete layer thickness
    const double staggerPhi    = layerStagger*cellPhiSize;   
    const double lcalThickness = Layering(xmlLumiCal).totalThickness();
    const double lcalCentreZ   = lcalInnerZ+lcalThickness*0.5;
    const double lcalXoffset = lcalCentreZ * std::tan( fullCrossingAngle/2. );

    // inner/outer radii are not the sensor dims, these we have to compute
    const double sensInnerR = lcalInnerR/cos( lcalSectors*cellPhiSize/2. ) + lcalRGap;
    const double sensOuterR = lcalOuterR*cos( lcalSectors*cellPhiSize/2. ) - lcalRGap ;
    const double cellRsize  = (sensOuterR - sensInnerR )/lcalRSects;

    //========== fill data for reconstruction ============================
    LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
    caloData->layoutType = LayeredCalorimeterData::EndcapLayout ;
    caloData->inner_symmetry = 0  ; // hardcoded tube
    caloData->outer_symmetry = 0  ; 
     
    /// extent of the calorimeter sensor in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
    caloData->extent[0] = sensInnerR; //lcalInnerR ;
    caloData->extent[1] = sensOuterR; //lcalOuterR ;
    caloData->extent[2] = lcalInnerZ ;
    caloData->extent[3] = lcalInnerZ + lcalThickness ;
    caloData->outer_phi0 = lcalPhi0;
    // this is rather FYIO, not really needed for reco
    caloData->gap0 = lcalPhiGap;    
    caloData->gap2 = lcalRGap;    

    std::cout << " The crossing angle is:  " << fullCrossingAngle << " radian"  << std::endl;
    std::cout << "     "<< detName << " z-begin : " << lcalInnerZ/dd4hep::mm << " mm" << std::endl;
    std::cout << "     "<< detName << " z-end   : " << (lcalInnerZ+lcalThickness)/dd4hep::mm << " mm" << std::endl;
    std::cout << "   (x,y,z)-center : " << "( " << lcalXoffset/dd4hep::mm << ",0 ,+- " << lcalCentreZ/dd4hep::mm << " ) mm" << std::endl;
    std::cout << "     sensorInnerR : " << sensInnerR/dd4hep::mm << " mm" << std::endl;
    std::cout << "     sensorouterR : " << sensOuterR/dd4hep::mm << " mm" << std::endl;
    std::cout << "        thickness : " << lcalThickness/dd4hep::mm << " mm" << std::endl;
    std::cout << "          cell_dR : " << cellRsize/dd4hep::mm <<  " mm" << std::endl;
    std::cout << "        cell_dPhi : " << cellPhiSize << "  radian" << std::endl;
    std::cout << "       offset Phi : " << lcalPhi0/dd4hep::deg   << " deg"  << std::endl;
    std::cout << "    stagger angle : " << staggerPhi/dd4hep::deg   << " deg"  << std::endl;
    std::cout << "    tile gap size : " << 2.*lcalRGap/dd4hep::mm << " mm" << std::endl;
    std::cout << "    tile span phi : " << lcalTileDphi/dd4hep::deg << " deg" << std::endl;
    std::cout << "   services space : " << lcalExtraS/dd4hep::mm << " mm" << std::endl;
    
    
    //Whole Lcal Envelope to place the layers in
    Tube envelopeTube (lcalInnerR, lcalOuterR + lcalExtraS, lcalThickness*0.5 );
    Volume     envelopeVol(detName+"_module",envelopeTube,air);
    envelopeVol.setVisAttributes(theDetector,xmlLumiCal.visStr());
    //
    // build services space containing support and electronics
    // simplified NOT layered for compactness
    double bolt_radius = 0.6 *dd4hep::cm; 
    double ear_height = lcalExtraS - 0.1*dd4hep::cm;
    double ear_dy    = lcalOuterR + ear_height;
   if( lcalExtraS > 0. ){
      Tube envelopeExtras (lcalOuterR, lcalOuterR + lcalExtraS-0.1, lcalThickness*0.5 );
      Volume     ServicesVol(detName+"_Services",envelopeExtras,air);
      envelopeVol.placeVolume( ServicesVol );
      // mounting ears
      Material earMaterial = theDetector.material("TungstenDens24");
      int n_ears = 3;  // ears come in pairs
      double ear_width = ear_height/2.;
      double ear_dx    = ear_width/sqrt(1. - ( lcalOuterR*lcalOuterR - ear_width*ear_width)/(ear_dy*ear_dy) );
      double ear_hdz   = lcalThickness*0.5;
      double earsDPhi = 180./double( n_ears ) * dd4hep::deg;
      // anonymous volumes - support
      EllipticalTube EarEllipse( ear_dx, ear_dy, ear_hdz );
      Tube earClip ( 0., lcalOuterR,  (ear_hdz + 0.1*dd4hep::cm), 0., 360.*dd4hep::deg );
      Tube boltHole( 0., bolt_radius, (ear_hdz + 0.1*dd4hep::cm) );
      // mixed elcetronics
      double phiSpanFE = earsDPhi - 2.* atan( ear_width/sqrt(lcalOuterR*lcalOuterR-ear_width*ear_width) ); 
      Material mixFEmat = theDetector.material("siPCBMix");
      Tube mixFE( lcalOuterR, lcalOuterR+lcalExtraS-0.2*dd4hep::cm, 0.5*lcalThickness, -phiSpanFE/2., phiSpanFE/2.);
      Volume mixFEvol(detName+"FEmix", mixFE, mixFEmat );
      // clipping
      // pair of ears
      SubtractionSolid EarShape0( EarEllipse, earClip );
      // holes for bolts 
      Position boltPos1( 0.,  (lcalOuterR+3.*bolt_radius), 0.);
      Position boltPos2( 0., -(lcalOuterR+3.*bolt_radius), 0.);
      SubtractionSolid EarShape1( EarShape0, boltHole, boltPos1);
      SubtractionSolid EarShape2( EarShape1, boltHole, boltPos2);
      // final pair of ears
      Volume EarsVolume(detName+"Ears", EarShape2, earMaterial);
      EarsVolume.setVisAttributes( theDetector,"GrayVis");
      // place mountings and electronics  in Services volume
      double earsAng = 0.;
      double feAng = 0.;
      for( int ie=0; ie<n_ears; ie++ ){
	RotationZ earsRotZ( earsAng );
	RotationZ feRotZ1( feAng );
	RotationZ feRotZ2( feAng + M_PI );
	Position  earsPos( 0., 0., 0.);
	ServicesVol.placeVolume( EarsVolume, Transform3D( earsRotZ, earsPos ));
	ServicesVol.placeVolume( mixFEvol, Transform3D( feRotZ1, earsPos) );
	ServicesVol.placeVolume( mixFEvol, Transform3D( feRotZ2, earsPos) );
	earsAng += earsDPhi;
	feAng += earsDPhi;
      }


   }// end if (  lcalExtraS > 0. )
  
    // need this to make  tile gaps  
	double gapDy = lcalPhiGap/2.;
	double cs = sqrt( 1. - (gapDy/sensOuterR)*(gapDy/sensOuterR));
	double gapDx   = ( sensOuterR*cs - sensInnerR )*0.5; // must not stick out , so a bit shorter
	double gapPosX = sensInnerR + gapDx;  
	int nGaps = (int)lcalModules;
    ////////////////////////////////////////////////////////////////////////////////
    // Create all the layers
    ////////////////////////////////////////////////////////////////////////////////
    
    //  Loop over all the layer (repeat=NN) sections
    //  counter for the current layer to be placed
    int thisLayerId = 0;
    double mtotalRadLen = 0.;
    double mtotalDepthZ = 0.;
   //This is the starting point to place all layers, we need this when we have more than one layer block
    double referencePosition = -lcalThickness*0.5;
    for(dd4hep::xml::Collection_t coll(xmlLumiCal,_U(layer)); coll; ++coll)  {
        dd4hep::xml::Component xmlLayer(coll); //we know this thing is a layer
        
        double layerThickness = 0;
        for(dd4hep::xml::Collection_t l(xmlLayer,_U(slice)); l; ++l) layerThickness += xml_comp_t(l).thickness();
        
        //Loop for repeat=NN
        for(int i=0, repeat=xmlLayer.repeat(); i<repeat; ++i)  {
            
            std::string layer_name = detName + dd4hep::xml::_toString(thisLayerId,"_layer%d");
            Tube layer_base(lcalInnerR,lcalOuterR,layerThickness*0.5, lcalPhi0, lcalPhi0+360.*dd4hep::deg);            
            Volume layer_vol(layer_name,layer_base,air);
                       
            int sliceID=0;
            double inThisLayerPosition = -layerThickness*0.5;
            double nRadiationLengths=0.;
            double nInteractionLengths=0.;
            double thickness_sum=0;
            
            LayeredCalorimeterData::Layer caloLayer ;
            
            for( dd4hep::xml::Collection_t collSlice(xmlLayer,_U(slice)); collSlice; ++collSlice )  {
	      dd4hep::xml::Component compSlice = collSlice;
	      double slice_thickness = compSlice.thickness();
	      double slice_minR = lcalInnerR;
	      double slice_maxR = lcalOuterR;
	      // sensor dims differ
	      if ( compSlice.isSensitive() ){
		slice_minR = sensInnerR;
		slice_maxR = sensOuterR;
	      }

	      std::string sliceName = layer_name + dd4hep::xml::_toString(sliceID,"slice%d");
	      Material   slice_material  = theDetector.material(compSlice.materialStr());
                
	      Tube sliceBase(slice_minR,slice_maxR,slice_thickness/2);                
	      Volume slice_vol (sliceName,sliceBase,slice_material);
                
                nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
                nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
                thickness_sum += slice_thickness/2;
              
                if ( compSlice.isSensitive() )  {
		  if ( gapDy > 0. ) {      // put gaps in sensor slice
		    Box gapEnv (gapDx, gapDy , slice_thickness/2. );
		    Volume gapVol ("phiGap", gapEnv, slice_material );
		    gapVol.setVisAttributes(theDetector,"PhiGapVis");
		    double rotPhi = lcalPhi0 - cellPhiSize/2.;
		    for( int gap=0; gap < nGaps; gap++ ) {
		      RotationZ gapRot ( rotPhi );
		      Position gapPos ( gapPosX*cos( rotPhi ), gapPosX*sin( rotPhi ), 0.0 );
		      slice_vol.placeVolume( gapVol,Transform3D( gapRot, gapPos ));
		      rotPhi += lcalTileDphi;
		    }
		  }

		   mtotalRadLen += nRadiationLengths; 
		   mtotalDepthZ += thickness_sum; 
 #if DD4HEP_VERSION_GE( 0, 15 )
                   //Store "inner" quantities
                    caloLayer.inner_nRadiationLengths = nRadiationLengths;
                    caloLayer.inner_nInteractionLengths = nInteractionLengths;
                    caloLayer.inner_thickness = thickness_sum;
                    //Store scintillator thickness
                    caloLayer.sensitive_thickness = slice_thickness;
                    //Reset counters to measure "outside" quantitites
                    nRadiationLengths=0.;
                    nInteractionLengths=0.;
                    thickness_sum = 0.;
#endif                    
                    slice_vol.setSensitiveDetector(sens);
                }
                
                nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
                nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
                thickness_sum += slice_thickness/2;
                
                slice_vol.setAttributes(theDetector,compSlice.regionStr(),compSlice.limitsStr(),compSlice.visStr());

		Position slicePos(0,0,inThisLayerPosition+slice_thickness*0.5);
                layer_vol.placeVolume(slice_vol, slicePos);
                    
                inThisLayerPosition += slice_thickness;
                ++sliceID;
            }//For all slices in this layer

	    mtotalRadLen += nRadiationLengths; 
	    mtotalDepthZ  += thickness_sum; 
            
            ///Needs to be innermost face distance
            caloLayer.distance = lcalCentreZ + referencePosition;

#if DD4HEP_VERSION_GE( 0, 15 )
            caloLayer.outer_nRadiationLengths = nRadiationLengths;
            caloLayer.outer_nInteractionLengths = nInteractionLengths;
            caloLayer.outer_thickness = thickness_sum;
            if ( thisLayerId == 0 )
	      std::cout<<"  Layer thickness : "
		       << (caloLayer.inner_thickness+caloLayer.outer_thickness)/dd4hep::mm 
		       << "   ( rad. length X0:  "
		       << (caloLayer.outer_nRadiationLengths+ caloLayer.inner_nRadiationLengths)
		       << " )" << std::endl;
#else
            if ( thisLayerId == 0 )
	      std::cout<<"  Layer thickness : " <<  mtotalDepthZ/dd4hep::mm 
		       << "   ( rad. length X0:  " <<  mtotalRadLen  << " )" << std::endl;
#endif            
            caloLayer.cellSize0 = cellRsize ;
            caloLayer.cellSize1 = cellPhiSize ;

            Position  layerPos(0,0,referencePosition+0.5*layerThickness);
	    RotationZ layerRot( staggerPhi );
	    PlacedVolume pv;
	    // every other layer has gaps staggered 
	    if( thisLayerId%2 == 1 && layerStagger != 0. ) {
	      pv = envelopeVol.placeVolume(layer_vol, Transform3D( layerRot, layerPos ));
	      caloLayer.phi0 = staggerPhi;
	    }else{
	      pv = envelopeVol.placeVolume(layer_vol,layerPos);
	      caloLayer.phi0 = 0.;
	    }
	    pv.addPhysVolID("layer",thisLayerId);
              
            caloData->layers.push_back( caloLayer ) ;
            layer_vol.setVisAttributes(theDetector,xmlLayer.visStr());
             
            
            referencePosition += layerThickness;
            ++thisLayerId;
            
        }//for all layers
        
    }// for all layer collections
    std::cout<<"   LumiCal length : "<< mtotalDepthZ/dd4hep::mm 
	     << " ( rad. length X0: "<< mtotalRadLen << "   )"<<std::endl;
    std::cout << "-----------------------------------------------------------------"  << std::endl<<std::endl;;
    
    const Position bcForwardPos ( lcalXoffset,0.0, lcalCentreZ);
    const Position bcBackwardPos( lcalXoffset,0.0,-lcalCentreZ);
    const RotationY bcForwardRot ( fullCrossingAngle*0.5  );
    const RotationY bcBackwardRot( M_PI-fullCrossingAngle*0.5 );
    
    PlacedVolume pv =
    envelope.placeVolume(envelopeVol, Transform3D( bcForwardRot, bcForwardPos ) );
    pv.addPhysVolID("barrel", 1);
    lumiCalDE_1.setPlacement(pv);

    PlacedVolume pv2 =
    envelope.placeVolume(envelopeVol, Transform3D( bcBackwardRot, bcBackwardPos ) );
    pv2.addPhysVolID("barrel", 2);
    lumiCalDE_2.setPlacement(pv2);
    
    sdet.addExtension< LayeredCalorimeterData >( caloData ) ;
    
    return sdet;
}
                                               
DECLARE_DETELEMENT(LumiCal_o1_v03,create_detector)
