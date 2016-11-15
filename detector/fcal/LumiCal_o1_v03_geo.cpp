// $Id$
#include <DD4hep/DetFactoryHelper.h>
#include "DD4hep/DetType.h"
#include <XML/Layering.h>
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

#include <string>

using namespace DD4hep;
using namespace DD4hep::Geometry;

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

static DD4hep::Geometry::Ref_t create_detector(DD4hep::Geometry::LCDD& lcdd,
                                               xml_h element,
                                               DD4hep::Geometry::SensitiveDetector sens) {
    
  std::cout << __PRETTY_FUNCTION__  << std::endl << std::endl;
  std::cout << " Here is my LumiCal_o1_v03 :"  << std::endl;
  std::cout << "----------------------------"  << std::endl;
  sens.setType("calorimeter");
  //Materials
  DD4hep::Geometry::Material air = lcdd.air();
    
    //Access to the XML File
    xml_det_t     xmlLumiCal    = element;
    const std::string detName   = xmlLumiCal.nameStr();
    
    DD4hep::Geometry::DetElement sdet ( detName, xmlLumiCal.id() );
    
    // --- create an envelope volume and position it into the world ---------------------
    
    DD4hep::Geometry::Volume envelope = DD4hep::XML::createPlacedEnvelope( lcdd, element , sdet ) ;
    DD4hep::Geometry::DetElement lumiCalDE_1(sdet,"Calorimeter1",1);
    DD4hep::Geometry::DetElement lumiCalDE_2(sdet,"Calorimeter2",2);

    sdet.setTypeFlag( DetType::CALORIMETER |  DetType::ENDCAP  | DetType::ELECTROMAGNETIC |  DetType::FORWARD ) ;

    if( lcdd.buildType() == DD4hep::BUILD_ENVELOPE ) return sdet ;
    
    //-----------------------------------------------------------------------------------
    
    //Parameters we have to know about
    DD4hep::XML::Component xmlParameter = xmlLumiCal.child(_Unicode(parameter));
    const double fullCrossingAngle  = xmlParameter.attr< double >(_Unicode(crossingangle));

    
    //LumiCal Dimensions
    DD4hep::XML::Dimension dimensions =  xmlLumiCal.dimensions();
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
    const double lcalThickness = DD4hep::Layering(xmlLumiCal).totalThickness();
    const double lcalCentreZ   = lcalInnerZ+lcalThickness*0.5;

    // inner/outer radii are not the sensor dims, these we have to compute
    const double sensInnerR = lcalInnerR/cos( lcalSectors*cellPhiSize/2. ) + lcalRGap;
    const double sensOuterR = lcalOuterR*cos( lcalSectors*cellPhiSize/2. ) - lcalRGap ;
    const double cellRsize  = (sensOuterR - sensInnerR )/lcalRSects;

    //========== fill data for reconstruction ============================
    DD4hep::DDRec::LayeredCalorimeterData* caloData = new DD4hep::DDRec::LayeredCalorimeterData ;
    caloData->layoutType = DD4hep::DDRec::LayeredCalorimeterData::EndcapLayout ;
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
    std::cout << "     sensorInnerR : " << sensInnerR/dd4hep::mm << " mm" << std::endl;
    std::cout << "     sensorouterR : " << sensOuterR/dd4hep::mm << " mm" << std::endl;
    std::cout << "        thickness : " << lcalThickness/dd4hep::mm << " mm" << std::endl;
    std::cout << "          cell_dR : " << cellRsize/dd4hep::mm <<  " mm" << std::endl;
    std::cout << "        cell_dPhi : " << cellPhiSize << "  radian" << std::endl;
    std::cout << "    stagger angle : " << staggerPhi/dd4hep::deg   << " deg"  << std::endl;
    std::cout << "    tile gap size : " << 2.*lcalRGap/dd4hep::mm << " mm" << std::endl;
    std::cout << "    tile span phi : " << lcalTileDphi/dd4hep::deg << " deg" << std::endl;
    std::cout << "   services space : " << lcalExtraS/dd4hep::mm << " mm" << std::endl;
    
    
    //Whole Lcal Envelope to place the layers in
    DD4hep::Geometry::Tube envelopeTube (lcalInnerR, lcalOuterR + lcalExtraS, lcalThickness*0.5 );
    DD4hep::Geometry::Volume     envelopeVol(detName+"_module",envelopeTube,air);
    envelopeVol.setVisAttributes(lcdd,xmlLumiCal.visStr());
    //
    // build services space containing support and electronics
    // simplified NOT layered for compactness
    double bolt_radius = 0.4 *dd4hep::cm; 
    double ear_hight = 10.*bolt_radius;
    double ear_dx  = lcalOuterR + ear_hight;
   if( ear_dx <= serviceRmax ){
      DD4hep::Geometry::Tube envelopeExtras (lcalOuterR, lcalOuterR + lcalExtraS-0.1, lcalThickness*0.5 );
      DD4hep::Geometry::Volume     ServicesVol(detName+"_Services",envelopeExtras,air);
      envelopeVol.placeVolume( ServicesVol );
      // mounting ears
      DD4hep::Geometry::Material earMaterial = lcdd.material("TungstenDens24");
      int n_ears = 3;  // ears come in pairs
      double earsAng0 = 30.*dd4hep::deg;
      double ear_width = ear_hight;
      double ear_dy    = ear_width*ear_dx/sqrt(ear_dx*ear_dx-lcalOuterR*lcalOuterR);
      double ear_hdz   = lcalThickness*0.5;
      // anonymous volumes - support
      DD4hep::Geometry::EllipticalTube EarEllipse( ear_dx, ear_dy, ear_hdz );
      DD4hep::Geometry::Tube earClip ( 0., lcalOuterR,  (ear_hdz + 0.1*dd4hep::cm), 0., 360.*dd4hep::deg );
      DD4hep::Geometry::Tube boltHole( 0., bolt_radius, (ear_hdz + 0.1*dd4hep::cm) );
      // mixed elcetronics
      double phiSpanFE = lcalTileDphi - atan( ear_width/lcalOuterR); 
      DD4hep::Geometry::Material mixFEmat = lcdd.material("siPCBMix");
      DD4hep::Geometry::Tube mixFE( lcalOuterR, lcalOuterR+lcalExtraS-0.2*dd4hep::cm, 0.5*lcalThickness, -phiSpanFE, phiSpanFE);
      DD4hep::Geometry::Volume mixFEvol(detName+"FEmix", mixFE, mixFEmat );
      // clipping
      // pair of ears
      DD4hep::Geometry::SubtractionSolid EarShape0( EarEllipse, earClip );
      // holes for bolts 
      DD4hep::Geometry::Position boltPos1(  (lcalOuterR+3.*bolt_radius), 0., 0.);
      DD4hep::Geometry::Position boltPos2( -(lcalOuterR+3.*bolt_radius), 0., 0.);
      DD4hep::Geometry::SubtractionSolid EarShape1( EarShape0, boltHole, boltPos1);
      DD4hep::Geometry::SubtractionSolid EarShape2( EarShape1, boltHole, boltPos2);
      // final pair of ears
      DD4hep::Geometry::Volume EarsVolume(detName+"Ears", EarShape2, earMaterial);
      EarsVolume.setVisAttributes( lcdd,"GrayVis");
      // place mountings and electronics  in Services volume
      double earsAng = earsAng0;
      double earsDPhi = 180./double( n_ears ) * dd4hep::deg;
      double feAng = 0.;
      for( int ie=0; ie<n_ears; ie++ ){
	DD4hep::Geometry::RotationZ earsRotZ( earsAng );
	DD4hep::Geometry::RotationZ feRotZ1( feAng );
	DD4hep::Geometry::RotationZ feRotZ2( feAng + 180.*dd4hep::deg );
	DD4hep::Geometry::Position  earsPos( 0., 0., 0.);
	ServicesVol.placeVolume( EarsVolume, DD4hep::Geometry::Transform3D( earsRotZ, earsPos ));
	ServicesVol.placeVolume( mixFEvol, DD4hep::Geometry::Transform3D( feRotZ1, earsPos) );
	ServicesVol.placeVolume( mixFEvol, DD4hep::Geometry::Transform3D( feRotZ2, earsPos) );
	earsAng += earsDPhi;
	feAng += earsDPhi;
      }


    }
    // needed to make  tile gaps  
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
    for(DD4hep::XML::Collection_t coll(xmlLumiCal,_U(layer)); coll; ++coll)  {
        DD4hep::XML::Component xmlLayer(coll); //we know this thing is a layer
        
        
        //This just calculates the total size of a single layer
        //Why no convenience function for this?
        double layerThickness = 0;
        for(DD4hep::XML::Collection_t l(xmlLayer,_U(slice)); l; ++l) layerThickness += xml_comp_t(l).thickness();
        
        //Loop for repeat=NN
        for(int i=0, repeat=xmlLayer.repeat(); i<repeat; ++i)  {
            
            std::string layer_name = detName + DD4hep::XML::_toString(thisLayerId,"_layer%d");
            DD4hep::Geometry::Tube layer_base(lcalInnerR,lcalOuterR,layerThickness*0.5, lcalPhi0, lcalPhi0+360.*dd4hep::deg);
            
            DD4hep::Geometry::Volume layer_vol(layer_name,layer_base,air);
            
            
            int sliceID=0;
            double inThisLayerPosition = -layerThickness*0.5;
            double nRadiationLengths=0.;
            double nInteractionLengths=0.;
            double thickness_sum=0;
            
            DD4hep::DDRec::LayeredCalorimeterData::Layer caloLayer ;
            
            for( DD4hep::XML::Collection_t collSlice(xmlLayer,_U(slice)); collSlice; ++collSlice )  {
	      DD4hep::XML::Component compSlice = collSlice;
	      double slice_thickness = compSlice.thickness();
	      double slice_minR = lcalInnerR;
	      double slice_maxR = lcalOuterR;
	      // sensor dims differ
	      if ( compSlice.isSensitive() ){
		slice_minR = sensInnerR;
		slice_maxR = sensOuterR;
	      }

	      std::string sliceName = layer_name + DD4hep::XML::_toString(sliceID,"slice%d");
	      DD4hep::Geometry::Material   slice_material  = lcdd.material(compSlice.materialStr());
                
	      DD4hep::Geometry::Tube sliceBase(slice_minR,slice_maxR,slice_thickness/2);                
	      DD4hep::Geometry::Volume slice_vol (sliceName,sliceBase,slice_material);
                
                nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
                nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
                thickness_sum += slice_thickness/2;
              
                if ( compSlice.isSensitive() )  {
		  if ( gapDy > 0. ) {      // put gaps in sensor slice
		    DD4hep::Geometry::Box gapEnv (gapDx, gapDy , slice_thickness/2. );
		    DD4hep::Geometry::Volume gapVol ("phiGap", gapEnv, slice_material );
		    gapVol.setVisAttributes(lcdd,"PhiGapVis");
		    double rotPhi = lcalPhi0 - cellPhiSize/2.;
		    for( int gap=0; gap < nGaps; gap++ ) {
		      DD4hep::Geometry::RotationZ gapRot ( rotPhi );
		      DD4hep::Geometry::Position gapPos ( gapPosX*cos( rotPhi ), gapPosX*sin( rotPhi ), 0.0 );
		      slice_vol.placeVolume( gapVol,DD4hep::Geometry::Transform3D( gapRot, gapPos ));
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
                
                slice_vol.setAttributes(lcdd,compSlice.regionStr(),compSlice.limitsStr(),compSlice.visStr());

		DD4hep::Geometry::Position slicePos(0,0,inThisLayerPosition+slice_thickness*0.5);
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
#pragma message("FIXME: Temporary layerStaggerPhi is put in place of obsolete 'thickness'")
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
            
            caloData->layers.push_back( caloLayer ) ;
            layer_vol.setVisAttributes(lcdd,xmlLayer.visStr());
             
            DD4hep::Geometry::Position  layerPos(0,0,referencePosition+0.5*layerThickness);
	    DD4hep::Geometry::RotationZ layerRot( staggerPhi );
	    DD4hep::Geometry::PlacedVolume pv;
	    // every other layer has gaps staggered 
	    if( thisLayerId%2 == 1 && layerStagger != 0. ) {
	      pv = envelopeVol.placeVolume(layer_vol, DD4hep::Geometry::Transform3D( layerRot, layerPos ));
	    }else{
	      pv = envelopeVol.placeVolume(layer_vol,layerPos);
	    }
            pv.addPhysVolID("layer",thisLayerId);
            
            referencePosition += layerThickness;
            ++thisLayerId;
            
        }//for all layers
        
    }// for all layer collections
    std::cout<<"   LumiCal length : "<< mtotalDepthZ/dd4hep::mm 
	     << " ( rad. length X0: "<< mtotalRadLen << "   )"<<std::endl;
    std::cout << "-----------------------------------------------------------------"  << std::endl<<std::endl;;
    
    const DD4hep::Geometry::Position bcForwardPos (std::tan(0.5*fullCrossingAngle)*lcalCentreZ,0.0, lcalCentreZ);
    const DD4hep::Geometry::Position bcBackwardPos(std::tan(0.5*fullCrossingAngle)*lcalCentreZ,0.0,-lcalCentreZ);
    const DD4hep::Geometry::RotationY bcForwardRot ( fullCrossingAngle*0.5  );
    const DD4hep::Geometry::RotationY bcBackwardRot( M_PI-fullCrossingAngle*0.5 );
    
    DD4hep::Geometry::PlacedVolume pv =
    envelope.placeVolume(envelopeVol, DD4hep::Geometry::Transform3D( bcForwardRot, bcForwardPos ) );
    pv.addPhysVolID("barrel", 1);
    lumiCalDE_1.setPlacement(pv);

    DD4hep::Geometry::PlacedVolume pv2 =
    envelope.placeVolume(envelopeVol, DD4hep::Geometry::Transform3D( bcBackwardRot, bcBackwardPos ) );
    pv2.addPhysVolID("barrel", 2);
    lumiCalDE_2.setPlacement(pv2);
    
    sdet.addExtension< DD4hep::DDRec::LayeredCalorimeterData >( caloData ) ;
    
    return sdet;
}
                                               
DECLARE_DETELEMENT(LumiCal_o1_v03,create_detector)
