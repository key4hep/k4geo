// DD4hep
#include "DD4hep/DetFactoryHelper.h"

// todo: remove gaudi logging and properly capture output
#define endmsg std::endl
#define lLog std::cout
namespace MSG {
const std::string DEBUG = " Debug: ";
const std::string INFO  = " Info: ";
}

using dd4hep::Volume;
using dd4hep::DetElement;
using dd4hep::xml::Dimension;
using dd4hep::PlacedVolume;

namespace det {
void buildEC(dd4hep::Detector& aLcdd, dd4hep::SensitiveDetector& aSensDet, dd4hep::Volume& aEnvelope,
             dd4hep::DetElement& aHCal, xml_det_t aXmlElement, int sign) {

  dd4hep::SensitiveDetector sensDet = aSensDet;
  Dimension sensDetType = aXmlElement.child(_Unicode(sensitive));
  sensDet.setType(sensDetType.typeStr());

  Dimension dimensions(aXmlElement.child(_Unicode(dimensions)));
  xml_comp_t xEndPlate = aXmlElement.child(_Unicode(end_plate));
  double dZEndPlate = xEndPlate.thickness() / 2.;
  xml_comp_t xFacePlate = aXmlElement.child(_Unicode(face_plate));
  double dRhoFacePlate = xFacePlate.thickness() / 2.;
  xml_comp_t xSpace = aXmlElement.child(_Unicode(plate_space));  // to avoid overlaps
  double space = xSpace.thickness();
  xml_comp_t xSteelSupport = aXmlElement.child(_Unicode(steel_support));
  double dSteelSupport = xSteelSupport.thickness();
  lLog << MSG::DEBUG << "steel support thickness " << dSteelSupport << endmsg;
  lLog << MSG::DEBUG << "steel support material  " << xSteelSupport.materialStr() << endmsg;

  double sensitiveBarrel1Rmin = dimensions.rmin1() + 2 * dRhoFacePlate + space;
  double sensitiveBarrel2Rmin = dimensions.rmin2() + 2 * dRhoFacePlate + space;
  double sensitiveBarrel3Rmin = dimensions.rmin() + 2 * dRhoFacePlate + space;

  // Offset in z is given as distance from 0 to the middle of the Calorimeter volume
  double extBarrelOffset1 = dimensions.offset();
  double extBarrelOffset2 = dimensions.z_offset();
  double extBarrelOffset3 = dimensions.v_offset();

  // Hard-coded assumption that we have two different sequences for the modules
  std::vector<xml_comp_t> sequences = {aXmlElement.child(_Unicode(sequence_a)), aXmlElement.child(_Unicode(sequence_b))};
  // NOTE: This assumes that both have the same dimensions!
  Dimension sequenceDimensions(sequences[1].dimensions());
  double dzSequence = sequenceDimensions.dz();
  lLog << MSG::DEBUG << "sequence thickness " << dzSequence << endmsg;

  // calculate the number of modules fitting in  Z
  unsigned int numSequencesZ1 = static_cast<unsigned>((2 * dimensions.width() - 2 * dZEndPlate - space) / dzSequence);
  unsigned int numSequencesZ2 = static_cast<unsigned>((2 * dimensions.dz() - 2 * dZEndPlate - space) / dzSequence);
  unsigned int numSequencesZ3 = static_cast<unsigned>((2 * dimensions.z_length() - 2 * dZEndPlate - space) / dzSequence);

  unsigned int numSequencesR1 = 0;
  unsigned int numSequencesR2 = 0;
  unsigned int numSequencesR3 = 0;
  double moduleDepth1 = 0.;
  double moduleDepth2 = 0.;
  double moduleDepth3 = 0.;
  std::vector<double> layerDepths1 = std::vector<double>();
  std::vector<double> layerDepths2 = std::vector<double>();
  std::vector<double> layerDepths3 = std::vector<double>();

  // get all 'layer' children of the 'layers' tag
  std::vector<xml_comp_t> Layers;
  for (xml_coll_t xCompColl(aXmlElement.child(_Unicode(layers)), _Unicode(layer)); xCompColl;
       ++xCompColl) {
    Layers.push_back(xCompColl);
  }
  
  for (std::vector<xml_comp_t>::iterator it = Layers.begin(); it != Layers.end(); ++it) {
    xml_comp_t layer = *it;
    Dimension layerDimension(layer.dimensions());
    numSequencesR1 += layerDimension.nmodules();
    numSequencesR2 += layerDimension.nModules();
    numSequencesR3 += layerDimension.nPads();
    for (int nLayer = 0; nLayer < layerDimension.nmodules(); nLayer++) {
      moduleDepth1 += layerDimension.dr();
      layerDepths1.push_back(layerDimension.dr());
    }
    for (int nLayer = 0; nLayer < layerDimension.nModules(); nLayer++) {
      moduleDepth2 += layerDimension.dr();
      layerDepths2.push_back(layerDimension.dr());
    }
    for (int nLayer = 0; nLayer < layerDimension.nPads(); nLayer++) {
      moduleDepth3 += layerDimension.dr();
      layerDepths3.push_back(layerDimension.dr());
    }
  }

  lLog << MSG::DEBUG << "retrieved number of layers in first Endcap part:  " << numSequencesR1
       << " , which end up to a full module depth in rho of " << moduleDepth1 << endmsg;
  lLog << MSG::DEBUG << "retrieved number of layers in first Endcap part:  " << layerDepths1.size() << endmsg;
  lLog << MSG::DEBUG << "retrieved number of layers in second Endcap part:  " << numSequencesR2
       << " , which end up to a full module depth in rho of " << moduleDepth2 << endmsg;
  lLog << MSG::DEBUG << "retrieved number of layers in second Endcap part:  " << layerDepths2.size() << endmsg;
  lLog << MSG::DEBUG << "retrieved number of layers in third Endcap part:  " << numSequencesR3
       << " , which end up to a full module depth in rho of " << moduleDepth3 << endmsg;
  lLog << MSG::DEBUG << "retrieved number of layers in third Endcap part:  " << layerDepths3.size() << endmsg;

  lLog << MSG::INFO << "constructing first part EC: with offset " << extBarrelOffset1 << ": "<< numSequencesZ1
       << " rings in Z, " << numSequencesR1 << " layers in Rho, " << numSequencesR1 * numSequencesZ1
       << " tiles" << endmsg;

  lLog << MSG::INFO << "constructing second part EC: with offset " << extBarrelOffset2 << ": " << numSequencesZ2
       << " rings in Z, " << numSequencesR2 << " layers in Rho, "
       << layerDepths2.size() * numSequencesZ2 << " tiles" << endmsg;

   lLog << MSG::INFO << "constructing third part EC: with offset " << extBarrelOffset3 << ": " << numSequencesZ3
       << " rings in Z, " << numSequencesR3 << " layers in Rho, "
       << layerDepths3.size() * numSequencesZ3 << " tiles" << endmsg;

  lLog << MSG::INFO << "number of channels: "
       << (numSequencesR1 * numSequencesZ1) + (numSequencesR2 * numSequencesZ2) + (numSequencesR3 * numSequencesZ3)
       << endmsg;

  // Calculate correction along z based on the module size (can only have natural number of modules)
  double dzDetector1 = (numSequencesZ1 * dzSequence) / 2 + 2 * dZEndPlate + space;
  lLog << MSG::INFO
       << "correction of dz (negative = size reduced) first part EC :" << dzDetector1*2 - dimensions.width()*2
       << endmsg;
  double dzDetector2 = (numSequencesZ2 * dzSequence) / 2;
  lLog << MSG::INFO << "dz second part EC:" << dzDetector2 * 2
       << endmsg;
  lLog << MSG::INFO << "width second part EC:" << dimensions.dz() * 2
       << endmsg;
  lLog << MSG::INFO << "correction of dz (negative = size reduced) second part EB:" << dzDetector2*2 - dimensions.dz()*2
       << endmsg;

  double dzDetector3 = (numSequencesZ3 * dzSequence) / 2 + 2 * dZEndPlate + space;
  lLog << MSG::INFO << "dz third part EC:" << dzDetector2 * 2
       << endmsg;

  // Add structural support made of steel inside of HCal
  DetElement facePlate1(aHCal, "FacePlate_" + std::to_string(1 * sign), 0);
  dd4hep::Tube facePlateShape1(dimensions.rmin1(), (sensitiveBarrel1Rmin - space),
                               (dzDetector1 - 2 * dZEndPlate - space));
  Volume facePlateVol1("facePlateVol1", facePlateShape1, aLcdd.material(xFacePlate.materialStr()));
  facePlateVol1.setVisAttributes(aLcdd, xFacePlate.visStr());
  dd4hep::Position offsetFace1(0, 0, sign * extBarrelOffset1);

  // Faceplate for 2nd part of extended Barrel
  DetElement facePlate2(aHCal, "FacePlate_" + std::to_string(2 * sign), 0);
  dd4hep::Tube facePlateShape2(dimensions.rmin2(), (sensitiveBarrel2Rmin - space),
                               dzDetector2);
  Volume facePlateVol2("facePlateVol2", facePlateShape2, aLcdd.material(xFacePlate.materialStr()));
  facePlateVol2.setVisAttributes(aLcdd, xFacePlate.visStr());
  dd4hep::Position offsetFace2(0, 0, sign * extBarrelOffset2);

  // Faceplate for 3rd part of extended Barrel
  DetElement facePlate3(aHCal, "FacePlate_" + std::to_string(3 * sign), 0);
  dd4hep::Tube facePlateShape3(dimensions.rmin(), (sensitiveBarrel3Rmin - space),
                               (dzDetector3 - 2 * dZEndPlate - space));
  Volume facePlateVol3("facePlateVol3", facePlateShape3, aLcdd.material(xFacePlate.materialStr()));
  facePlateVol3.setVisAttributes(aLcdd, xFacePlate.visStr());
  dd4hep::Position offsetFace3(0, 0, sign * extBarrelOffset3);

  PlacedVolume placedFacePlate1 = aEnvelope.placeVolume(facePlateVol1, offsetFace1);
  facePlate1.setPlacement(placedFacePlate1);
  PlacedVolume placedFacePlate2 = aEnvelope.placeVolume(facePlateVol2, offsetFace2);
  facePlate2.setPlacement(placedFacePlate2);
  PlacedVolume placedFacePlate3 = aEnvelope.placeVolume(facePlateVol3, offsetFace3);
  facePlate3.setPlacement(placedFacePlate3);

  // Add structural support made of steel at both ends of extHCal
  dd4hep::Tube endPlateShape1(dimensions.rmin1(), (dimensions.rmax1() - dSteelSupport), dZEndPlate);
  Volume endPlateVol1("endPlateVol1", endPlateShape1, aLcdd.material(xEndPlate.materialStr()));
  endPlateVol1.setVisAttributes(aLcdd, xEndPlate.visStr());
   dd4hep::Tube endPlateShape3(dimensions.rmin(), (dimensions.rmax() - dSteelSupport), dZEndPlate);
  Volume endPlateVol3("endPlateVol3", endPlateShape3, aLcdd.material(xEndPlate.materialStr()));
  endPlateVol3.setVisAttributes(aLcdd, xEndPlate.visStr());

  // Endplates placed for the extended Barrels in front and in the back to the central Barrel
  DetElement endPlatePos(aHCal, "endPlate_" + std::to_string(1 * sign), 0);
  dd4hep::Position posOffset(0, 0, sign * (extBarrelOffset3 + dzDetector3 - dZEndPlate));
  PlacedVolume placedEndPlatePos = aEnvelope.placeVolume(endPlateVol3, posOffset);
  endPlatePos.setPlacement(placedEndPlatePos);

  DetElement endPlateNeg(aHCal, "endPlate_" + std::to_string(2 * sign), 0);
  dd4hep::Position negOffset(0, 0, sign * (extBarrelOffset1 - dzDetector1 + dZEndPlate));
  PlacedVolume placedEndPlateNeg = aEnvelope.placeVolume(endPlateVol1, negOffset);
  endPlateNeg.setPlacement(placedEndPlateNeg);

  std::vector<dd4hep::PlacedVolume> layers;
  layers.reserve(layerDepths1.size()+layerDepths2.size()+layerDepths3.size());
  std::vector<std::vector<dd4hep::PlacedVolume> > seqInLayers;
  seqInLayers.reserve(layerDepths1.size()+layerDepths2.size()+layerDepths3.size());
  std::vector<dd4hep::PlacedVolume> tilesPerLayer;
  tilesPerLayer.reserve(layerDepths1.size()+layerDepths2.size()+layerDepths3.size());

  // loop over R ("layers")
  double layerR = 0.;
  for (unsigned int idxLayer = 0; idxLayer < layerDepths1.size(); ++idxLayer) {
    // in Module rmin = 0  for first wedge, changed radius to the full radius starting at (0,0,0)
    double rminLayer = sensitiveBarrel1Rmin + layerR;
    double rmaxLayer = sensitiveBarrel1Rmin + layerR + layerDepths1.at(idxLayer);
    layerR += layerDepths1.at(idxLayer);

    //alternate: even layers consist of tile sequence b, odd layer of tile sequence a
    unsigned int sequenceIdx = (idxLayer+1) % 2;
    
    dd4hep::Tube tileSequenceShape(rminLayer, rmaxLayer, 0.5*dzSequence);
    Volume tileSequenceVolume("HCalECTileSequenceVol1", tileSequenceShape, aLcdd.air());

    lLog << MSG::DEBUG << "layer radii:  " << rminLayer << " - " << rmaxLayer << " [cm]" << endmsg;
    

    dd4hep::Tube layerShape(rminLayer, rmaxLayer, dzDetector1 );
    Volume layerVolume("HCalECLayerVol1", layerShape, aLcdd.air());

  
    layerVolume.setVisAttributes(aLcdd.invisible());
    unsigned int idxSubMod = 0;

    dd4hep::Position moduleOffset1 (0,0,sign * extBarrelOffset1);

    dd4hep::PlacedVolume placedLayerVolume = aEnvelope.placeVolume(layerVolume, moduleOffset1);
    unsigned int type1 = 0;
    if (sign<0) {
      type1 = 3;
    }
    placedLayerVolume.addPhysVolID("type", type1);  // First module type=0,3 in front of second +/-
    placedLayerVolume.addPhysVolID("layer", idxLayer);
    layers.push_back(placedLayerVolume);
   
    double tileZOffset = - 0.5* dzSequence;

    // first Z loop (tiles that make up a sequence)
    for (xml_coll_t xCompColl(sequences[sequenceIdx], _Unicode(module_component)); xCompColl;
	 ++xCompColl, ++idxSubMod) {
      xml_comp_t xComp = xCompColl;
      dd4hep::Tube tileShape(rminLayer, rmaxLayer, 0.5 * xComp.thickness());
      
      Volume tileVol("HCalECTileVol_"+ xComp.materialStr()
		     , tileShape,
		     aLcdd.material(xComp.materialStr()));
      tileVol.setVisAttributes(aLcdd, xComp.visStr());
      
      dd4hep::Position tileOffset(0, 0, tileZOffset + 0.5 * xComp.thickness() );
      dd4hep::PlacedVolume placedTileVol = tileSequenceVolume.placeVolume(tileVol, tileOffset);
      
      if (xComp.isSensitive()) {
        tileVol.setSensitiveDetector(sensDet);
        tilesPerLayer.push_back(placedTileVol);
      }
      tileZOffset += xComp.thickness();
    }

    // second z loop (place sequences in layer)
    std::vector<dd4hep::PlacedVolume> seqs; 
    double zOffset = - dzDetector1 + 0.5 * dzSequence; //2*dZEndPlate + space + 0.5 * dzSequence;
    
    for (uint numSeq=0; numSeq < numSequencesZ1; numSeq++){
      dd4hep::Position tileSequencePosition(0, 0, zOffset);
      dd4hep::PlacedVolume placedTileSequenceVolume = layerVolume.placeVolume(tileSequenceVolume, tileSequencePosition);
      placedTileSequenceVolume.addPhysVolID("row", numSeq);
      seqs.push_back(placedTileSequenceVolume);
      zOffset += dzSequence;
    }
    seqInLayers.push_back(seqs);
  }


  layerR = 0.;
  // Placement of subWedges in Wedge, 2nd part
  for (unsigned int idxLayer = 0; idxLayer < layerDepths2.size(); ++idxLayer) {
    // in Module rmin = 0  for first wedge, changed radius to the full radius starting at (0,0,0)                                                                      
    double rminLayer = sensitiveBarrel2Rmin + layerR;
    double rmaxLayer = sensitiveBarrel2Rmin + layerR + layerDepths2.at(idxLayer);
    layerR += layerDepths2.at(idxLayer);

    //alternate: even layers consist of tile sequence b, odd layer of tile sequence a                                                                                  
    unsigned int sequenceIdx = (idxLayer+1) % 2;

    dd4hep::Tube tileSequenceShape(rminLayer, rmaxLayer, 0.5*dzSequence);
    Volume tileSequenceVolume("HCalECTileSequenceVol2", tileSequenceShape, aLcdd.air());

    lLog << MSG::DEBUG << "layer radii:  " << rminLayer << " - " << rmaxLayer << " [cm]" << endmsg;


    dd4hep::Tube layerShape(rminLayer, rmaxLayer, dzDetector2);
    Volume layerVolume("HCalECLayerVol2", layerShape, aLcdd.air());

    layerVolume.setVisAttributes(aLcdd.invisible());
    unsigned int idxSubMod = 0;

    double tileZOffset = - 0.5* dzSequence;

    // first Z loop (tiles that make up a sequence)                                                                                                                    
    for (xml_coll_t xCompColl(sequences[sequenceIdx], _Unicode(module_component)); xCompColl;
         ++xCompColl, ++idxSubMod) {
      xml_comp_t xComp = xCompColl;
      dd4hep::Tube tileShape(rminLayer, rmaxLayer, 0.5 * xComp.thickness());

      Volume tileVol("HCalECTileVol_"
		     , tileShape,
                     aLcdd.material(xComp.materialStr()));
      tileVol.setVisAttributes(aLcdd, xComp.visStr());

      dd4hep::Position tileOffset(0, 0, tileZOffset + 0.5 * xComp.thickness() );
      dd4hep::PlacedVolume placedTileVol = tileSequenceVolume.placeVolume(tileVol, tileOffset);

      if (xComp.isSensitive()) {
        tileVol.setSensitiveDetector(sensDet);
        tilesPerLayer.push_back(placedTileVol);
      }
      tileZOffset += xComp.thickness();
    }

    // second z loop (place sequences in layer)                                                                                                                        
    std::vector<dd4hep::PlacedVolume> seqs;
    double zOffset = - dzDetector2 + 0.5 * dzSequence; //(dzSequence * 0.5);

    for (uint numSeq=0; numSeq < numSequencesZ2; numSeq++){
      dd4hep::Position tileSequencePosition(0, 0, zOffset);
      dd4hep::PlacedVolume placedTileSequenceVolume = layerVolume.placeVolume(tileSequenceVolume, tileSequencePosition);
      placedTileSequenceVolume.addPhysVolID("row", numSeq);
      seqs.push_back(placedTileSequenceVolume);
      zOffset += dzSequence;
    }
    seqInLayers.push_back(seqs);

    dd4hep::Position moduleOffset2 (0, 0, sign * extBarrelOffset2);
    dd4hep::PlacedVolume placedLayerVolume = aEnvelope.placeVolume(layerVolume, moduleOffset2);
    unsigned int type2 = 1;
    if (sign<0) {
      type2 = 4;
    }
    placedLayerVolume.addPhysVolID("type", type2);  // Second module type=1,4 behind the first +/-
    placedLayerVolume.addPhysVolID("layer", layerDepths1.size() + idxLayer);
    layers.push_back(placedLayerVolume);
  }

  layerR = 0.;
  // Placement of subWedges in Wedge, 3th part
  for (unsigned int idxLayer = 0; idxLayer < layerDepths3.size(); ++idxLayer) {
    // in Module rmin = 0  for first wedge, changed radius to the full radius starting at (0,0,0)                                                                      
    double rminLayer = sensitiveBarrel3Rmin + layerR;
    double rmaxLayer = sensitiveBarrel3Rmin + layerR + layerDepths3.at(idxLayer);
    layerR += layerDepths3.at(idxLayer);
    
    //alternate: even layers consist of tile sequence b, odd layer of tile sequence a                                                                                  
    unsigned int sequenceIdx = (idxLayer+1) % 2;

    dd4hep::Tube tileSequenceShape(rminLayer, rmaxLayer, 0.5*dzSequence);
    Volume tileSequenceVolume("HCalECTileSequenceVol3", tileSequenceShape, aLcdd.air());

    lLog << MSG::DEBUG << "layer radii:  " << rminLayer << " - " << rmaxLayer << " [cm]" << endmsg;

    dd4hep::Tube layerShape(rminLayer, rmaxLayer, dzDetector3);
    Volume layerVolume("HCalECLayerVol3", layerShape, aLcdd.air());

    layerVolume.setVisAttributes(aLcdd.invisible());
    unsigned int idxSubMod = 0;

    double tileZOffset = - 0.5* dzSequence;

    // first Z loop (tiles that make up a sequence)                                                                                                                    
    for (xml_coll_t xCompColl(sequences[sequenceIdx], _Unicode(module_component)); xCompColl;
	 ++xCompColl, ++idxSubMod) {
      xml_comp_t xComp = xCompColl;
      dd4hep::Tube tileShape(rminLayer, rmaxLayer, 0.5 * xComp.thickness());
      
      Volume tileVol("HCalECTileVol_"
		     , tileShape,
		     aLcdd.material(xComp.materialStr()));
      tileVol.setVisAttributes(aLcdd, xComp.visStr());

      dd4hep::Position tileOffset(0, 0, tileZOffset + 0.5 * xComp.thickness() );
      dd4hep::PlacedVolume placedTileVol = tileSequenceVolume.placeVolume(tileVol, tileOffset);

      if (xComp.isSensitive()) {
	tileVol.setSensitiveDetector(sensDet);
	tilesPerLayer.push_back(placedTileVol);
      }
      tileZOffset += xComp.thickness();
    }
       
    // second z loop (place sequences in layer)                                                                                                                        
    std::vector<dd4hep::PlacedVolume> seqs;
    double zOffset = - dzDetector3 + 0.5 * dzSequence; //2*dZEndPlate + space + (dzSequence * 0.5);

    for (uint numSeq=0; numSeq < numSequencesZ3; numSeq++){
      dd4hep::Position tileSequencePosition(0, 0, zOffset);
      dd4hep::PlacedVolume placedTileSequenceVolume = layerVolume.placeVolume(tileSequenceVolume, tileSequencePosition);
      placedTileSequenceVolume.addPhysVolID("row", numSeq);
      seqs.push_back(placedTileSequenceVolume);
      zOffset += dzSequence;
    }
    seqInLayers.push_back(seqs);

    dd4hep::Position moduleOffset3 (0, 0, sign * extBarrelOffset3);
    dd4hep::PlacedVolume placedLayerVolume = aEnvelope.placeVolume(layerVolume, moduleOffset3);
    unsigned int type3 = 2;
    if (sign<0) {
      type3 = 5;
    }
    placedLayerVolume.addPhysVolID("type", type3);  // Second module type=2,5 behind the first +/-
    placedLayerVolume.addPhysVolID("layer", layerDepths1.size() + layerDepths2.size() + idxLayer);
    layers.push_back(placedLayerVolume);
  }

  // Placement of DetElements
  lLog << MSG::DEBUG << "Layers in r :    " << layers.size() << std::endl;
  lLog << MSG::DEBUG << "Tiles in layers :" << tilesPerLayer.size() << std::endl;
  
  for (uint iLayer = 0; iLayer < (layerDepths1.size()+layerDepths2.size()+layerDepths3.size()); iLayer++) {
    DetElement layerDet(aHCal, dd4hep::xml::_toString(sign*(iLayer+1), "layer%d"), sign*(iLayer+1));
    layerDet.setPlacement(layers[iLayer]);
    
    for (uint iSeq = 0; iSeq < seqInLayers[iLayer].size(); iSeq++){
      DetElement seqDet(layerDet, dd4hep::xml::_toString(iSeq, "seq%d"), sign*(iSeq+1));
      seqDet.setPlacement(seqInLayers[iLayer][iSeq]);

      DetElement tileDet(seqDet, dd4hep::xml::_toString(iSeq, "tile%d"), sign*(iSeq+1));
      tileDet.setPlacement(tilesPerLayer[iLayer]);
    }
  }  

}

static dd4hep::Ref_t createHCalEC(dd4hep::Detector& lcdd, xml_h xmlElement, dd4hep::SensitiveDetector sensDet) {


  xml_det_t xmlDet = xmlElement;
  std::string detName = xmlDet.nameStr();

  // Make DetElement
  dd4hep::DetElement hCalEC(detName, xmlDet.id());

  // Make volume that envelopes the whole barrel; set material to air
  Dimension dimensions(xmlDet.dimensions());

  dd4hep::Tube envelope(dimensions.rmin(), dimensions.rmax1(), (dimensions.v_offset() + dimensions.z_length()));
  dd4hep::Tube negative1(dimensions.rmin(), dimensions.rmax1(), (dimensions.offset() - dimensions.width()));
  dd4hep::Tube negative2(dimensions.rmin(), dimensions.rmin1(), (dimensions.z_offset() - dimensions.dz()));
  dd4hep::Tube negative3(dimensions.rmin(), dimensions.rmin2(), (dimensions.v_offset() - dimensions.z_length()));
  dd4hep::SubtractionSolid envelopeShapeTmp1(envelope, negative1);
  dd4hep::SubtractionSolid envelopeShapeTmp2(envelopeShapeTmp1, negative2);
  dd4hep::SubtractionSolid envelopeShape(envelopeShapeTmp2, negative3);

  Volume envelopeVolume(detName + "_volume", envelopeShape, lcdd.air());
  envelopeVolume.setVisAttributes(lcdd, dimensions.visStr());

  lLog << MSG::DEBUG << "Placing detector on the positive side: (cm) " << (dimensions.offset() + dimensions.dz())
       << endmsg;
  buildEC(lcdd, sensDet, envelopeVolume, hCalEC, xmlElement, 1);
  lLog << MSG::DEBUG << "Placing detector on the negative side: (cm) " << -(dimensions.offset() + dimensions.dz())
      << endmsg;
  buildEC(lcdd, sensDet, envelopeVolume, hCalEC, xmlElement, -1);

  // Place envelope volume
  Volume motherVol = lcdd.pickMotherVolume(hCalEC);

  PlacedVolume placedHCal = motherVol.placeVolume(envelopeVolume);
  placedHCal.addPhysVolID("system", xmlDet.id());
  hCalEC.setPlacement(placedHCal);

  return hCalEC;
}
}  // namespace hcal

DECLARE_DETELEMENT(CaloThreePartsEndcap_o1_v01, det::createHCalEC)
