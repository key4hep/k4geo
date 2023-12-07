// DD4hep
#include "DD4hep/DetFactoryHelper.h"


using dd4hep::Volume;
using dd4hep::DetElement;
using dd4hep::xml::Dimension;
using dd4hep::PlacedVolume;


// todo: remove gaudi logging and properly capture output
#define endmsg std::endl
#define lLog std::cout
namespace MSG {
const std::string ERROR = " Error: ";
const std::string DEBUG = " Debug: ";
const std::string INFO  = " Info: ";
}


namespace det {

static dd4hep::Ref_t createHCal(dd4hep::Detector& lcdd, xml_det_t xmlDet, dd4hep::SensitiveDetector sensDet) {
  
  /////////////////// config parsing ///////////////////////////////////

  // Make volume that envelopes the whole barrel; set material to air
  Dimension xDimensions(xmlDet.dimensions());


  // sensitive detector type read from xml (for example "SimpleCalorimeterSD")
  Dimension xSensitive = xmlDet.child(_U(sensitive));
  sensDet.setType(xSensitive.typeStr());


  xml_comp_t xEndPlate = xmlDet.child(_Unicode(end_plate));
  double dZEndPlate = xEndPlate.thickness();
  xml_comp_t xFacePlate = xmlDet.child(_Unicode(face_plate));
  xml_comp_t xSpace = xmlDet.child(_Unicode(plate_space));  // to avoid overlaps
  double space = xSpace.thickness();
  xml_comp_t xSteelSupport = xmlDet.child(_Unicode(steel_support));
  double dSteelSupport = xSteelSupport.thickness();

  lLog << MSG::DEBUG << "steel support thickness: " << dSteelSupport << " [cm]" << endmsg;
  lLog << MSG::DEBUG << "steel support material:  " << xSteelSupport.materialStr() << endmsg;

  double sensitiveBarrelRmin = xDimensions.rmin() + xFacePlate.thickness() + space;

  // Hard-coded assumption that we have two different sequences for the modules
  std::vector<xml_comp_t> sequences = {xmlDet.child(_Unicode(sequence_a)), xmlDet.child(_Unicode(sequence_b))};
  // NOTE: This assumes that both have the same dimensions!
  Dimension sequenceDimensions(sequences[1].dimensions());
  double dzSequence = sequenceDimensions.dz();
  lLog << MSG::DEBUG << "sequence thickness " << dzSequence << endmsg;

  // calculate the number of modules fitting in Z
  unsigned int numSequencesZ = static_cast<unsigned>((2 * xDimensions.dz() - 2 * dZEndPlate - 2 * space) / dzSequence);


  // get all 'layer' children of the 'layers' tag
  std::vector<xml_comp_t> Layers;
  for (xml_coll_t xCompColl(xmlDet.child(_Unicode(layers)), _Unicode(layer)); xCompColl;
	        ++xCompColl) {
        Layers.push_back(xCompColl);
  }
  unsigned int numSequencesR = 0;
  double moduleDepth = 0.;
  std::vector<double> layerDepths = std::vector<double>();
  for (std::vector<xml_comp_t>::iterator it = Layers.begin(); it != Layers.end(); ++it) {
    xml_comp_t layer = *it;
    Dimension layerDimension(layer.dimensions());
    numSequencesR += layerDimension.nModules();
    for (int nLayer = 0; nLayer < layerDimension.nModules(); nLayer++) {
      moduleDepth += layerDimension.dr();
      layerDepths.push_back(layerDimension.dr());
    }
  }
  lLog << MSG::DEBUG << "retrieved number of layers:  " << numSequencesR
       << " , which end up to a full module depth in rho of " << moduleDepth << endmsg;
  lLog << MSG::DEBUG << "retrieved number of layers:  " << layerDepths.size() << endmsg;

  lLog << MSG::INFO << "constructing: " << numSequencesZ << " rings in Z, " << numSequencesR
       << " layers in Rho, " << numSequencesR * numSequencesZ << " tiles" << endmsg;

  // Calculate correction along z based on the module size (can only have natural number of modules)
  double dzDetector = (numSequencesZ * dzSequence) / 2 + dZEndPlate + space;
  lLog << MSG::DEBUG << "dzDetector:  " <<  dzDetector << endmsg;
  lLog << MSG::INFO << "correction of dz (negative = size reduced):" << dzDetector - xDimensions.dz() << endmsg;
  
  double rminSupport = sensitiveBarrelRmin + moduleDepth;
  double rmaxSupport = sensitiveBarrelRmin + moduleDepth + dSteelSupport;


  ////////////////////// detector building //////////////////////


  std::vector<dd4hep::PlacedVolume> layers;
  layers.reserve(layerDepths.size());
  std::vector<std::vector<dd4hep::PlacedVolume> > seqInLayers;
  seqInLayers.reserve(layerDepths.size());
  std::vector<dd4hep::PlacedVolume> tilesPerLayer;
  tilesPerLayer.reserve(layerDepths.size());


  // top level det element representing whole hcal barrel
  DetElement hCal(xmlDet.nameStr(), xmlDet.id());

  /// envelope shape
  dd4hep::Tube envelopeShape(xDimensions.rmin(), xDimensions.rmax(), xDimensions.dz());
  Volume envelopeVolume("HCalEnvelopeVolume", envelopeShape, lcdd.air());
  envelopeVolume.setVisAttributes(lcdd, xDimensions.visStr());

  // Add structural support made of steel inside of HCal
  dd4hep::Tube facePlateShape(xDimensions.rmin(), sensitiveBarrelRmin, (dzDetector - dZEndPlate - space));
  Volume facePlateVol("HCalFacePlateVol", facePlateShape, lcdd.material(xFacePlate.materialStr()));
  facePlateVol.setVisAttributes(lcdd, xFacePlate.visStr());
  PlacedVolume placedFacePlate = envelopeVolume.placeVolume(facePlateVol);
  DetElement facePlate_det(hCal, "HCalFacePlate", 0);
  facePlate_det.setPlacement(placedFacePlate);

  // Add structural support made of steel at both ends of HCal
  dd4hep::Tube endPlateShape(xDimensions.rmin(), (xDimensions.rmax() - dSteelSupport), dZEndPlate / 2);
  Volume endPlateVol("HCalEndPlateVol", endPlateShape, lcdd.material(xEndPlate.materialStr()));
  endPlateVol.setVisAttributes(lcdd, xEndPlate.visStr());

  DetElement endPlatePos(hCal, "HCalEndPlatePos", 0);
  dd4hep::Position posOffset(0, 0, dzDetector - (dZEndPlate / 2));
  PlacedVolume placedEndPlatePos = envelopeVolume.placeVolume(endPlateVol, posOffset);
  endPlatePos.setPlacement(placedEndPlatePos);

  DetElement endPlateNeg(hCal, "HCalEndPlateNeg", 1);
  dd4hep::Position negOffset(0, 0, -dzDetector + (dZEndPlate / 2));
  PlacedVolume placedEndPlateNeg = envelopeVolume.placeVolume(endPlateVol, negOffset);
  endPlateNeg.setPlacement(placedEndPlateNeg);

  dd4hep::Tube supportShape(rminSupport, rmaxSupport, (dzDetector - dZEndPlate - space));
  Volume steelSupportVolume("HCalSteelSupportVol", supportShape,
			    lcdd.material(xSteelSupport.materialStr()));
  steelSupportVolume.setVisAttributes(lcdd.invisible());
  PlacedVolume placedSupport = envelopeVolume.placeVolume(steelSupportVolume);
  DetElement support(hCal, "HCalSteelSupport", 0);
  support.setPlacement(placedSupport);

  //  double sensitiveBarrelDz = dzDetector - dZEndPlate;

  // loop over R ("layers")
  double layerR = 0.;
  for (unsigned int idxLayer = 0; idxLayer < layerDepths.size(); ++idxLayer) {
    std::string layerName = "HCalLayer" + std::to_string(idxLayer);

    // in Module rmin = 0  for first wedge, changed radius to the full radius starting at (0,0,0)
    double rminLayer = sensitiveBarrelRmin + layerR;
    double rmaxLayer = sensitiveBarrelRmin + layerR + layerDepths.at(idxLayer);
    layerR += layerDepths.at(idxLayer);

    //alternate: even layers consist of tile sequence b, odd layer of tile sequence a
    unsigned int sequenceIdx = idxLayer % 2;
    
    dd4hep::Tube tileSequenceShape(rminLayer, rmaxLayer, 0.5*dzSequence);
    Volume tileSequenceVolume("HCalTileSequenceVol", tileSequenceShape, lcdd.air());

    lLog << MSG::DEBUG << "layer radii:  " << rminLayer << " - " << rmaxLayer << " [cm]" << endmsg;
    

    dd4hep::Tube layerShape(rminLayer, rmaxLayer, dzDetector - dZEndPlate - space );
    Volume layerVolume("HCalLayerVol", layerShape, lcdd.air());

    
    layerVolume.setVisAttributes(lcdd.invisible());
    unsigned int idxSubMod = 0;
    unsigned int idxActMod = 0;
    

    dd4hep::PlacedVolume placedLayerVolume = envelopeVolume.placeVolume(layerVolume);
    placedLayerVolume.addPhysVolID("layer", idxLayer);
    layers.push_back(placedLayerVolume);
   
    double tileZOffset = - 0.5* dzSequence;
    // first Z loop (tiles that make up a sequence)
    for (xml_coll_t xCompColl(sequences[sequenceIdx], _Unicode(module_component)); xCompColl;
	        ++xCompColl, ++idxSubMod) {
      xml_comp_t xComp = xCompColl;
      dd4hep::Tube tileShape(rminLayer, rmaxLayer, 0.5 * xComp.thickness());
      
      Volume tileVol("HCalTileVol_"+xComp.nameStr(), tileShape,
            lcdd.material(xComp.materialStr()));
      tileVol.setVisAttributes(lcdd, xComp.visStr());
      
      dd4hep::Position tileOffset(0, 0, tileZOffset + 0.5 * xComp.thickness() );
      dd4hep::PlacedVolume placedTileVol = tileSequenceVolume.placeVolume(tileVol, tileOffset);
      
      if (xComp.isSensitive()) {
        tileVol.setSensitiveDetector(sensDet);
        tilesPerLayer.push_back(placedTileVol);
	idxActMod++;
      }
      tileZOffset += xComp.thickness();
    }

    // second z loop (place sequences in layer)
    std::vector<dd4hep::PlacedVolume> sq_vector; 

    for (uint numSeq=0; numSeq < numSequencesZ; numSeq++){
      double zOffset = - dzDetector + (2 * numSeq + 1) * (dzSequence * 0.5);
      dd4hep::Position tileSequencePosition(0, 0, zOffset);
      dd4hep::PlacedVolume placedTileSequenceVolume = layerVolume.placeVolume(tileSequenceVolume, tileSequencePosition);
      placedTileSequenceVolume.addPhysVolID("row", numSeq);
      sq_vector.push_back(placedTileSequenceVolume);
    }
    seqInLayers.push_back(sq_vector);

  }

  // Place det elements wihtin each other to recover volume positions later via cellID  
  for (uint iLayer = 0; iLayer < numSequencesR; iLayer++) {
    DetElement layerDet(hCal, dd4hep::xml::_toString(iLayer, "layer%d"), iLayer);
    layerDet.setPlacement(layers[iLayer]);
    
    for (uint iSeq = 0; iSeq < seqInLayers[iLayer].size(); iSeq++){
      DetElement seqDet(layerDet, dd4hep::xml::_toString(iSeq, "seq%d"), iSeq);
      seqDet.setPlacement(seqInLayers[iLayer][iSeq]);

      DetElement tileDet(seqDet, dd4hep::xml::_toString(iSeq, "tile%d"), iSeq);
      tileDet.setPlacement(tilesPerLayer[iLayer]);
    }
  }
  
  // Place envelope (or barrel) volume
  Volume motherVol = lcdd.pickMotherVolume(hCal);
  motherVol.setVisAttributes(lcdd.invisible());
  PlacedVolume placedHCal = motherVol.placeVolume(envelopeVolume);
  placedHCal.addPhysVolID("system", hCal.id());
  hCal.setPlacement(placedHCal);
  return hCal;
}
}  // namespace hcal

DECLARE_DETELEMENT(HCalTileBarrel_o1_v01, det::createHCal)
