// DD4hep
#include "DD4hep/DetFactoryHelper.h"
#include <DDRec/DetectorData.h>
#include "DD4hep/Printout.h"

using dd4hep::Volume;
using dd4hep::DetElement;
using dd4hep::xml::Dimension;
using dd4hep::PlacedVolume;

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

  dd4hep::printout(dd4hep::DEBUG, "HCalTileBarrel_o1_v02", "steel support thickness (cm): %.2f", dSteelSupport); 

  double sensitiveBarrelRmin = xDimensions.rmin() + xFacePlate.thickness() + space;

  // Hard-coded assumption that we have two different sequences for the modules
  std::vector<xml_comp_t> sequences = {xmlDet.child(_Unicode(sequence_a)), xmlDet.child(_Unicode(sequence_b))};
  // Check if both sequences are present
  if (!sequences[0] || !sequences[1]) {
    dd4hep::printout(dd4hep::ERROR, "HCalTileBarrel_o1_v02", "The two sequences 'sequence_a' and 'sequence_b' must be present in the xml file."); 
    throw std::runtime_error("Missing sequence_a or sequence_b in the xml file.");
  }
  // Check if both sequences have the same dimensions
  Dimension dimensionsA(sequences[0].dimensions());
  Dimension dimensionsB(sequences[1].dimensions());
  if (dimensionsA.dz() != dimensionsB.dz()) {
    dd4hep::printout(dd4hep::ERROR, "HCalTileBarrel_o1_v02", "The dimensions of sequence_a and sequence_b do not match."); 
    throw std::runtime_error("The dimensions of the sequence_a and sequence_b do not match.");
  }
  double dzSequence = dimensionsB.dz();
  dd4hep::printout(dd4hep::DEBUG, "HCalTileBarrel_o1_v02", "sequence thickness %.2f", dzSequence); 

  // calculate the number of sequences fitting in Z
  unsigned int numSequencesZ = static_cast<unsigned>((2 * xDimensions.dz() - 2 * dZEndPlate - 2 * space) / dzSequence);

  // get all 'layer' children of the 'layers' tag
  std::vector<xml_comp_t> Layers;
  for (xml_coll_t xCompColl(xmlDet.child(_Unicode(layers)), _Unicode(layer)); xCompColl; ++xCompColl) {
    Layers.push_back(xCompColl);
  }
  unsigned int numLayersR = 0;
  double moduleDepth = 0.;
  std::vector<double> layerDepths = std::vector<double>();
  std::vector<double> layerInnerRadii = std::vector<double>();
  for (std::vector<xml_comp_t>::iterator it = Layers.begin(); it != Layers.end(); ++it) {
    xml_comp_t layer = *it;
    Dimension layerDimension(layer.dimensions());
    numLayersR += layerDimension.nModules();
    for (int nLayer = 0; nLayer < layerDimension.nModules(); nLayer++) {
      moduleDepth += layerDimension.dr();
      layerDepths.push_back(layerDimension.dr());
    }
  }
  // Calculate correction along z based on the module size (can only have natural number of modules)
  double dzDetector = (numSequencesZ * dzSequence) / 2 + dZEndPlate + space;
  
  dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "dzDetector (cm): %.2f", dzDetector); 
  dd4hep::printout(dd4hep::DEBUG, "HCalTileBarrel_o1_v02", "correction of dz in cm (negative = size reduced): %.2f", dzDetector - xDimensions.dz()); 

  dd4hep::printout(dd4hep::DEBUG, "HCalTileBarrel_o1_v02", "retrieved number of radial layers: %d , which end up to a full module depth in rho of %.2f cm", numLayersR, moduleDepth); 
  dd4hep::printout(dd4hep::DEBUG, "HCalTileBarrel_o1_v02", "retrieved number of radial layers: %d", layerDepths.size()); 
  dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "constructing: %d sequences in Z, %d radial layers, in total %d tiles", numSequencesZ, numLayersR, numLayersR * numSequencesZ); 
  
  double rminSupport = sensitiveBarrelRmin + moduleDepth;
  double rmaxSupport = sensitiveBarrelRmin + moduleDepth + dSteelSupport;

  double sensitiveBarrelRmax = sensitiveBarrelRmin + moduleDepth;


  ////////////////////// detector building //////////////////////


  std::vector<dd4hep::PlacedVolume> layers;
  layers.reserve(layerDepths.size());
  std::vector<std::vector<dd4hep::PlacedVolume> > seqInLayers;
  seqInLayers.reserve(layerDepths.size());
  std::vector<dd4hep::PlacedVolume> tilesPerLayer;
  tilesPerLayer.reserve(layerDepths.size());


  // top level det element representing whole hcal barrel
  DetElement caloDetElem(xmlDet.nameStr(), xmlDet.id());

  /// envelope shape
  dd4hep::Tube envelopeShape(xDimensions.rmin(), xDimensions.rmax(), xDimensions.dz());

  Volume envelopeVolume("HCalEnvelopeVolume", envelopeShape, lcdd.air());
  envelopeVolume.setVisAttributes(lcdd, xDimensions.visStr());

  // Add structural support made of steel inside of HCal
  dd4hep::Tube facePlateShape(xDimensions.rmin(), sensitiveBarrelRmin, (dzDetector - dZEndPlate - space));
  Volume facePlateVol("HCalFacePlateVol", facePlateShape, lcdd.material(xFacePlate.materialStr()));
  facePlateVol.setVisAttributes(lcdd, xFacePlate.visStr());
  PlacedVolume placedFacePlate = envelopeVolume.placeVolume(facePlateVol);
  DetElement facePlate_det(caloDetElem, "HCalFacePlate", 0);
  facePlate_det.setPlacement(placedFacePlate);

  // Add structural support made of steel at both ends of HCal
  dd4hep::Tube endPlateShape(xDimensions.rmin(), (xDimensions.rmax() - dSteelSupport), dZEndPlate / 2);
  Volume endPlateVol("HCalEndPlateVol", endPlateShape, lcdd.material(xEndPlate.materialStr()));
  endPlateVol.setVisAttributes(lcdd, xEndPlate.visStr());

  DetElement endPlatePos(caloDetElem, "HCalEndPlatePos", 0);
  dd4hep::Position posOffset(0, 0, dzDetector - (dZEndPlate / 2));
  PlacedVolume placedEndPlatePos = envelopeVolume.placeVolume(endPlateVol, posOffset);
  endPlatePos.setPlacement(placedEndPlatePos);

  DetElement endPlateNeg(caloDetElem, "HCalEndPlateNeg", 1);
  dd4hep::Position negOffset(0, 0, -dzDetector + (dZEndPlate / 2));
  PlacedVolume placedEndPlateNeg = envelopeVolume.placeVolume(endPlateVol, negOffset);
  endPlateNeg.setPlacement(placedEndPlateNeg);

  dd4hep::Tube supportShape(rminSupport, rmaxSupport, (dzDetector - dZEndPlate - space));
  Volume steelSupportVolume("HCalSteelSupportVol", supportShape, lcdd.material(xSteelSupport.materialStr()));
  steelSupportVolume.setVisAttributes(lcdd.invisible());
  PlacedVolume placedSupport = envelopeVolume.placeVolume(steelSupportVolume);
  DetElement support(caloDetElem, "HCalSteelSupport", 0);
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
    layerInnerRadii.push_back(rminLayer);
    dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "layer %d (cm): %.2f - %.2f", idxLayer, rminLayer, rmaxLayer); 

    //alternate: even layers consist of tile sequence b, odd layer of tile sequence a
    unsigned int sequenceIdx = idxLayer % 2;
    
    dd4hep::Tube tileSequenceShape(rminLayer, rmaxLayer, 0.5*dzSequence);
    Volume tileSequenceVolume("HCalTileSequenceVol", tileSequenceShape, lcdd.air());    

    dd4hep::Tube layerShape(rminLayer, rmaxLayer, dzDetector - dZEndPlate - space );
    Volume layerVolume("HCalLayerVol", layerShape, lcdd.air());

    
    layerVolume.setVisAttributes(lcdd.invisible());
    unsigned int idxSubMod = 0;
    

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
  for (uint iLayer = 0; iLayer < numLayersR; iLayer++) {
    DetElement layerDet(caloDetElem, dd4hep::xml::_toString(iLayer, "layer%d"), iLayer);
    layerDet.setPlacement(layers[iLayer]);
    
    for (uint iSeq = 0; iSeq < seqInLayers[iLayer].size(); iSeq++){
      DetElement seqDet(layerDet, dd4hep::xml::_toString(iSeq, "seq%d"), iSeq);
      seqDet.setPlacement(seqInLayers[iLayer][iSeq]);

      DetElement tileDet(seqDet, dd4hep::xml::_toString(iSeq, "tile%d"), iSeq);
      tileDet.setPlacement(tilesPerLayer[iLayer]);
    }
  }
  
  // Place envelope (or barrel) volume
  Volume motherVol = lcdd.pickMotherVolume(caloDetElem);
  motherVol.setVisAttributes(lcdd.invisible());
  PlacedVolume envelopePhysVol = motherVol.placeVolume(envelopeVolume);
  envelopePhysVol.addPhysVolID("system", caloDetElem.id());
  caloDetElem.setPlacement(envelopePhysVol);

  
  // Create caloData object
  auto caloData = new dd4hep::rec::LayeredCalorimeterData;
  caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::BarrelLayout;
  caloDetElem.addExtension<dd4hep::rec::LayeredCalorimeterData>(caloData);

  caloData->extent[0] = sensitiveBarrelRmin;
  caloData->extent[1] = sensitiveBarrelRmax;
  caloData->extent[2] = 0.;      // NN: for barrel detectors this is 0
  caloData->extent[3] = dzDetector; 

  dd4hep::rec::LayeredCalorimeterData::Layer caloLayer;

  for (unsigned int idxLayer = 0; idxLayer < layerDepths.size(); ++idxLayer) {
        const double difference_bet_r1r2 = layerDepths.at(idxLayer); 

        caloLayer.distance                  = layerInnerRadii.at(idxLayer); // radius of the current layer   
        caloLayer.sensitive_thickness       = difference_bet_r1r2;  // radial dimension of the current layer 
        caloLayer.inner_thickness           = difference_bet_r1r2 / 2.0;
        caloLayer.outer_thickness           = difference_bet_r1r2 / 2.0;

        caloData->layers.push_back(caloLayer);
  }


  return caloDetElem;

}
}  // namespace hcal

DECLARE_DETELEMENT(HCalTileBarrel_o1_v02, det::createHCal)
