/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.

// Author     : S.Sasikumar
//
//==========================================================================
//
// Pandora Layered Calorimeter Plugin
//
// The plugin builds LayeredCalorimeterData dynamically by
// configuring detector's material properties, interaction lengths, and type classifications
// ensuring proper adaptibility of a calorimeter with PandoraPFA irrespective of its geometry.
//
//==========================================================================
*/

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Detector.h>
#include <DD4hep/Handle.h>
#include <DD4hep/Printout.h>
#include <DD4hep/DetElement.h>
#include <XML/Utilities.h>
#include <DDRec/MaterialManager.h>
#include <DD4hep/Factories.h>
#include <DDRec/Vector3D.h>
#include <DDRec/DetectorData.h>

#include <string>
#include <vector>
#include <stdexcept>

using namespace dd4hep;

namespace {

  static long addLayeredCalorimeterData(dd4hep::Detector& description, int argc, char** argv) {
    const std::string LOG_SOURCE("LayeredCalorimeterPlugin");

    // Ensure there are at least 4 arguments
    if (argc < 7) {
      printout(PrintLevel::ERROR, LOG_SOURCE, "Insufficient arguments. Expected: path, Rmin, Rmax, half_length,angle, and layer heights.");
      return 1;
    }

    // Parse basic parameters
    std::string subDetectorName = argv[0];
    std::string type_flags = argv[1];
    double Rmin = dd4hep::_toDouble(argv[2]);
    double Rmax = dd4hep::_toDouble(argv[3]);
    double half_length = dd4hep::_toDouble(argv[4]);
    double angle = dd4hep::_toDouble(argv[5]);
    std::string name_of_the_material = argv[6];

    // Parse layer heights (all remaining arguments)
    std::vector<double> layerHeights;
    for (int i = 7; i < argc; ++i) {
      layerHeights.push_back(dd4hep::_toDouble(argv[i]));
    }

    dd4hep::printout(PrintLevel::DEBUG, LOG_SOURCE, "Parsed: path=%s, Rmin=%f, Rmax=%f, dR=%f, angle=%f,layers=%zu",
                    subDetectorName.c_str(), Rmin, Rmax, half_length, angle,layerHeights.size());

    // Retrieve the detector element
    auto caloDetElem = description.detector(subDetectorName);
    if (!caloDetElem.isValid()) {
      dd4hep::printout(PrintLevel::ERROR, LOG_SOURCE, "Invalid DetElement path: %s", subDetectorName.c_str());
      return 1;
    }
     
    unsigned int typeFlag = dd4hep::_toInt(type_flags);
    dd4hep::printout(PrintLevel::INFO, LOG_SOURCE,"+++ setDetectorTypeFlags for detector :%s set to 0x%x",
		     subDetectorName.c_str(),
		     typeFlag);

    caloDetElem.setTypeFlag(typeFlag);


    // Create LayeredCalorimeterData
    auto caloData = new dd4hep::rec::LayeredCalorimeterData;
    caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::BarrelLayout;
    caloData->extent[0] = Rmin;
    caloData->extent[1] = Rmax;
    caloData->extent[2] = 0.0; // Barrel-specific value
    caloData->extent[3] = half_length;

    caloDetElem.addExtension<dd4hep::rec::LayeredCalorimeterData>(caloData);

    // Material manager
    rec::MaterialManager matMgr(caloDetElem.volume());
    double rad_first = Rmin;
    double rad_last = 0;
    double dR =  Rmax - Rmin;
    double scale_fact = dR / (-Rmin * cos(angle) + sqrt(pow(Rmax, 2) - pow(Rmin * sin(angle), 2)));

    // Loop through provided layer heights
    for (double layerHeight : layerHeights) {
      rad_last = rad_first + (layerHeight * scale_fact);
      dd4hep::rec::Vector3D ivr1 = dd4hep::rec::Vector3D(0., rad_first, 0); // defining starting vector points of the given layer
      dd4hep::rec::Vector3D ivr2 = dd4hep::rec::Vector3D(0., rad_last, 0);  // defining end vector points of the given layer
     
      const dd4hep::rec::MaterialVec &materials = matMgr.materialsBetween(ivr1, ivr2); // calling material manager(MM) to get material info between two points
      auto mat = matMgr.createAveragedMaterial(materials);                             // creating average of all the material between two points to calculate X0 and lambda of averaged material
      const double nRadiationLengths = mat.radiationLength();
      const double nInteractionLengths = mat.interactionLength();
      const double difference_bet_r1r2 = (ivr1 - ivr2).r();
      const double value_of_x0 = layerHeight / nRadiationLengths;
      const double value_of_lambda = layerHeight / nInteractionLengths;

      double thickness_sen = 0.;
      double absorberThickness = 0.;
      for (size_t imat = 0; imat < materials.size(); imat++) {

	std::string material_found_by_MM(materials.at(imat).first.name());
	if (name_of_the_material.compare(material_found_by_MM) == 0){
	  thickness_sen += materials.at(imat).second;
	}
	else {
	  absorberThickness += materials.at(imat).second;
	}
      }
      rad_first = rad_last;
     

      // Define and add the layer
      rec::LayeredCalorimeterData::Layer caloLayer;
      caloLayer.distance = rad_first;
      caloLayer.sensitive_thickness       = thickness_sen;
      caloLayer.inner_nRadiationLengths   = value_of_x0 / 2.0;
      caloLayer.inner_nInteractionLengths = value_of_lambda / 2.0;
      caloLayer.inner_thickness           = difference_bet_r1r2 / 2.0;

      caloLayer.outer_nRadiationLengths   = value_of_x0 / 2.0;
      caloLayer.outer_nInteractionLengths = value_of_lambda / 2.0;
      caloLayer.outer_thickness           = difference_bet_r1r2 / 2;
      caloLayer.absorberThickness         = absorberThickness;
      
      caloLayer.cellSize0 = 20 * dd4hep::mm;
      caloLayer.cellSize1 = 20 * dd4hep::mm;
     

      caloData->layers.push_back(caloLayer);
     
    }

    dd4hep::printout(PrintLevel::INFO, LOG_SOURCE, "Successfully added LayeredCalorimeterData to %s", subDetectorName.c_str());
    return 1;
  }
}// namespace

DECLARE_APPLY(Pandora_LayeredCalorimeterPlugin, ::addLayeredCalorimeterData)

