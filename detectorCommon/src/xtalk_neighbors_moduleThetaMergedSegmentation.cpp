#include "detectorCommon/xtalk_neighbors_moduleThetaMergedSegmentation.h"

#include <iostream>

namespace det {
namespace crosstalk {

  // use it for module-theta merged readout (FCCSWGridModuleThetaMerged_k4geo)
  std::vector<std::pair<uint64_t, double>>
  getNeighboursModuleThetaMerged(const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo& aSeg,
                                 const dd4hep::DDSegmentation::BitFieldCoder& aDecoder,
                                 const std::vector<std::string>& aFieldNames,
                                 const std::vector<std::vector<std::pair<int, int>>>& aFieldExtremes_layer,
                                 uint64_t aCellId, const double aXtalkCoefRadial, const double aXtalkCoefTheta,
                                 const double aXtalkCoefDiagonal, const double aXtalkCoefTower) {

    std::vector<std::pair<uint64_t, double>> xtalk_neighbours;

    // check that field names and extremes have the proper length
    if (aFieldNames.size() != 3 || aFieldExtremes_layer[0].size() != 3) {
      std::cout << "ERROR: the vectors aFieldNames and aFieldSizes should be of length = 3, corresponding to the "
                   "theta/module/layer fields"
                << std::endl;
      std::cout << "ERROR: will return empty neighbour map" << std::endl;
      return xtalk_neighbours;
    }
    // find index of layer, module and theta in the field name vector
    int idModuleField(-1);
    int idThetaField(-1);
    int idLayerField(-1);
    for (size_t itField = 0; itField < aFieldNames.size(); itField++) {
      if (aFieldNames[itField] == aSeg.fieldNameModule())
        idModuleField = itField;
      else if (aFieldNames[itField] == aSeg.fieldNameTheta())
        idThetaField = itField;
      else if (aFieldNames[itField] == aSeg.fieldNameLayer())
        idLayerField = itField;
    }
    if (idModuleField < 0) {
      std::cout << "WARNING: module field " << aSeg.fieldNameModule() << " not found in aFieldNames vector"
                << std::endl;
      std::cout << "WARNING: will return empty neighbour map" << std::endl;
      return xtalk_neighbours;
    }
    if (idThetaField < 0) {
      std::cout << "WARNING: theta field " << aSeg.fieldNameTheta() << " not found in aFieldNames vector" << std::endl;
      std::cout << "WARNING: will return empty neighbour map" << std::endl;
      return xtalk_neighbours;
    }
    if (idLayerField < 0) {
      std::cout << "WARNING: layer field " << aSeg.fieldNameLayer() << " not found in aFieldNames vector" << std::endl;
      std::cout << "WARNING: will return empty neighbour map" << std::endl;
      return xtalk_neighbours;
    }

    // retrieve layer/module/theta of cell under study
    int layer_id = aDecoder.get(aCellId, aFieldNames[idLayerField]);
    // int module_id = aDecoder.get(aCellId, aFieldNames[idModuleField]);
    int theta_id = aDecoder.get(aCellId, aFieldNames[idThetaField]);

    int layer_indice_offset = layer_id - aFieldExtremes_layer[0][idLayerField].first;

    // now find crosstalk neighbours
    dd4hep::DDSegmentation::CellID cID = aCellId;

    // for neighbours across different layers, we have to take into
    // account that merging along module and/or theta could be different
    // so one cell in layer N could be neighbour to several in layer N+-1
    // The cells are classified in a different way whether they are
    // direct neighbours (common surface), diagonal neighbours (common edge or vertex)
    // or neither.
    // To decide this, we need to check how the cells are related in both directions:
    // neighbours (edge at least partially in common), diagonal neigbours (common vertex),
    // none
    for (int deltaLayer = -1; deltaLayer < 2; deltaLayer += 2) {

      // no neighbours in layer N-1 for innermost layer
      if (layer_id == aFieldExtremes_layer[layer_indice_offset][idLayerField].first && deltaLayer < 0)
        continue;
      // and in layer N+1 for outermost layer
      if (layer_id == aFieldExtremes_layer[layer_indice_offset][idLayerField].second && deltaLayer > 0)
        continue;

      // set layer field of neighbour cell
      aDecoder.set(cID, aSeg.fieldNameLayer(), layer_id + deltaLayer);

      // find the neighbour(s) in module and theta
      // if the numbers of module (theta) merged cells across 2 layers are the
      // same then we just take the same module (theta) ID
      // otherwise, we need to do some math to account for the different mergings
      // note: even if number of merged cells in layer-1 is larger, a cell
      // in layer could neighbour more than one cell in layer-1 if the merged
      // cells are not aligned, for example if cells are grouped by 3 in a layer
      // and by 4 in the next one, cell 435 in the former (which groups together
      // 435-436-437) will be neighbour to cells 432 and 436 of the latter
      // this might introduce duplicates, we will remove them later
      // another issue is that it could add spurious cells beyond the maximum module number
      // to prevent this we would need to know the max module number in layer -1
      // which would require modifying this function passing the extrema for all layers
      // instead of the extrema only for a certain layer
      // this border effect is also present in the original method..
      for (int i = -1; i <= aSeg.mergedThetaCells(layer_id); i++) {
        double this_xtalk =
            0.0; // will be set to the crosstalk coefficient for this neighbour, either radial or diagonal
        int theta_id_neighbour = (theta_id + i) - ((theta_id + i) % (aSeg.mergedThetaCells(layer_id + deltaLayer)
                                                                         ? aSeg.mergedThetaCells(layer_id + deltaLayer)
                                                                         : 1));
        // Do we need to check if the index neighbour theta is valid?
        if ((theta_id_neighbour >= theta_id && theta_id_neighbour < (theta_id + aSeg.mergedThetaCells(layer_id))) or
            (theta_id_neighbour < theta_id &&
             theta_id_neighbour >
                 (theta_id - aSeg.mergedThetaCells(layer_id + deltaLayer)))) { // crosstalk neighbour type 1
          this_xtalk = aXtalkCoefRadial;
        } else if ((theta_id_neighbour == (theta_id + aSeg.mergedThetaCells(layer_id))) or
                   (theta_id_neighbour ==
                    (theta_id - aSeg.mergedThetaCells(layer_id + deltaLayer)))) { // crosstalk neighbour type 3
          this_xtalk = aXtalkCoefDiagonal;
        } else
          continue;

        // check if neighbour theta index is valid
        if (theta_id_neighbour < aFieldExtremes_layer[layer_indice_offset + deltaLayer][idThetaField].first ||
            theta_id_neighbour > aFieldExtremes_layer[layer_indice_offset + deltaLayer][idThetaField].second)
          continue;
        aDecoder.set(cID, aSeg.fieldNameTheta(), theta_id_neighbour);
        xtalk_neighbours.emplace_back(cID, this_xtalk);
      }
    }

    // cross talk neighbour type 2
    // for neighbours in module/theta direction at same layer_id, do +-nMergedCells instead of +-1
    // reset cellID
    aDecoder.set(cID, aSeg.fieldNameLayer(), layer_id);
    aDecoder.set(cID, aSeg.fieldNameTheta(), theta_id);
    // loop over theta cells
    for (int i = -1; i <= 1; i += 2) {
      // calculate theta_id of neighbour
      int theta_id_neighbour = theta_id + i * aSeg.mergedThetaCells(layer_id);
      // check if it is within the theta range of the layer
      if ((theta_id_neighbour < aFieldExtremes_layer[layer_indice_offset][idThetaField].first) ||
          (theta_id_neighbour > aFieldExtremes_layer[layer_indice_offset][idThetaField].second))
        continue;
      // set theta_id of cell ID
      aDecoder[aSeg.fieldNameTheta()].set(cID, theta_id + i * aSeg.mergedThetaCells(layer_id));
      xtalk_neighbours.emplace_back(cID, aXtalkCoefTheta);
    }

    // cross talk type 4
    // reset cellID
    aDecoder.set(cID, aSeg.fieldNameLayer(), layer_id);
    aDecoder.set(cID, aSeg.fieldNameTheta(), theta_id);
    for (int loopLayer = aFieldExtremes_layer[0][idLayerField].first;
         loopLayer <= aFieldExtremes_layer[0][idLayerField].second; loopLayer++) {
      // cells in adjacent layers have been studied, skip
      if ((loopLayer - layer_id) <= 1 && (loopLayer - layer_id) >= -1)
        continue;
      aDecoder.set(cID, aSeg.fieldNameLayer(), loopLayer);
      int loopLayer_indice_offset = loopLayer - aFieldExtremes_layer[0][idLayerField].first;
      // find crosstalk neighbours in the same theta tower
      for (int i = -1; i <= aSeg.mergedThetaCells(layer_id); i++) {
        int theta_id_neighbour =
            (theta_id + i) -
            ((theta_id + i) % (aSeg.mergedThetaCells(loopLayer) ? aSeg.mergedThetaCells(loopLayer) : 1));
        // check if the neighbour theta index is valid
        if (theta_id_neighbour >= theta_id && theta_id_neighbour < (theta_id + aSeg.mergedThetaCells(layer_id)) &&
            theta_id_neighbour >= aFieldExtremes_layer[loopLayer_indice_offset][idThetaField].first &&
            theta_id_neighbour <=
                aFieldExtremes_layer[loopLayer_indice_offset][idThetaField].second) { // crosstalk neighbour type 4
          aDecoder.set(cID, aSeg.fieldNameTheta(), theta_id_neighbour);
          xtalk_neighbours.emplace_back(cID, aXtalkCoefTower);
        } else if (theta_id_neighbour < theta_id &&
                   theta_id_neighbour > (theta_id - aSeg.mergedThetaCells(loopLayer)) &&
                   theta_id_neighbour >= aFieldExtremes_layer[loopLayer_indice_offset][idThetaField].first &&
                   theta_id_neighbour <= aFieldExtremes_layer[loopLayer_indice_offset][idThetaField]
                                             .second) { // crosstalk neighbour type 4
          aDecoder.set(cID, aSeg.fieldNameTheta(), theta_id_neighbour);
          xtalk_neighbours.emplace_back(cID, aXtalkCoefTower);
        } else
          continue;
      }
    }

    // remove duplicates
    std::vector<uint64_t> sort_neighbours;
    std::vector<std::pair<uint64_t, double>> new_neighbours;
    for (auto& i : xtalk_neighbours) {
      bool IsNewNeighbour = true;
      for (unsigned int i_loop = 0; i_loop < sort_neighbours.size(); i_loop++) {
        if (sort_neighbours[i_loop] == i.first) {
          IsNewNeighbour = false;
          break;
        }
      }
      if (IsNewNeighbour) {
        new_neighbours.emplace_back(i);
        sort_neighbours.emplace_back(i.first);
      }
    }
    xtalk_neighbours = new_neighbours;

    return xtalk_neighbours;
  }

  // return indices of layer/module/theta fields of a give cell. This is a function to be used for debug purpose
  std::vector<int> getCellIndices(const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo& aSeg,
                                  const dd4hep::DDSegmentation::BitFieldCoder& aDecoder,
                                  const std::vector<std::string>& aFieldNames, uint64_t aCellId) {

    // find index of layer, module and theta in the field names vector
    int idModuleField(-1);
    int idThetaField(-1);
    int idLayerField(-1);
    for (size_t itField = 0; itField < aFieldNames.size(); itField++) {
      if (aFieldNames[itField] == aSeg.fieldNameModule())
        idModuleField = itField;
      else if (aFieldNames[itField] == aSeg.fieldNameTheta())
        idThetaField = itField;
      else if (aFieldNames[itField] == aSeg.fieldNameLayer())
        idLayerField = itField;
    }
    // retrieve layer/module/theta of the cell under study
    int layer_id = aDecoder.get(aCellId, aFieldNames[idLayerField]);
    int module_id = aDecoder.get(aCellId, aFieldNames[idModuleField]);
    int theta_id = aDecoder.get(aCellId, aFieldNames[idThetaField]);
    return {layer_id, module_id, theta_id};
  }

} // namespace crosstalk
} // namespace det
