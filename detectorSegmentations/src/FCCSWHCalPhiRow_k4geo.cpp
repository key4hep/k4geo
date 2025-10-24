#include "detectorSegmentations/FCCSWHCalPhiRow_k4geo.h"
#include "DD4hep/Printout.h"

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  FCCSWHCalPhiRow_k4geo::FCCSWHCalPhiRow_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) {
    commonSetup();
  }

  FCCSWHCalPhiRow_k4geo::FCCSWHCalPhiRow_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) { commonSetup(); }

  /// Initialization common to all ctors.
  void FCCSWHCalPhiRow_k4geo::commonSetup() {
    // define type and description
    _type = "FCCSWHCalPhiRow_k4geo";
    _description = "Phi-theta segmentation in the global coordinates";

    // register all necessary parameters
    registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
    registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
    registerParameter("grid_size_row", "Number of rows combined in a cell", m_gridSizeRow, std::vector<int>());
    registerParameter("dz_row", "dz of row", m_dz_row, 0.);
    registerParameter("detLayout", "The detector layout (0 = Barrel; 1 = Endcap)", m_detLayout, -1);
    registerParameter("offset_z", "Offset in z-axis of the layer center", m_offsetZ, std::vector<double>());
    registerParameter("width_z", "Width in z of the layer", m_widthZ, std::vector<double>());
    registerParameter("offset_r", "Offset in radius of the layer (Rmin)", m_offsetR, std::vector<double>());
    registerParameter("numLayers", "Number of layers", m_numLayers, std::vector<int>());
    registerParameter("dRlayer", "dR of the layer", m_dRlayer, std::vector<double>());
    registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
    registerIdentifier("identifier_row", "Cell ID identifier for row", m_rowID, "row");
    registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");

    m_layerIndex = decoder()->index(m_layerID);
    m_rowIndex = decoder()->index(m_rowID);
    m_phiIndex = decoder()->index(m_phiID);

    // Only endcap has "type" --- but it's too early to look at m_detLayout.
    for (const dd4hep::DDSegmentation::BitFieldElement& bfe : decoder()->fields()) {
      if (bfe.name() == "type") {
        m_typeIndex = decoder()->index("type");
        break;
      }
    }
  }

  /// determine the global position based on the cell ID
  Vector3D FCCSWHCalPhiRow_k4geo::position(const CellID& cID) const {
    uint layer = decoder()->get(cID, m_layerIndex);

    const LayerInfo& li = getLayerInfo(layer);

    double radius = li.radius;
    double minLayerZ = li.zmin;

    // get index of the cell in the layer (index starts from 1!)
    int idx = decoder()->get(cID, m_rowIndex);
    // calculate z-coordinate of the cell center
    double zpos = minLayerZ + (idx - 1) * m_dz_row * m_gridSizeRow[layer] + 0.5 * m_dz_row * m_gridSizeRow[layer];

    // for negative-z Endcap, the index is negative (starts from -1!)
    if (idx < 0)
      zpos = -minLayerZ + (idx + 1) * m_dz_row * m_gridSizeRow[layer] - 0.5 * m_dz_row * m_gridSizeRow[layer];

    return Vector3D(radius * std::cos(phi(cID)), radius * std::sin(phi(cID)), zpos);
  }

  auto FCCSWHCalPhiRow_k4geo::getLayerInfo(const unsigned layer) const -> const LayerInfo& {
    // If the LayerInfo vector hasn't been made yet, calculate it now.
    const std::vector<LayerInfo>* liv = m_layerInfo.load();
    if (!liv) {
      auto liv_new = new std::vector<LayerInfo>(initLayerInfo());
      if (m_layerInfo.compare_exchange_strong(liv, liv_new)) {
        liv = liv_new;
      } else {
        delete liv_new;
      }
    }

    return liv->at(layer);
  }

  // Initialize derived derived layer information.
  std::vector<FCCSWHCalPhiRow_k4geo::LayerInfo> FCCSWHCalPhiRow_k4geo::initLayerInfo() const {
    std::vector<LayerInfo> out;
    if (!checkParameters()) {
      out.resize((*decoder())[m_layerIndex].maxValue() + 1);
      return out;
    }

    if (m_detLayout == 0)
      dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiRow_k4geo", "Barrel configuration found!");
    else
      dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiRow_k4geo", "EndCap configuration found!");

    // calculate the radius for each layer
    uint N_dR = m_numLayers.size() / m_offsetZ.size();
    std::vector<double> moduleDepth = m_offsetR;
    for (uint i_section = 0; i_section < m_offsetZ.size(); i_section++) {
      // lower and upper edges in z-axis
      double zmin = m_offsetZ[i_section] - 0.5 * m_widthZ[i_section];
      double zmax = m_offsetZ[i_section] + 0.5 * m_widthZ[i_section];

      // Loop over groups of layers.
      for (uint i_dR = 0; i_dR < N_dR; i_dR++) {
        // Loop over individual layers.
        for (int i_lay = 0; i_lay < m_numLayers[i_dR + i_section * N_dR]; i_lay++) {
          moduleDepth[i_section] += m_dRlayer[i_dR];
          out.push_back(LayerInfo{.radius = moduleDepth[i_section] - m_dRlayer[i_dR] * 0.5,
                                  .halfDepth = m_dRlayer[i_dR] / 2,
                                  .zmin = zmin,
                                  .zmax = zmax});
        }
      }
    }

    // print info of calculated radii and edges
    for (uint i_layer = 0; const LayerInfo& li : out) {
      dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiRow_k4geo", "layer %d radius: %.2f, z range: %.2f - %.2f cm",
                       i_layer++, li.radius, li.zmin, li.zmax);
    }

    // determine theta bins and cell edges for each layer
    for (uint i_layer = 0; LayerInfo& li : out) {
      defineCellIndexes(li, i_layer);
      ++i_layer;
    }
    return out;
  }

  /*
   *  This function fills the cellIndexes vector per layer with the cell indexes.
   *  The cell index is encoded in CellID with "row" field.
   *  In case of a cell with single row/sequence, the index is directly the number of row in the layer.
   *  In case of a cell with several rows/sequences merged, the index is the number of cell in the layer.
   *  For the layers of negative-z Endcap, indexes of cells are negative (-1, -2, ..., -N).
   *
   *  cellIndexes vector can be used to define the neighbours using CreateFCCeeCaloNeighbours tool.
   */
  void FCCSWHCalPhiRow_k4geo::defineCellIndexes(LayerInfo& li, const uint layer) const {
    double minLayerZ = li.zmin;
    double maxLayerZ = li.zmax;

    // find rows/sequences that fit within the given layer range along z
    int irow = 0;
    while ((minLayerZ + (irow + 1) * m_dz_row) < (maxLayerZ + 0.0001)) {
      // define the cell index
      int idx = floor(irow / m_gridSizeRow[layer]) + 1;
      // add the index if it is not already there
      if (std::find(li.cellIndexes.begin(), li.cellIndexes.end(), idx) == li.cellIndexes.end())
        li.cellIndexes.push_back(idx);
      irow++;
    }

    // for the EndCap, do it again but for negative-z part
    if (m_detLayout == 1) {
      irow = 0;
      while ((minLayerZ + (irow + 1) * m_dz_row) < (maxLayerZ + 0.0001)) {
        // define the cell index with negative sign
        int idx = -(floor(irow / m_gridSizeRow[layer]) + 1);
        // add the index if it is not already there
        if (std::find(li.cellIndexes.begin(), li.cellIndexes.end(), idx) == li.cellIndexes.end())
          li.cellIndexes.push_back(idx);
        irow++;
      }
    }

    // find edges of each cell in the given layer along z axis
    for (auto idx : li.cellIndexes) {
      // calculate z-coordinates of the cell edges
      double z1 = minLayerZ + (idx - 1) * m_dz_row * m_gridSizeRow[layer]; // lower edge
      double z2 = minLayerZ + (idx)*m_dz_row * m_gridSizeRow[layer];       // upper edge

      // for negative-z Endcap, the index is negative (starts from -1!)
      if (idx < 0) {
        z1 = -minLayerZ + (idx)*m_dz_row * m_gridSizeRow[layer];       // lower edge
        z2 = -minLayerZ + (idx + 1) * m_dz_row * m_gridSizeRow[layer]; // upper edge
      }

      li.cellEdges[idx] = std::make_pair(z1, z2);
    }

    dd4hep::printout(dd4hep::DEBUG, "FCCSWHCalPhiRow_k4geo", "Number of cells in layer %d: %d", layer,
                     li.cellIndexes.size());
  }

  // Check consistency of input geometric variables.
  bool FCCSWHCalPhiRow_k4geo::checkParameters() const {
    // check if all necessary variables are available
    if (m_detLayout == -1 || m_offsetZ.empty() || m_widthZ.empty() || m_offsetR.empty() || m_numLayers.empty() ||
        m_dRlayer.empty()) {
      dd4hep::printout(
          dd4hep::ERROR, "FCCSWHCalPhiRow_k4geo", "Please check the readout description in the XML file!\n%s",
          "One of the variables is missing: detLayout | offset_z | width_z | offset_r | numLayers | dRlayer");
      return false;
    }

    // some sanity checks of the xml
    if (m_offsetZ.size() != m_offsetR.size()) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiRow_k4geo",
                       "Please check the readout description in the XML file!\n%s",
                       "Number of elements in offsetZ and offsetR must be the same!");
      return false;
    }

    if (m_widthZ.size() != m_offsetR.size()) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiRow_k4geo",
                       "Please check the readout description in the XML file!\n%s",
                       "Number of elements in widthZ and offsetR must be the same!");
      return false;
    }

    if (m_detLayout == 0 && m_offsetZ.size() != 1) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiRow_k4geo",
                       "Please check the readout description in the XML file!\n%s",
                       "Number of elements in offsetZ/offsetR/widthZ must be 1 for the Barrel!");
      return false;
    }

    if (m_numLayers.size() % m_offsetZ.size() != 0) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiRow_k4geo",
                       "Please check the readout description in the XML file!\n%s",
                       "Number of elements in numLayers must be multiple of offsetZ.size()!");
      return false;
    }

    if (m_dRlayer.size() != m_numLayers.size() / m_offsetZ.size()) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiRow_k4geo",
                       "Please check the readout description in the XML file!\n%s",
                       "Number of elements in dRlayer must be equal to numLayers.size()/offsetZ.size()!");
      return false;
    }

    uint nlayers = 0;
    for (auto n : m_numLayers)
      nlayers += n;
    if (m_gridSizeRow.size() != nlayers) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiRow_k4geo",
                       "Please check the readout description in the XML file!\n%s",
                       "Number of elements in gridSizeRow must be equal to sum of contents of numLayers!");
      return false;
    }

    return true;
  }

  /// create the cell ID based on the position
  CellID FCCSWHCalPhiRow_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                                       const VolumeID& vID) const {

    // get the row number from volumeID (starts from 0!)
    int nrow = decoder()->get(vID, m_rowIndex);
    // get the layer number from volumeID
    uint layer = decoder()->get(vID, m_layerIndex);

    CellID cID = vID;

    // get the cell index (start from 1!)
    int idx = floor(nrow / m_gridSizeRow[layer]) + 1;

    // if the hit is in the negative-z part of the Endcap then assign negative index
    if (m_detLayout == 1 && globalPosition.z() < 0)
      idx *= -1;

    decoder()->set(cID, m_rowIndex, idx);
    decoder()->set(cID, m_phiIndex,
                   positionToBin(dd4hep::DDSegmentation::Util::phiFromXYZ(globalPosition), 2 * M_PI / (double)m_phiBins,
                                 m_offsetPhi));

    // For endcap, the volume ID comes with "type" field information which would screw up the topo-clustering,
    // therefore, lets set it to zero, as it is for the cell IDs in the neighbours map.
    if (m_detLayout == 1)
      decoder()->set(cID, m_typeIndex, 0);

    return cID;
  }

  /// determine the azimuthal angle phi based on the cell ID
  double FCCSWHCalPhiRow_k4geo::phi(const CellID cID) const {
    CellID phiValue = decoder()->get(cID, m_phiIndex);
    return binToPosition(phiValue, 2. * M_PI / (double)m_phiBins, m_offsetPhi);
  }

  /// Get the min and max layer indexes for each part of the HCal
  std::vector<std::pair<uint, uint>> FCCSWHCalPhiRow_k4geo::getMinMaxLayerId() const {

    std::vector<std::pair<uint, uint>> minMaxLayerId;

    std::vector<uint> minLayerId(m_offsetZ.size(), 0);
    std::vector<uint> maxLayerId(m_offsetZ.size(), 0);

    uint Nl = m_numLayers.size() / m_offsetZ.size();

    for (uint i_section = 0; i_section < m_offsetZ.size(); i_section++) {
      if (i_section > 0) {
        minLayerId[i_section] = maxLayerId[i_section - 1] + 1;
        maxLayerId[i_section] = maxLayerId[i_section - 1];
      }

      for (uint i = 0; i < Nl; i++) {
        maxLayerId[i_section] += m_numLayers[i + i_section * Nl];
      }

      if (i_section == 0)
        maxLayerId[0] -= 1;

      minMaxLayerId.push_back(std::make_pair(minLayerId[i_section], maxLayerId[i_section]));
    }

    return minMaxLayerId;
  }

  /// Calculates the neighbours of the given cell ID and adds them to the list of neighbours
  std::vector<uint64_t> FCCSWHCalPhiRow_k4geo::neighbours(const CellID cID) const {
    std::vector<uint64_t> cellNeighbours;

    uint EndcapPart = 0;
    int minLayerId = -1;
    int maxLayerId = -1;

    int currentLayerId = decoder()->get(cID, m_layerIndex);
    int currentCellId = decoder()->get(cID, m_rowIndex);

    const LayerInfo& li = getLayerInfo(currentLayerId);

    int minCellId = li.cellIndexes.front();
    int maxCellId = li.cellIndexes.back();

    //--------------------------------
    // Determine min and max layer Id
    //--------------------------------
    std::vector<int> minLayerIdEndcap(m_offsetZ.size(), 0);
    std::vector<int> maxLayerIdEndcap(m_offsetZ.size(), 0);
    if (m_detLayout == 1) // Endcap
    {
      std::vector<std::pair<uint, uint>> minMaxLayerId(getMinMaxLayerId());
      if (minMaxLayerId.empty()) {
        dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiRow_k4geo",
                         "Can not get Endcap min and max layer indexes! --> returning empty neighbours");
        return cellNeighbours;
      }

      // determine min and max layer Ids and which section/part it is
      for (uint i_section = 0; i_section < minMaxLayerId.size(); i_section++) {
        minLayerIdEndcap[i_section] = minMaxLayerId[i_section].first;
        maxLayerIdEndcap[i_section] = minMaxLayerId[i_section].second;

        if (currentLayerId >= minLayerIdEndcap[i_section] && currentLayerId <= maxLayerIdEndcap[i_section]) {
          minLayerId = minLayerIdEndcap[i_section];
          maxLayerId = maxLayerIdEndcap[i_section];
          EndcapPart = i_section;
        }
      }

      // correct the min and max CellId for endcap
      if (currentCellId < 0) // negative-z part (cell index is negative (-1, -2, ... -N))
      {
        // second half of elements in li.cellIndexes vector corresponds to the negative-z layer cells
        minCellId = li.cellIndexes.back();
        maxCellId = li.cellIndexes[li.cellIndexes.size() / 2];
      } else // positive-z part (cell index is positive (1, 2, ... N))
      {
        // first half of elements in li.cellIndexes vector corresponds to the positive-z layer cells
        minCellId = li.cellIndexes.front();
        maxCellId = li.cellIndexes[li.cellIndexes.size() / 2 - 1];
      }
    } else // for Barrel
    {
      minLayerId = 0;
      maxLayerId = m_layerInfo.load()->size() - 1;
    }
    //--------------------------------

    //--------------------------------
    // Find neighbours
    //--------------------------------

    // if this is not the first layer then look for neighbours in the previous layer
    if (currentLayerId > minLayerId) {
      CellID nID = cID;
      int prevLayerId = currentLayerId - 1;
      decoder()->set(nID, m_layerIndex, prevLayerId);

      // if the granularity is the same for the previous layer then take the cells with currentCellId, currentCellId -
      // 1, and currentCellId + 1
      if (m_gridSizeRow[prevLayerId] == m_gridSizeRow[currentLayerId]) {
        decoder()->set(nID, m_rowIndex, currentCellId);
        cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
        if (currentCellId > minCellId) {
          decoder()->set(nID, m_rowIndex, currentCellId - 1);
          cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
        }
        if (currentCellId < maxCellId) {
          decoder()->set(nID, m_rowIndex, currentCellId + 1);
          cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
        }
      }
      // if the granularity is different
      else {
        // determine the cell index in the previous layer that is below of the current cell
        int idx = (currentCellId > 0) ? ((currentCellId - 1) / m_gridSizeRow[prevLayerId] + 1)
                                      : ((currentCellId + 1) / m_gridSizeRow[prevLayerId] - 1);
        decoder()->set(nID, m_rowIndex, idx);
        cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module

        //
        if ((m_gridSizeRow[prevLayerId] - abs(currentCellId) % m_gridSizeRow[prevLayerId]) ==
                (m_gridSizeRow[prevLayerId] - 1) &&
            currentCellId > minCellId) {
          decoder()->set(nID, m_rowIndex, (idx > 0) ? (idx - 1) : (idx + 1));
          cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
        }

        //
        if (abs(currentCellId) % m_gridSizeRow[prevLayerId] == 0 && currentCellId < maxCellId) {
          decoder()->set(nID, m_rowIndex, (idx > 0) ? (idx + 1) : (idx - 1));
          cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
        }
      }
    }

    // if this is not the last layer then look for neighbours in the next layer
    if (currentLayerId < maxLayerId) {
      CellID nID = cID;
      int nextLayerId = currentLayerId + 1;
      decoder()->set(nID, m_layerIndex, nextLayerId);

      // if the granularity is the same for the next layer then take the cells with currentCellId, currentCellId - 1,
      // and currentCellId + 1
      if (m_gridSizeRow[nextLayerId] == m_gridSizeRow[currentLayerId]) {
        decoder()->set(nID, m_rowIndex, currentCellId);
        cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
        if (currentCellId > minCellId) {
          decoder()->set(nID, m_rowIndex, currentCellId - 1);
          cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
        }
        if (currentCellId < maxCellId) {
          decoder()->set(nID, m_rowIndex, currentCellId + 1);
          cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
        }
      }
      // if the granularity is different
      else {
        // determine the cell index in the next layer that is below of the current cell
        int idx = (currentCellId > 0) ? ((currentCellId - 1) / m_gridSizeRow[nextLayerId] + 1)
                                      : ((currentCellId + 1) / m_gridSizeRow[nextLayerId] - 1);
        decoder()->set(nID, m_rowIndex, idx);
        cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module

        //
        if ((m_gridSizeRow[nextLayerId] - abs(currentCellId) % m_gridSizeRow[nextLayerId]) ==
                (m_gridSizeRow[nextLayerId] - 1) &&
            currentCellId > minCellId) {
          decoder()->set(nID, m_rowIndex, (idx > 0) ? (idx - 1) : (idx + 1));
          cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
        }

        //
        if (abs(currentCellId) % m_gridSizeRow[nextLayerId] == 0 && currentCellId < maxCellId) {
          decoder()->set(nID, m_rowIndex, (idx > 0) ? (idx + 1) : (idx - 1));
          cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
        }
      }
    }

    // if this is not the first cell in the given layer then add the previous cell
    if (currentCellId > minCellId) {
      CellID nID = cID;
      decoder()->set(nID, m_rowIndex, currentCellId - 1);
      cellNeighbours.push_back(nID); // add the previous cell from current layer of the same phi module
    }
    // if this is not the last cell in the given layer then add the next cell
    if (currentCellId < maxCellId) {
      CellID nID = cID;
      decoder()->set(nID, m_rowIndex, currentCellId + 1);
      cellNeighbours.push_back(nID); // add the next cell from current layer of the same phi module
    }

    // if this is the Endcap (and consists of more than 1 part/section) then look for neighbours in different parts as
    // well
    if (m_detLayout == 1 && m_offsetZ.size() > 1) {
      double currentLayerRmin = li.radius - li.halfDepth;
      double currentLayerRmax = li.radius + li.halfDepth;

      // if the cell is in negative-z part, then swap min and max cell indexes
      if (currentCellId < 0) {
        minCellId = li.cellIndexes[li.cellIndexes.size() / 2]; // this should be -1
        maxCellId = li.cellIndexes.back();                     // this should be -N
      }

      // if it is the last cell in the first part
      if (EndcapPart == 0 && currentCellId == maxCellId) {
        // find the layers in the part2 that share a border with the current layer
        for (int part2layerId = minLayerIdEndcap[1]; part2layerId <= maxLayerIdEndcap[1]; part2layerId++) {
          const LayerInfo& part2li = getLayerInfo(part2layerId);
          double Rmin = part2li.radius - part2li.halfDepth;
          double Rmax = part2li.radius + part2li.halfDepth;

          if ((Rmin >= currentLayerRmin && Rmin <= currentLayerRmax) ||
              (Rmax >= currentLayerRmin && Rmax <= currentLayerRmax) ||
              (currentLayerRmin >= Rmin && currentLayerRmin <= Rmax) ||
              (currentLayerRmax >= Rmin && currentLayerRmax <= Rmax)) {
            CellID nID = cID;
            decoder()->set(nID, m_layerIndex, part2layerId);
            decoder()->set(nID, m_rowIndex, minCellId);
            cellNeighbours.push_back(nID); // add the first cell from part2 layer
          }
        }
      }

      // if the Endcap consists of more than 2 parts:
      for (uint i_section = 1; i_section < (m_offsetZ.size() - 1); i_section++) {
        if (i_section != EndcapPart)
          continue;

        // if it is the first cell then look for neighbours in previous part
        if (currentCellId == minCellId) {
          // find the layers in previous part that share a border with the current layer
          for (int prevPartLayerId = minLayerIdEndcap[i_section - 1];
               prevPartLayerId <= maxLayerIdEndcap[i_section - 1]; prevPartLayerId++) {
            const LayerInfo& prevPartli = getLayerInfo(prevPartLayerId);
            double Rmin = prevPartli.radius - prevPartli.halfDepth;
            double Rmax = prevPartli.radius + prevPartli.halfDepth;

            if ((Rmin >= currentLayerRmin && Rmin <= currentLayerRmax) ||
                (Rmax >= currentLayerRmin && Rmax <= currentLayerRmax) ||
                (currentLayerRmin >= Rmin && currentLayerRmin <= Rmax) ||
                (currentLayerRmax >= Rmin && currentLayerRmax <= Rmax)) {
              CellID nID = cID;
              decoder()->set(nID, m_layerIndex, prevPartLayerId);
              decoder()->set(nID, m_rowIndex, maxCellId);
              cellNeighbours.push_back(nID); // add the last cell from the previous part layer
            }
          }
        }
        // if it is the last cell then look for neighbours in the next part
        if (currentCellId == maxCellId) {
          // find the layers in the next that share a border with the current layer
          for (int nextPartLayerId = minLayerIdEndcap[i_section + 1];
               nextPartLayerId <= maxLayerIdEndcap[i_section + 1]; nextPartLayerId++) {
            const LayerInfo& nextPartli = getLayerInfo(nextPartLayerId);
            double Rmin = nextPartli.radius - nextPartli.halfDepth;
            double Rmax = nextPartli.radius + nextPartli.halfDepth;

            if ((Rmin >= currentLayerRmin && Rmin <= currentLayerRmax) ||
                (Rmax >= currentLayerRmin && Rmax <= currentLayerRmax) ||
                (currentLayerRmin >= Rmin && currentLayerRmin <= Rmax) ||
                (currentLayerRmax >= Rmin && currentLayerRmax <= Rmax)) {
              CellID nID = cID;
              decoder()->set(nID, m_layerIndex, nextPartLayerId);
              decoder()->set(nID, m_rowIndex, minCellId);
              cellNeighbours.push_back(nID); // add the first cell from the next part layer
            }
          }
        }
      }

      // if it is the first cell in the last part of the Endcap
      if (EndcapPart == (m_offsetZ.size() - 1) && currentCellId == minCellId) {
        // find the layers in the previous part that share a border with the current layer
        for (int prevPartLayerId = minLayerIdEndcap[m_offsetZ.size() - 2];
             prevPartLayerId <= maxLayerIdEndcap[m_offsetZ.size() - 2]; prevPartLayerId++) {
          const LayerInfo& prevPartli = getLayerInfo(prevPartLayerId);
          double Rmin = prevPartli.radius - prevPartli.halfDepth;
          double Rmax = prevPartli.radius + prevPartli.halfDepth;

          if ((Rmin >= currentLayerRmin && Rmin <= currentLayerRmax) ||
              (Rmax >= currentLayerRmin && Rmax <= currentLayerRmax) ||
              (currentLayerRmin >= Rmin && currentLayerRmin <= Rmax) ||
              (currentLayerRmax >= Rmin && currentLayerRmax <= Rmax)) {
            CellID nID = cID;
            decoder()->set(nID, m_layerIndex, prevPartLayerId);
            decoder()->set(nID, m_rowIndex, maxCellId);
            cellNeighbours.push_back(nID); // add the last cell from the part2 layer
          }
        }
      }
    }

    // Now loop over the neighbours and add the cells from next/previous phi module
    std::vector<uint64_t> cellNeighboursCopy(cellNeighbours);
    for (auto nID : cellNeighboursCopy) {
      CellID newID = nID;
      // previous: if the current is 0 then previous is the last bin (id = m_phiBins - 1) else current - 1
      decoder()->set(newID, m_phiIndex,
                     (decoder()->get(nID, m_phiIndex) == 0) ? m_phiBins - 1 : decoder()->get(nID, m_phiIndex) - 1);
      cellNeighbours.push_back(newID);
      // next: if the current is the last bin (id = m_phiBins - 1) then the next is the first bin (id = 0) else current
      // + 1
      decoder()->set(newID, m_phiIndex,
                     (decoder()->get(nID, m_phiIndex) == (m_phiBins - 1)) ? 0 : decoder()->get(nID, m_phiIndex) + 1);
      cellNeighbours.push_back(newID);
    }

    // At the end, find neighbours with the same layer/row in next/previous phi module
    CellID nID = cID;
    // previous: if the current is 0 then previous is the last bin (id = m_phiBins - 1) else current - 1
    decoder()->set(nID, m_phiIndex,
                   (decoder()->get(cID, m_phiIndex) == 0) ? m_phiBins - 1 : decoder()->get(cID, m_phiIndex) - 1);
    cellNeighbours.push_back(nID);
    // next: if the current is the last bin (id = m_phiBins - 1) then the next is the first bin (id = 0) else current +
    // 1
    decoder()->set(nID, m_phiIndex,
                   (decoder()->get(cID, m_phiIndex) == (m_phiBins - 1)) ? 0 : decoder()->get(cID, m_phiIndex) + 1);
    cellNeighbours.push_back(nID);

    return cellNeighbours;
  }

  // Implement the signature from the Segmentations base class.
  void FCCSWHCalPhiRow_k4geo::neighbours(const CellID& cellID, std::set<CellID>& neighbours) const {
    std::vector<uint64_t> neigh = this->neighbours(cellID);
    neighbours.clear();
    neighbours.insert(neigh.begin(), neigh.end());
  }

  /// Determine minimum and maximum polar angle of the cell
  std::array<double, 2> FCCSWHCalPhiRow_k4geo::cellTheta(const CellID cID) const {
    std::array<double, 2> cTheta = {M_PI, M_PI};

    // get the cell index
    int idx = decoder()->get(cID, m_rowIndex);
    // get the layer index
    uint layer = decoder()->get(cID, m_layerIndex);

    const LayerInfo& li = getLayerInfo(layer);

    const auto& edges = li.cellEdges.find(idx)->second;
    double zlow = edges.first;
    double zhigh = edges.second;

    double Rmin = li.radius - li.halfDepth;
    double Rmax = li.radius + li.halfDepth;

    if (theta(cID) < M_PI / 2.) {
      cTheta[0] = std::atan2(Rmin, zhigh); // theta min
      cTheta[1] = std::atan2(Rmax, zlow);  // theta max
    } else {
      cTheta[0] = std::atan2(Rmax, zhigh); // theta min
      cTheta[1] = std::atan2(Rmin, zlow);  // theta max
    }

    return cTheta;
  }

  /// determine maximum theta value of the detector. This is used by SW clustering
  double FCCSWHCalPhiRow_k4geo::thetaMax() const {
    std::vector<std::pair<uint, uint>> minMaxLayerId(getMinMaxLayerId());
    if (minMaxLayerId.empty())
      return 0.;

    // get the first layerId in the Barrel or in the last part of the Endcap
    uint layer = minMaxLayerId[minMaxLayerId.size() - 1].first;

    const LayerInfo& li = getLayerInfo(layer);

    // get the last cell index (which is in the positive-z side)
    int idx = abs(li.cellIndexes.back());

    // get the z-coordinate of the right-hand edge of the last cell
    const auto& edges = li.cellEdges.find(idx)->second;
    double zhigh = edges.second;

    // get the inner radius of the first layer
    double Rmin = li.radius - 0.5 * li.halfDepth;

    // calculate the minimum theta of the last cell in the first layer -> this is the minimum theta of the detector
    // (Barrel or Endcap)
    double thetaMin = std::atan2(Rmin, zhigh); // theta min

    return (M_PI - thetaMin); // theta max
  }

} // namespace DDSegmentation
} // namespace dd4hep
