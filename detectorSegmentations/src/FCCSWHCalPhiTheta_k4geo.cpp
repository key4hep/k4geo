#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"
#include "DD4hep/Printout.h"

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  FCCSWHCalPhiTheta_k4geo::FCCSWHCalPhiTheta_k4geo(const std::string& cellEncoding) : GridTheta_k4geo(cellEncoding) {
    // define type and description
    _type = "FCCSWHCalPhiTheta_k4geo";
    _description = "Phi-theta segmentation in the global coordinates";

    // register all necessary parameters (additional to those registered in GridTheta_k4geo)
    registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
    registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
    registerParameter("detLayout", "The detector layout (0 = Barrel; 1 = Endcap)", m_detLayout, -1);
    registerParameter("offset_z", "Offset in z-axis of the layer center", m_offsetZ, std::vector<double>());
    registerParameter("width_z", "Width in z of the layer", m_widthZ, std::vector<double>());
    registerParameter("offset_r", "Offset in radius of the layer (Rmin)", m_offsetR, std::vector<double>());
    registerParameter("numLayers", "Number of layers", m_numLayers, std::vector<int>());
    registerParameter("dRlayer", "dR of the layer", m_dRlayer, std::vector<double>());
    registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
    registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");
  }

  FCCSWHCalPhiTheta_k4geo::FCCSWHCalPhiTheta_k4geo(const BitFieldCoder* decoder) : GridTheta_k4geo(decoder) {
    // define type and description
    _type = "FCCSWHCalPhiTheta_k4geo";
    _description = "Phi-theta segmentation in the global coordinates";

    // register all necessary parameters (additional to those registered in GridTheta_k4geo)
    registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
    registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
    registerParameter("detLayout", "The detector layout (0 = Barrel; 1 = Endcap)", m_detLayout, -1);
    registerParameter("offset_z", "Offset in z-axis of the layer center", m_offsetZ, std::vector<double>());
    registerParameter("width_z", "Width in z of the layer", m_widthZ, std::vector<double>());
    registerParameter("offset_r", "Offset in radius of the layer (Rmin)", m_offsetR, std::vector<double>());
    registerParameter("numLayers", "Number of layers", m_numLayers, std::vector<int>());
    registerParameter("dRlayer", "dR of the layer", m_dRlayer, std::vector<double>());
    registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
    registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");
  }

  /** /// determine the global position based on the cell ID
  Vector3D FCCSWHCalPhiTheta_k4geo::position(const CellID& cID) const {
    uint layer = _decoder->get(cID,m_layerID);
    double radius = 1.0;

    if(m_radii.empty()) defineCellsInRZplan();
    if(!m_radii.empty()) radius = m_radii[layer];

    return positionFromRThetaPhi(radius, theta(cID), phi(cID));
  }
  **/

  /// determine the global position based on the cell ID
  /// returns the geometric center of the cell
  Vector3D FCCSWHCalPhiTheta_k4geo::position(const CellID& cID) const {
    uint layer = _decoder->get(cID, m_layerID);
    int thetaID = _decoder->get(cID, m_thetaID);
    double zpos = 0.;
    double radius = 1.0;

    if (m_radii.empty())
      defineCellsInRZplan();
    if (!m_radii.empty())
      radius = m_radii[layer];
    if (!m_cellEdges.empty())
      zpos = m_cellEdges[layer][thetaID].first +
             (m_cellEdges[layer][thetaID].second - m_cellEdges[layer][thetaID].first) * 0.5;

    auto pos = positionFromRThetaPhi(radius, theta(cID), phi(cID));

    // return the position with corrected z coordinate to match to the geometric center
    return Vector3D(pos.x(), pos.y(), zpos);
  }

  void FCCSWHCalPhiTheta_k4geo::defineCellsInRZplan() const {
    if (m_radii.empty()) {
      // check if all necessary variables are available
      if (m_detLayout == -1 || m_offsetZ.empty() || m_widthZ.empty() || m_offsetR.empty() || m_numLayers.empty() ||
          m_dRlayer.empty()) {
        dd4hep::printout(
            dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo", "Please check the readout description in the XML file!\n%s",
            "One of the variables is missing: detLayout | offset_z | width_z | offset_r | numLayers | dRlayer");
        return;
      }

      // some sanity checks of the xml
      if (m_offsetZ.size() != m_offsetR.size()) {
        dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo",
                         "Please check the readout description in the XML file!\n%s",
                         "Number of elements in offsetZ and offsetR must be the same!");
        return;
      }
      if (m_widthZ.size() != m_offsetR.size()) {
        dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo",
                         "Please check the readout description in the XML file!\n%s",
                         "Number of elements in widthZ and offsetR must be the same!");
        return;
      }
      if (m_detLayout == 0 && m_offsetZ.size() != 1) {
        dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo",
                         "Please check the readout description in the XML file!\n%s",
                         "Number of elements in offsetZ/offsetR/widthZ must be 1 for the Barrel!");
        return;
      }
      if (m_numLayers.size() % m_offsetZ.size() != 0) {
        dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo",
                         "Please check the readout description in the XML file!\n%s",
                         "Number of elements in numLayers must be multiple of offsetZ.size()!");
        return;
      }
      if (m_dRlayer.size() != m_numLayers.size() / m_offsetZ.size()) {
        dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo",
                         "Please check the readout description in the XML file!\n%s",
                         "Number of elements in dRlayer must be equal to numLayers.size()/offsetZ.size()!");
        return;
      }

      if (m_detLayout == 0)
        dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiTheta_k4geo", "Barrel configuration found!");
      else
        dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiTheta_k4geo", "EndCap configuration found!");

      // calculate the radius for each layer
      uint N_dR = m_numLayers.size() / m_offsetZ.size();
      std::vector<double> moduleDepth(m_offsetZ.size());
      for (uint i_section = 0; i_section < m_offsetZ.size(); i_section++) {
        for (uint i_dR = 0; i_dR < N_dR; i_dR++) {
          for (int i_row = 1; i_row <= m_numLayers[i_dR + i_section * N_dR]; i_row++) {
            moduleDepth[i_section] += m_dRlayer[i_dR];
            m_radii.push_back(m_offsetR[i_section] + moduleDepth[i_section] - m_dRlayer[i_dR] * 0.5);
            // layer lower and upper edges in z-axis
            m_layerEdges.push_back(std::make_pair(m_offsetZ[i_section] - 0.5 * m_widthZ[i_section],
                                                  m_offsetZ[i_section] + 0.5 * m_widthZ[i_section]));
            m_layerDepth.push_back(m_dRlayer[i_dR]);
          }
        }
      }

      // print info of calculated radii and edges
      for (uint i_layer = 0; i_layer < m_radii.size(); i_layer++) {
        dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiTheta_k4geo", "layer %d radius: %.2f, z range: %.2f - %.2f cm",
                         i_layer, m_radii[i_layer], m_layerEdges[i_layer].first, m_layerEdges[i_layer].second);
      }

      // allocate thetaBins vector for each layer
      m_thetaBins.resize(m_radii.size());
      // allocate cellEdges vector for each layer
      m_cellEdges.resize(m_radii.size());

      // determine theta bins and cell edges for each layer
      for (uint i_layer = 0; i_layer < m_radii.size(); i_layer++)
        defineCellEdges(i_layer);
    }
  }

  void FCCSWHCalPhiTheta_k4geo::defineCellEdges(const uint layer) const {
    if (m_thetaBins[layer].size() == 0 && m_radii.size() > 0) {
      // find theta bins that fit within the given layer
      int ibin =
          positionToBin(0.02, m_gridSizeTheta, m_offsetTheta); // <--- start from theta bin outside the HCal theta range
      while (m_radii[layer] * std::cos(m_offsetTheta + ibin * m_gridSizeTheta) /
                 std::sin(m_offsetTheta + ibin * m_gridSizeTheta) >
             m_layerEdges[layer].first) {
        if (m_radii[layer] * std::cos(m_offsetTheta + ibin * m_gridSizeTheta) /
                std::sin(m_offsetTheta + ibin * m_gridSizeTheta) <
            m_layerEdges[layer].second) {
          m_thetaBins[layer].push_back(ibin);
        }
        ibin++;
      }

      // find edges of each cell (theta bin) in the given layer
      auto prevBin = m_thetaBins[layer][0];
      // set the upper edge of the first cell in the given layer (starting from positive z part)
      m_cellEdges[layer][prevBin] = std::make_pair(0., m_layerEdges[layer].second);
      for (auto bin : m_thetaBins[layer]) {
        if (bin != prevBin) {
          double z1 = m_radii[layer] * std::cos(m_offsetTheta + bin * m_gridSizeTheta) /
                      std::sin(m_offsetTheta + bin * m_gridSizeTheta);
          double z2 = m_radii[layer] * std::cos(m_offsetTheta + prevBin * m_gridSizeTheta) /
                      std::sin(m_offsetTheta + prevBin * m_gridSizeTheta);
          // set the lower edge of the prevBin cell
          m_cellEdges[layer][prevBin].first = z1 + 0.5 * (z2 - z1);
          // set the upper edge of current bin cell
          m_cellEdges[layer][bin] = std::make_pair(0., m_cellEdges[layer][prevBin].first);
          prevBin = bin;
        }
      }
      // set the lower edge of the last cell in the given layer
      m_cellEdges[layer][prevBin].first = m_layerEdges[layer].first;

      // for the EndCap, do it again but for negative z part
      if (m_detLayout == 1) {
        while (m_radii[layer] * std::cos(m_offsetTheta + ibin * m_gridSizeTheta) /
                   std::sin(m_offsetTheta + ibin * m_gridSizeTheta) >
               (-m_layerEdges[layer].second)) {
          if (m_radii[layer] * std::cos(m_offsetTheta + ibin * m_gridSizeTheta) /
                  std::sin(m_offsetTheta + ibin * m_gridSizeTheta) <
              (-m_layerEdges[layer].first)) {
            m_thetaBins[layer].push_back(ibin);
          }
          ibin++;
        }

        // Create a span view over the theta bins corresponding to the Endcap in negative z part
        std::span<int> thetaBins(m_thetaBins[layer].begin() + m_thetaBins[layer].size() / 2, m_thetaBins[layer].end());
        prevBin = thetaBins[0];

        // set the upper edge of the first cell in the given layer at negative z part
        m_cellEdges[layer][prevBin] = std::make_pair(0., -m_layerEdges[layer].first);
        for (auto bin : thetaBins) {
          if (bin != prevBin) {
            double z1 = m_radii[layer] * std::cos(m_offsetTheta + bin * m_gridSizeTheta) /
                        std::sin(m_offsetTheta + bin * m_gridSizeTheta);
            double z2 = m_radii[layer] * std::cos(m_offsetTheta + prevBin * m_gridSizeTheta) /
                        std::sin(m_offsetTheta + prevBin * m_gridSizeTheta);
            // set the lower edge of the prevBin cell
            m_cellEdges[layer][prevBin].first = z1 + 0.5 * (z2 - z1);
            // set the upper edge of current bin cell
            m_cellEdges[layer][bin] = std::make_pair(0., m_cellEdges[layer][prevBin].first);
            prevBin = bin;
          }
        }
        // set the lower edge of the last cell in the given layer
        m_cellEdges[layer][prevBin].first = (-m_layerEdges[layer].second);
      } // negative-z endcap

      dd4hep::printout(dd4hep::DEBUG, "FCCSWHCalPhiTheta_k4geo", "Number of cells in layer %d: %d", layer,
                       m_thetaBins[layer].size());
      for (auto bin : m_thetaBins[layer])
        dd4hep::printout(dd4hep::DEBUG, "FCCSWHCalPhiTheta_k4geo", "Layer %d cell theta bin: %d, edges: %.2f - %.2f cm",
                         layer, bin, m_cellEdges[layer][bin]);
    }
  }

  /// create the cell ID based on the position
  CellID FCCSWHCalPhiTheta_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                                         const VolumeID& vID) const {

    CellID cID = vID;

    // The volume ID comes with "row" field information (number of sequences) that would screw up the topo-clustering
    // using cell neighbours map produced with RecFCCeeCalorimeter/src/components/CreateFCCeeCaloNeighbours.cpp,
    // therefore, lets set it to zero, as it is for the cell IDs in the neighbours map.
    _decoder->set(cID, "row", 0);

    // For endcap, the volume ID comes with "type" field information which would screw up the topo-clustering as the
    // "row" field, therefore, lets set it to zero, as it is for the cell IDs in the neighbours map.
    if (m_detLayout == 1)
      _decoder->set(cID, "type", 0);

    double lTheta = thetaFromXYZ(globalPosition);
    double lPhi = phiFromXYZ(globalPosition);
    uint layer = _decoder->get(vID, m_layerID);

    // define cell boundaries in R-z plan
    if (m_radii.empty())
      defineCellsInRZplan();

    // check if the cells are defined for the given layer
    if (m_thetaBins[layer].empty())
      dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo", "No cells are defined for layer %d", layer);

    // find the cell (theta bin) corresponding to the hit and return the cellID
    for (auto bin : m_thetaBins[layer]) {
      double posz = globalPosition.z();
      if (posz > m_cellEdges[layer][bin].first && posz < m_cellEdges[layer][bin].second) {
        _decoder->set(cID, m_thetaID, bin);
        _decoder->set(cID, m_phiID, positionToBin(lPhi, 2 * M_PI / (double)m_phiBins, m_offsetPhi));
        return cID;
      }
    }

    dd4hep::printout(dd4hep::WARNING, "FCCSWHCalPhiTheta_k4geo", "The hit is outside the defined range of the layer %d",
                     layer);

    _decoder->set(cID, m_thetaID, positionToBin(lTheta, m_gridSizeTheta, m_offsetTheta));
    _decoder->set(cID, m_phiID, positionToBin(lPhi, 2 * M_PI / (double)m_phiBins, m_offsetPhi));
    return cID;
  }

  /// determine the azimuthal angle phi based on the cell ID
  double FCCSWHCalPhiTheta_k4geo::phi(const CellID& cID) const {
    CellID phiValue = _decoder->get(cID, m_phiID);
    return binToPosition(phiValue, 2. * M_PI / (double)m_phiBins, m_offsetPhi);
  }

  /// Get the min and max layer indexes for each part of the HCal
  std::vector<std::pair<uint, uint>> FCCSWHCalPhiTheta_k4geo::getMinMaxLayerId() const {
    std::vector<std::pair<uint, uint>> minMaxLayerId;

    if (m_radii.empty())
      defineCellsInRZplan();
    if (m_radii.empty())
      return minMaxLayerId;

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
  std::vector<uint64_t> FCCSWHCalPhiTheta_k4geo::neighbours(const CellID& cID, bool aDiagonal) const {
    std::vector<uint64_t> cellNeighbours;

    if (m_radii.empty())
      defineCellsInRZplan();
    if (m_thetaBins.empty())
      return cellNeighbours;

    uint EndcapPart = 0;
    int minLayerId = -1;
    int maxLayerId = -1;

    int currentLayerId = _decoder->get(cID, m_layerID);
    int currentCellThetaBin = _decoder->get(cID, m_thetaID);

    int minCellThetaBin = m_thetaBins[currentLayerId].front();
    int maxCellThetaBin = m_thetaBins[currentLayerId].back();

    //--------------------------------
    // Determine min and max layer Id
    //--------------------------------
    std::vector<int> minLayerIdEndcap(m_offsetZ.size(), 0);
    std::vector<int> maxLayerIdEndcap(m_offsetZ.size(), 0);
    if (m_detLayout == 1) {
      std::vector<std::pair<uint, uint>> minMaxLayerId(getMinMaxLayerId());
      if (minMaxLayerId.empty()) {
        dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo",
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

      // correct the min and max theta bin for endcap
      if (theta(cID) > M_PI / 2) // negative-z part
      {
        // second half of elements in m_thetaBins[currentLayerId] vector corresponds to the negative-z layer cells
        minCellThetaBin = m_thetaBins[currentLayerId][m_thetaBins[currentLayerId].size() / 2];
        maxCellThetaBin = m_thetaBins[currentLayerId].back();
      } else // positive-z part
      {
        // first half of elements in m_thetaBins[currentLayerId] vector corresponds to the positive-z layer cells
        minCellThetaBin = m_thetaBins[currentLayerId].front();
        maxCellThetaBin = m_thetaBins[currentLayerId][m_thetaBins[currentLayerId].size() / 2 - 1];
      }
    } else // for Barrel
    {
      minLayerId = 0;
      maxLayerId = m_radii.size() - 1;
    }
    //--------------------------------

    //--------------------------------
    // Find neighbours
    //--------------------------------

    //---------------------------------------------
    // this part is same for both Barrel and Endcap
    //---------------------------------------------
    // if this is not the first cell in the given layer then add the previous cell
    if (currentCellThetaBin > minCellThetaBin) {
      CellID nID = cID;
      _decoder->set(nID, m_thetaID, currentCellThetaBin - 1);
      cellNeighbours.push_back(nID); // add the previous cell from current layer of the same phi module
    }
    // if this is not the last cell in the given layer then add the next cell
    if (currentCellThetaBin < maxCellThetaBin) {
      CellID nID = cID;
      _decoder->set(nID, m_thetaID, currentCellThetaBin + 1);
      cellNeighbours.push_back(nID); // add the next cell from current layer of the same phi module
    }
    //----------------------------------------------

    // deal with the Barrel
    if (m_detLayout == 0) {
      double currentCellZmin = m_cellEdges[currentLayerId][currentCellThetaBin].first;
      double currentCellZmax = m_cellEdges[currentLayerId][currentCellThetaBin].second;

      // if this is not the first layer then look for neighbours in the previous layer
      if (currentLayerId > minLayerId) {
        CellID nID = cID;
        int prevLayerId = currentLayerId - 1;
        _decoder->set(nID, m_layerID, prevLayerId);

        _decoder->set(nID, m_thetaID, currentCellThetaBin);
        cellNeighbours.push_back(
            nID); // add the cell with the same theta bin from the previous layer of the same phi module

        // if the cID is in the positive-z side and prev layer cell is not in the first theta bin then add the cell from
        // previous theta bin
        if (theta(cID) < M_PI / 2. && currentCellThetaBin > m_thetaBins[prevLayerId].front()) {
          _decoder->set(nID, m_thetaID, currentCellThetaBin - 1);
          cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
          if (aDiagonal && currentCellThetaBin > (m_thetaBins[prevLayerId].front() + 1)) {
            // add the previous layer cell from the prev to prev theta bin if it overlaps with the current cell in
            // z-coordinate
            double zmin = m_cellEdges[prevLayerId][currentCellThetaBin - 2].first;
            if (zmin <= currentCellZmax) {
              // add the previous layer cell from the prev to prev theta bin
              _decoder->set(nID, m_thetaID, currentCellThetaBin - 2);
              cellNeighbours.push_back(nID);
            }
          }
        }
        // if the cID is in the negative-z side and prev layer cell is not in the last theta bin then add the cell from
        // previous theta bin
        if (theta(cID) > M_PI / 2. && currentCellThetaBin < m_thetaBins[prevLayerId].back()) {
          _decoder->set(nID, m_thetaID, currentCellThetaBin + 1);
          cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
          if (aDiagonal && currentCellThetaBin < (m_thetaBins[prevLayerId].back() - 1)) {
            // add the previous layer cell from the next to next theta bin if it overlaps with the current cell in
            // z-coordinate
            double zmax = m_cellEdges[prevLayerId][currentCellThetaBin + 2].second;
            if (zmax >= currentCellZmin) {
              // add the previous layer cell from the next to next theta bin
              _decoder->set(nID, m_thetaID, currentCellThetaBin + 2);
              cellNeighbours.push_back(nID);
            }
          }
        }
      }

      // if this is not the last layer then look for neighbours in the next layer
      if (currentLayerId < maxLayerId) {
        CellID nID = cID;
        int nextLayerId = currentLayerId + 1;
        _decoder->set(nID, m_layerID, nextLayerId);

        _decoder->set(nID, m_thetaID, currentCellThetaBin);
        cellNeighbours.push_back(
            nID); // add the cell with the same theta bin from the next layer of the same phi module

        // if the cID is in the positive-z side
        if (theta(cID) < M_PI / 2.) {
          // add the next layer cell from the next theta bin
          _decoder->set(nID, m_thetaID, currentCellThetaBin + 1);
          cellNeighbours.push_back(nID);

          if (aDiagonal) {
            // add the next layer cell from the next-to-next theta bin if it overlaps with the current cell in
            // z-coordinate
            double zmax = m_cellEdges[nextLayerId][currentCellThetaBin + 2].second;
            if (zmax >= currentCellZmin) {
              // add the next layer cell from the next to next theta bin
              _decoder->set(nID, m_thetaID, currentCellThetaBin + 2);
              cellNeighbours.push_back(nID);
            }
          }
        }
        // if the cID is in the negative-z side
        if (theta(cID) > M_PI / 2.) {
          // add the next layer cell from the previous theta bin
          _decoder->set(nID, m_thetaID, currentCellThetaBin - 1);
          cellNeighbours.push_back(nID);

          if (aDiagonal) {
            // add the next layer cell from the prev to prev theta bin if it overlaps with the current cell in
            // z-coordinate
            double zmin = m_cellEdges[nextLayerId][currentCellThetaBin - 2].first;
            if (zmin <= currentCellZmax) {
              // add the next layer cell from the prev to prev theta bin
              _decoder->set(nID, m_thetaID, currentCellThetaBin - 2);
              cellNeighbours.push_back(nID);
            }
          }
        }
      }
    }

    // Endcap
    if (m_detLayout == 1) {
      double currentCellZmin = m_cellEdges[currentLayerId][currentCellThetaBin].first;
      double currentCellZmax = m_cellEdges[currentLayerId][currentCellThetaBin].second;

      // if this is not the first layer then look for neighbours in the previous layer
      if (currentLayerId > minLayerId) {
        CellID nID = cID;
        int prevLayerId = currentLayerId - 1;
        _decoder->set(nID, m_layerID, prevLayerId);
        // find the ones that share at least part of a border with the current cell
        for (auto bin : m_thetaBins[prevLayerId]) {
          double zmin = m_cellEdges[prevLayerId][bin].first;
          double zmax = m_cellEdges[prevLayerId][bin].second;

          // if the cID is in the positive-z side
          if (theta(cID) < M_PI / 2.) {
            if ((zmin >= currentCellZmin && zmin < currentCellZmax) ||
                (zmax >= currentCellZmin && zmax <= currentCellZmax) ||
                (currentCellZmin >= zmin && currentCellZmax <= zmax)) {
              _decoder->set(nID, m_thetaID, bin);
              cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
            }
            if (aDiagonal && zmin == currentCellZmax) {
              _decoder->set(nID, m_thetaID, bin);
              cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
            }
          }
          // if the cID is in the negative-z side
          if (theta(cID) > M_PI / 2.) {
            if ((zmin >= currentCellZmin && zmin <= currentCellZmax) ||
                (zmax > currentCellZmin && zmax <= currentCellZmax) ||
                (currentCellZmin >= zmin && currentCellZmax <= zmax)) {
              _decoder->set(nID, m_thetaID, bin);
              cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
            }
            if (aDiagonal && zmax == currentCellZmin) {
              _decoder->set(nID, m_thetaID, bin);
              cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
            }
          }
        }
      }
      // if this is not the last layer then look for neighbours in the next layer
      if (currentLayerId < maxLayerId) {
        CellID nID = cID;
        int nextLayerId = currentLayerId + 1;
        _decoder->set(nID, m_layerID, nextLayerId);
        // find the ones that share at least part of a border with the current cell
        for (auto bin : m_thetaBins[nextLayerId]) {
          double zmin = m_cellEdges[nextLayerId][bin].first;
          double zmax = m_cellEdges[nextLayerId][bin].second;
          // if the cID is in the positive-z side
          if (theta(cID) < M_PI / 2.) {
            if ((zmin >= currentCellZmin && zmin <= currentCellZmax) ||
                (zmax > currentCellZmin && zmax <= currentCellZmax) ||
                (currentCellZmin >= zmin && currentCellZmax <= zmax)) {
              _decoder->set(nID, m_thetaID, bin);
              cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
            }
            if (aDiagonal && zmax == currentCellZmin) {
              _decoder->set(nID, m_thetaID, bin);
              cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
            }
          }
          // if the cID is in the negative-z side
          if (theta(cID) > M_PI / 2.) {
            if ((zmin >= currentCellZmin && zmin < currentCellZmax) ||
                (zmax >= currentCellZmin && zmax <= currentCellZmax) ||
                (currentCellZmin >= zmin && currentCellZmax <= zmax)) {
              _decoder->set(nID, m_thetaID, bin);
              cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
            }
            if (aDiagonal && zmin == currentCellZmax) {
              _decoder->set(nID, m_thetaID, bin);
              cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
            }
          }
        }
      }

      // if the Endcap consists of more than 1 part/section then look for neighbours in different parts as well
      if (m_offsetZ.size() > 1) {
        //
        double currentLayerRmin = m_radii[currentLayerId] - 0.5 * m_layerDepth[currentLayerId];
        double currentLayerRmax = m_radii[currentLayerId] + 0.5 * m_layerDepth[currentLayerId];

        // if the cell is in negative-z part, then swap min and max theta bins
        if (theta(cID) > M_PI / 2.) {
          minCellThetaBin = m_thetaBins[currentLayerId].back();
          maxCellThetaBin = m_thetaBins[currentLayerId][m_thetaBins[currentLayerId].size() / 2];
        }

        // if it is the last cell in the part1
        if (EndcapPart == 0 && currentCellThetaBin == minCellThetaBin) {
          // find the layers in the part2 that share a border with the current layer
          for (int part2layerId = minLayerIdEndcap[1]; part2layerId <= maxLayerIdEndcap[1]; part2layerId++) {
            double Rmin = m_radii[part2layerId] - 0.5 * m_layerDepth[part2layerId];
            double Rmax = m_radii[part2layerId] + 0.5 * m_layerDepth[part2layerId];

            if ((Rmin >= currentLayerRmin && Rmin <= currentLayerRmax) ||
                (Rmax > currentLayerRmin && Rmax <= currentLayerRmax) ||
                (currentLayerRmin >= Rmin && currentLayerRmin < Rmax) ||
                (currentLayerRmax >= Rmin && currentLayerRmax <= Rmax)) {
              CellID nID = cID;
              _decoder->set(nID, m_layerID, part2layerId);
              _decoder->set(nID, m_thetaID, maxCellThetaBin);
              cellNeighbours.push_back(nID); // add the last theta bin cell from part2 layer
            }
            if (aDiagonal && Rmax == currentLayerRmin) {
              CellID nID = cID;
              _decoder->set(nID, m_layerID, part2layerId);
              _decoder->set(nID, m_thetaID, maxCellThetaBin);
              cellNeighbours.push_back(nID); // add the last theta bin cell from part2 layer
            }
          }
        }

        // if the Endcap consists of more than 2 parts:
        for (uint i_section = 1; i_section < (m_offsetZ.size() - 1); i_section++) {
          if (i_section != EndcapPart)
            continue;

          // if it is the last theta bin cell then look for neighbours in previous part
          if (currentCellThetaBin == maxCellThetaBin) {
            // find the layers in the previous part that share a border with the current layer
            for (int prevPartLayerId = minLayerIdEndcap[i_section - 1];
                 prevPartLayerId <= maxLayerIdEndcap[i_section - 1]; prevPartLayerId++) {
              double Rmin = m_radii[prevPartLayerId] - 0.5 * m_layerDepth[prevPartLayerId];
              double Rmax = m_radii[prevPartLayerId] + 0.5 * m_layerDepth[prevPartLayerId];

              if ((Rmin >= currentLayerRmin && Rmin < currentLayerRmax) ||
                  (Rmax >= currentLayerRmin && Rmax <= currentLayerRmax) ||
                  (currentLayerRmin >= Rmin && currentLayerRmin <= Rmax) ||
                  (currentLayerRmax > Rmin && currentLayerRmax <= Rmax)) {
                CellID nID = cID;
                _decoder->set(nID, m_layerID, prevPartLayerId);
                _decoder->set(nID, m_thetaID, minCellThetaBin);
                cellNeighbours.push_back(nID); // add the first theta bin cell from the part1 layer
              }
              if (aDiagonal && Rmin == currentLayerRmax) {
                CellID nID = cID;
                _decoder->set(nID, m_layerID, prevPartLayerId);
                _decoder->set(nID, m_thetaID, minCellThetaBin);
                cellNeighbours.push_back(nID); // add the first theta bin cell from the part1 layer
              }
            }
          }

          // if it is the first theta bin cell then look for neighbours in the next part
          if (currentCellThetaBin == minCellThetaBin) {
            // find the layers in the next part that share a border with the current layer
            for (int nextPartLayerId = minLayerIdEndcap[i_section + 1];
                 nextPartLayerId <= maxLayerIdEndcap[i_section + 1]; nextPartLayerId++) {
              double Rmin = m_radii[nextPartLayerId] - 0.5 * m_layerDepth[nextPartLayerId];
              double Rmax = m_radii[nextPartLayerId] + 0.5 * m_layerDepth[nextPartLayerId];

              if ((Rmin >= currentLayerRmin && Rmin <= currentLayerRmax) ||
                  (Rmax > currentLayerRmin && Rmax <= currentLayerRmax) ||
                  (currentLayerRmin >= Rmin && currentLayerRmin < Rmax) ||
                  (currentLayerRmax >= Rmin && currentLayerRmax <= Rmax)) {
                CellID nID = cID;
                _decoder->set(nID, m_layerID, nextPartLayerId);
                _decoder->set(nID, m_thetaID, maxCellThetaBin);
                cellNeighbours.push_back(nID); // add the first cell from the part3 layer
              }
              if (aDiagonal && Rmax == currentLayerRmin) {
                CellID nID = cID;
                _decoder->set(nID, m_layerID, nextPartLayerId);
                _decoder->set(nID, m_thetaID, maxCellThetaBin);
                cellNeighbours.push_back(nID); // add the first cell from the part3 layer
              }
            }
          }
        }

        // if it is the last theta bin cell in the last part of the Endcap
        if (EndcapPart == (m_offsetZ.size() - 1) && currentCellThetaBin == maxCellThetaBin) {
          // find the layers in the previous part that share a border with the current layer
          for (int prevPartLayerId = minLayerIdEndcap[m_offsetZ.size() - 2];
               prevPartLayerId <= maxLayerIdEndcap[m_offsetZ.size() - 2]; prevPartLayerId++) {
            double Rmin = m_radii[prevPartLayerId] - 0.5 * m_layerDepth[prevPartLayerId];
            double Rmax = m_radii[prevPartLayerId] + 0.5 * m_layerDepth[prevPartLayerId];

            if ((Rmin >= currentLayerRmin && Rmin < currentLayerRmax) ||
                (Rmax >= currentLayerRmin && Rmax <= currentLayerRmax) ||
                (currentLayerRmin >= Rmin && currentLayerRmin <= Rmax) ||
                (currentLayerRmax > Rmin && currentLayerRmax <= Rmax)) {
              CellID nID = cID;
              _decoder->set(nID, m_layerID, prevPartLayerId);
              _decoder->set(nID, m_thetaID, minCellThetaBin);
              cellNeighbours.push_back(nID); // add the first theta bin cell from the part2 layer
            }
            if (aDiagonal && Rmin == currentLayerRmax) {
              CellID nID = cID;
              _decoder->set(nID, m_layerID, prevPartLayerId);
              _decoder->set(nID, m_thetaID, minCellThetaBin);
              cellNeighbours.push_back(nID); // add the first theta bin cell from the part2 layer
            }
          }
        }
      }
    }

    // Now loop over the neighbours and add the cells from next/previous phi module
    std::vector<uint64_t> cellNeighboursCopy(cellNeighbours);
    for (auto nID : cellNeighboursCopy) {
      CellID newID = nID;
      // previous: if the current is 0 then previous is the last bin (id = m_phiBins - 1) else current - 1
      _decoder->set(newID, m_phiID,
                    (_decoder->get(nID, m_phiID) == 0) ? m_phiBins - 1 : _decoder->get(nID, m_phiID) - 1);
      cellNeighbours.push_back(newID);
      // next: if the current is the last bin (id = m_phiBins - 1) then the next is the first bin (id = 0) else current
      // + 1
      _decoder->set(newID, m_phiID,
                    (_decoder->get(nID, m_phiID) == (m_phiBins - 1)) ? 0 : _decoder->get(nID, m_phiID) + 1);
      cellNeighbours.push_back(newID);
    }

    // At the end, find neighbours with the same layer/row in next/previous phi module
    CellID nID = cID;
    // previous: if the current is 0 then previous is the last bin (id = m_phiBins - 1) else current - 1
    _decoder->set(nID, m_phiID, (_decoder->get(cID, m_phiID) == 0) ? m_phiBins - 1 : _decoder->get(cID, m_phiID) - 1);
    cellNeighbours.push_back(nID);
    // next: if the current is the last bin (id = m_phiBins - 1) then the next is the first bin (id = 0) else current +
    // 1
    _decoder->set(nID, m_phiID, (_decoder->get(cID, m_phiID) == (m_phiBins - 1)) ? 0 : _decoder->get(cID, m_phiID) + 1);
    cellNeighbours.push_back(nID);

    return cellNeighbours;
  }

  /// Determine minimum and maximum polar angle of the cell
  std::array<double, 2> FCCSWHCalPhiTheta_k4geo::cellTheta(const CellID& cID) const {
    std::array<double, 2> cTheta = {M_PI, M_PI};

    // get the cell index
    int idx = _decoder->get(cID, m_thetaID);
    // get the layer index
    uint layer = _decoder->get(cID, m_layerID);

    if (m_radii.empty())
      defineCellsInRZplan();
    if (m_cellEdges.empty())
      return cTheta;

    double zlow = m_cellEdges[layer][idx].first;
    double zhigh = m_cellEdges[layer][idx].second;

    double Rmin = m_radii[layer] - 0.5 * m_layerDepth[layer];
    double Rmax = m_radii[layer] + 0.5 * m_layerDepth[layer];

    if (theta(cID) < M_PI / 2.) {
      cTheta[0] = std::atan2(Rmin, zhigh); // theta min
      cTheta[1] = std::atan2(Rmax, zlow);  // theta max
    } else {
      cTheta[0] = std::atan2(Rmax, zhigh); // theta min
      cTheta[1] = std::atan2(Rmin, zlow);  // theta max
    }

    return cTheta;
  }

} // namespace DDSegmentation
} // namespace dd4hep
