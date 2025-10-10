#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"
#include "DD4hep/Printout.h"

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  FCCSWHCalPhiTheta_k4geo::FCCSWHCalPhiTheta_k4geo(const std::string& cellEncoding) : GridTheta_k4geo(cellEncoding) {
    commonSetup();
  }

  FCCSWHCalPhiTheta_k4geo::FCCSWHCalPhiTheta_k4geo(const BitFieldCoder* decoder) : GridTheta_k4geo(decoder) {
    commonSetup();
  }

  /// Initialization common to all ctors.
  void FCCSWHCalPhiTheta_k4geo::commonSetup() {
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

    m_layerIndex = decoder()->index(m_layerID);
    m_rowIndex = decoder()->index("row");
    m_thetaIndex = decoder()->index(fieldNameTheta());
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
  /// returns the geometric center of the cell
  Vector3D FCCSWHCalPhiTheta_k4geo::position(const CellID& cID) const {
    uint layer = decoder()->get(cID, m_layerIndex);
    int thetaID = decoder()->get(cID, m_thetaIndex);
    const LayerInfo& li = getLayerInfo(layer);
    double radius = li.radius;
    const auto& edges = li.cellEdges.find(thetaID)->second;
    double zpos = (edges.first + edges.second) * 0.5;

    auto pos = positionFromRThetaPhi(radius, theta(cID), phi(cID));

    // return the position with corrected z coordinate to match to the geometric center
    return Vector3D(pos.x(), pos.y(), zpos);
  }

  auto FCCSWHCalPhiTheta_k4geo::getLayerInfo(const unsigned layer) const -> const LayerInfo& {
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
  std::vector<FCCSWHCalPhiTheta_k4geo::LayerInfo> FCCSWHCalPhiTheta_k4geo::initLayerInfo() const {
    std::vector<LayerInfo> out;
    if (!checkParameters()) {
      out.resize((*decoder())[m_layerIndex].maxValue() + 1);
      return out;
    }

    if (m_detLayout == 0)
      dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiTheta_k4geo", "Barrel configuration found!");
    else
      dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiTheta_k4geo", "EndCap configuration found!");

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
      dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiTheta_k4geo", "layer %d radius: %.2f, z range: %.2f - %.2f cm",
                       i_layer++, li.radius, li.zmin, li.zmax);
    }

    // determine theta bins and cell edges for each layer
    for (uint i_layer = 0; LayerInfo& li : out) {
      defineCellEdges(li, i_layer);
      ++i_layer;
    }
    return out;
  }

  void FCCSWHCalPhiTheta_k4geo::defineCellEdges(LayerInfo& li, const unsigned int layer) const {
    // find theta bins that fit within the given layer
    int ibin =
        positionToBin(0.02, gridSizeTheta(), offsetTheta()); // <--- start from theta bin outside the HCal theta range
    while (li.radius * std::cos(offsetTheta() + ibin * gridSizeTheta()) /
               std::sin(offsetTheta() + ibin * gridSizeTheta()) >
           li.zmin) {
      if (li.radius * std::cos(offsetTheta() + ibin * gridSizeTheta()) /
              std::sin(offsetTheta() + ibin * gridSizeTheta()) <
          li.zmax) {
        li.thetaBins.push_back(ibin);
      }
      ibin++;
    }

    // find edges of each cell (theta bin) in the given layer
    auto prevBin = li.thetaBins[0];
    // set the upper edge of the first cell in the given layer (starting from positive z part)
    li.cellEdges[prevBin] = std::make_pair(0., li.zmax);
    for (auto bin : li.thetaBins) {
      if (bin != prevBin) {
        double z1 = li.radius * std::cos(offsetTheta() + bin * gridSizeTheta()) /
                    std::sin(offsetTheta() + bin * gridSizeTheta());
        double z2 = li.radius * std::cos(offsetTheta() + prevBin * gridSizeTheta()) /
                    std::sin(offsetTheta() + prevBin * gridSizeTheta());
        // set the lower edge of the prevBin cell
        li.cellEdges[prevBin].first = z1 + 0.5 * (z2 - z1);
        // set the upper edge of current bin cell
        li.cellEdges[bin] = std::make_pair(0., li.cellEdges[prevBin].first);
        prevBin = bin;
      }
    }

    // set the lower edge of the last cell in the given layer
    li.cellEdges[prevBin].first = li.zmin;

    // for the EndCap, do it again but for negative z part
    if (m_detLayout == 1) {
      while (li.radius * std::cos(offsetTheta() + ibin * gridSizeTheta()) /
                 std::sin(offsetTheta() + ibin * gridSizeTheta()) >
             (-li.zmax)) {
        if (li.radius * std::cos(offsetTheta() + ibin * gridSizeTheta()) /
                std::sin(offsetTheta() + ibin * gridSizeTheta()) <
            (-li.zmin)) {
          li.thetaBins.push_back(ibin);
        }
        ibin++;
      }

      // Create a span view over the theta bins corresponding to the Endcap in negative z part
      std::span<int> thetaBins(li.thetaBins.begin() + li.thetaBins.size() / 2, li.thetaBins.end());
      prevBin = thetaBins[0];

      // set the upper edge of the first cell in the given layer at negative z part
      li.cellEdges[prevBin] = std::make_pair(0., -li.zmin);
      for (auto bin : thetaBins) {
        if (bin != prevBin) {
          double z1 = li.radius * std::cos(offsetTheta() + bin * gridSizeTheta()) /
                      std::sin(offsetTheta() + bin * gridSizeTheta());
          double z2 = li.radius * std::cos(offsetTheta() + prevBin * gridSizeTheta()) /
                      std::sin(offsetTheta() + prevBin * gridSizeTheta());
          // set the lower edge of the prevBin cell
          li.cellEdges[prevBin].first = z1 + 0.5 * (z2 - z1);
          // set the upper edge of current bin cell
          li.cellEdges[bin] = std::make_pair(0., li.cellEdges[prevBin].first);
          prevBin = bin;
        }
      }
      // set the lower edge of the last cell in the given layer
      li.cellEdges[prevBin].first = (-li.zmax);
    } // negative-z endcap

    dd4hep::printout(dd4hep::DEBUG, "FCCSWHCalPhiTheta_k4geo", "Number of cells in layer %d: %d", layer,
                     li.thetaBins.size());
    for (auto bin : li.thetaBins)
      dd4hep::printout(dd4hep::DEBUG, "FCCSWHCalPhiTheta_k4geo", "Layer %d cell theta bin: %d, edges: %.2f - %.2f cm",
                       layer, bin, li.cellEdges[bin].first, li.cellEdges[bin].second);
  }

  // Check consistency of input geometric variables.
  bool FCCSWHCalPhiTheta_k4geo::checkParameters() const {
    // check if all necessary variables are available
    if (m_detLayout == -1 || m_offsetZ.empty() || m_widthZ.empty() || m_offsetR.empty() || m_numLayers.empty() ||
        m_dRlayer.empty()) {
      dd4hep::printout(
          dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo", "Please check the readout description in the XML file!\n%s",
          "One of the variables is missing: detLayout | offset_z | width_z | offset_r | numLayers | dRlayer");
      return false;
    }

    // some sanity checks of the xml
    if (m_offsetZ.size() != m_offsetR.size()) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo",
                       "Please check the readout description in the XML file!\n%s",
                       "Number of elements in offsetZ and offsetR must be the same!");
      return false;
    }

    if (m_widthZ.size() != m_offsetR.size()) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo",
                       "Please check the readout description in the XML file!\n%s",
                       "Number of elements in widthZ and offsetR must be the same!");
      return false;
    }

    if (m_detLayout == 0 && m_offsetZ.size() != 1) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo",
                       "Please check the readout description in the XML file!\n%s",
                       "Number of elements in offsetZ/offsetR/widthZ must be 1 for the Barrel!");
      return false;
    }

    if (m_numLayers.size() % m_offsetZ.size() != 0) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo",
                       "Please check the readout description in the XML file!\n%s",
                       "Number of elements in numLayers must be multiple of offsetZ.size()!");
      return false;
    }

    if (m_dRlayer.size() != m_numLayers.size() / m_offsetZ.size()) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo",
                       "Please check the readout description in the XML file!\n%s",
                       "Number of elements in dRlayer must be equal to numLayers.size()/offsetZ.size()!");
      return false;
    }

    return true;
  }

  /// create the cell ID based on the position
  CellID FCCSWHCalPhiTheta_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                                         const VolumeID& vID) const {

    CellID cID = vID;

    // The volume ID comes with "row" field information (number of sequences) that would screw up the topo-clustering
    // using cell neighbours map produced with RecFCCeeCalorimeter/src/components/CreateFCCeeCaloNeighbours.cpp,
    // therefore, lets set it to zero, as it is for the cell IDs in the neighbours map.
    decoder()->set(cID, m_rowIndex, 0);

    // For endcap, the volume ID comes with "type" field information which would screw up the topo-clustering as the
    // "row" field, therefore, lets set it to zero, as it is for the cell IDs in the neighbours map.
    if (m_detLayout == 1)
      decoder()->set(cID, m_typeIndex, 0);

    double lTheta = thetaFromXYZ(globalPosition);
    double lPhi = phiFromXYZ(globalPosition);
    uint layer = decoder()->get(vID, m_layerIndex);
    const LayerInfo& li = getLayerInfo(layer);

    // find the cell (theta bin) corresponding to the hit and return the cellID
    for (auto bin : li.thetaBins) {
      double posz = globalPosition.z();
      const auto& edges = li.cellEdges.find(bin)->second;
      if (posz > edges.first && posz < edges.second) {
        decoder()->set(cID, m_thetaIndex, bin);
        decoder()->set(cID, m_phiIndex, positionToBin(lPhi, 2 * M_PI / (double)m_phiBins, m_offsetPhi));
        return cID;
      }
    }

    dd4hep::printout(dd4hep::WARNING, "FCCSWHCalPhiTheta_k4geo", "The hit is outside the defined range of the layer %d",
                     layer);

    decoder()->set(cID, m_thetaIndex, positionToBin(lTheta, gridSizeTheta(), offsetTheta()));
    decoder()->set(cID, m_phiIndex, positionToBin(lPhi, 2 * M_PI / (double)m_phiBins, m_offsetPhi));
    return cID;
  }

  /// determine the azimuthal angle phi based on the cell ID
  double FCCSWHCalPhiTheta_k4geo::phi(const CellID cID) const {
    CellID phiValue = decoder()->get(cID, m_phiIndex);
    return binToPosition(phiValue, 2. * M_PI / (double)m_phiBins, m_offsetPhi);
  }

  /// Get the min and max layer indexes for each part of the HCal
  std::vector<std::pair<uint, uint>> FCCSWHCalPhiTheta_k4geo::getMinMaxLayerId() const {
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
  std::vector<uint64_t> FCCSWHCalPhiTheta_k4geo::neighbours(const CellID cID, bool aDiagonal) const {
    std::vector<uint64_t> cellNeighbours;

    uint EndcapPart = 0;
    int minLayerId = -1;
    int maxLayerId = -1;

    int currentLayerId = decoder()->get(cID, m_layerIndex);
    int currentCellThetaBin = decoder()->get(cID, m_thetaIndex);

    const LayerInfo& li = getLayerInfo(currentLayerId);

    int minCellThetaBin = li.thetaBins.front();
    int maxCellThetaBin = li.thetaBins.back();

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
        // second half of elements in li.thetaBins vector corresponds to the negative-z layer cells
        minCellThetaBin = li.thetaBins[li.thetaBins.size() / 2];
        maxCellThetaBin = li.thetaBins.back();
      } else // positive-z part
      {
        // first half of elements in li.thetaBins vector corresponds to the positive-z layer cells
        minCellThetaBin = li.thetaBins.front();
        maxCellThetaBin = li.thetaBins[li.thetaBins.size() / 2 - 1];
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

    //---------------------------------------------
    // this part is same for both Barrel and Endcap
    //---------------------------------------------
    // if this is not the first cell in the given layer then add the previous cell
    if (currentCellThetaBin > minCellThetaBin) {
      CellID nID = cID;
      decoder()->set(nID, m_thetaIndex, currentCellThetaBin - 1);
      cellNeighbours.push_back(nID); // add the previous cell from current layer of the same phi module
    }
    // if this is not the last cell in the given layer then add the next cell
    if (currentCellThetaBin < maxCellThetaBin) {
      CellID nID = cID;
      decoder()->set(nID, m_thetaIndex, currentCellThetaBin + 1);
      cellNeighbours.push_back(nID); // add the next cell from current layer of the same phi module
    }
    //----------------------------------------------

    // deal with the Barrel
    if (m_detLayout == 0) {
      const auto& edges = li.cellEdges.find(currentCellThetaBin)->second;
      double currentCellZmin = edges.first;
      double currentCellZmax = edges.second;

      // if this is not the first layer then look for neighbours in the previous layer
      if (currentLayerId > minLayerId) {
        CellID nID = cID;
        int prevLayerId = currentLayerId - 1;
        const LayerInfo& prevLI = getLayerInfo(prevLayerId);
        decoder()->set(nID, m_layerIndex, prevLayerId);

        decoder()->set(nID, m_thetaIndex, currentCellThetaBin);
        cellNeighbours.push_back(
            nID); // add the cell with the same theta bin from the previous layer of the same phi module

        // if the cID is in the positive-z side and prev layer cell is not in the first theta bin then add the cell from
        // previous theta bin
        if (theta(cID) < M_PI / 2. && currentCellThetaBin > prevLI.thetaBins.front()) {
          decoder()->set(nID, m_thetaIndex, currentCellThetaBin - 1);
          cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
          if (aDiagonal && currentCellThetaBin > (prevLI.thetaBins.front() + 1)) {
            // add the previous layer cell from the prev to prev theta bin if it overlaps with the current cell in
            // z-coordinate
            const auto& edgesm2 = prevLI.cellEdges.find(currentCellThetaBin - 2)->second;
            double zmin = edgesm2.first;
            if (zmin <= currentCellZmax) {
              // add the previous layer cell from the prev to prev theta bin
              decoder()->set(nID, m_thetaIndex, currentCellThetaBin - 2);
              cellNeighbours.push_back(nID);
            }
          }
        }
        // if the cID is in the negative-z side and prev layer cell is not in the last theta bin then add the cell from
        // previous theta bin
        if (theta(cID) > M_PI / 2. && currentCellThetaBin < prevLI.thetaBins.back()) {
          decoder()->set(nID, m_thetaIndex, currentCellThetaBin + 1);
          cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
          if (aDiagonal && currentCellThetaBin < (prevLI.thetaBins.back() - 1)) {
            // add the previous layer cell from the next to next theta bin if it overlaps with the current cell in
            // z-coordinate
            const auto& edgesp2 = prevLI.cellEdges.find(currentCellThetaBin + 2)->second;
            double zmax = edgesp2.second;
            if (zmax >= currentCellZmin) {
              // add the previous layer cell from the next to next theta bin
              decoder()->set(nID, m_thetaIndex, currentCellThetaBin + 2);
              cellNeighbours.push_back(nID);
            }
          }
        }
      }

      // if this is not the last layer then look for neighbours in the next layer
      if (currentLayerId < maxLayerId) {
        CellID nID = cID;
        int nextLayerId = currentLayerId + 1;
        const LayerInfo& nextLI = getLayerInfo(nextLayerId);
        decoder()->set(nID, m_layerIndex, nextLayerId);

        decoder()->set(nID, m_thetaIndex, currentCellThetaBin);
        cellNeighbours.push_back(
            nID); // add the cell with the same theta bin from the next layer of the same phi module

        // if the cID is in the positive-z side
        if (theta(cID) < M_PI / 2.) {
          // add the next layer cell from the next theta bin
          decoder()->set(nID, m_thetaIndex, currentCellThetaBin + 1);
          cellNeighbours.push_back(nID);

          if (aDiagonal) {
            // add the next layer cell from the next-to-next theta bin if it overlaps with the current cell in
            // z-coordinate
            const auto& edgesp2 = nextLI.cellEdges.find(currentCellThetaBin + 2)->second;
            double zmax = edgesp2.second;
            if (zmax >= currentCellZmin) {
              // add the next layer cell from the next to next theta bin
              decoder()->set(nID, m_thetaIndex, currentCellThetaBin + 2);
              cellNeighbours.push_back(nID);
            }
          }
        }
        // if the cID is in the negative-z side
        if (theta(cID) > M_PI / 2.) {
          // add the next layer cell from the previous theta bin
          decoder()->set(nID, m_thetaIndex, currentCellThetaBin - 1);
          cellNeighbours.push_back(nID);

          if (aDiagonal) {
            // add the next layer cell from the prev to prev theta bin if it overlaps with the current cell in
            // z-coordinate
            const auto& edgesm2 = nextLI.cellEdges.find(currentCellThetaBin - 2)->second;
            double zmin = edgesm2.first;
            if (zmin <= currentCellZmax) {
              // add the next layer cell from the prev to prev theta bin
              decoder()->set(nID, m_thetaIndex, currentCellThetaBin - 2);
              cellNeighbours.push_back(nID);
            }
          }
        }
      }
    }

    // Endcap
    if (m_detLayout == 1) {
      const auto& edges = li.cellEdges.find(currentCellThetaBin)->second;
      double currentCellZmin = edges.first;
      double currentCellZmax = edges.second;

      // if this is not the first layer then look for neighbours in the previous layer
      if (currentLayerId > minLayerId) {
        CellID nID = cID;
        int prevLayerId = currentLayerId - 1;
        const LayerInfo& prevLI = getLayerInfo(prevLayerId);
        decoder()->set(nID, m_layerIndex, prevLayerId);
        // find the ones that share at least part of a border with the current cell
        for (auto bin : prevLI.thetaBins) {
          const auto& edgesp = prevLI.cellEdges.find(bin)->second;
          double zmin = edgesp.first;
          double zmax = edgesp.second;

          // if the cID is in the positive-z side
          if (theta(cID) < M_PI / 2.) {
            if ((zmin >= currentCellZmin && zmin < currentCellZmax) ||
                (zmax >= currentCellZmin && zmax <= currentCellZmax) ||
                (currentCellZmin >= zmin && currentCellZmax <= zmax)) {
              decoder()->set(nID, m_thetaIndex, bin);
              cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
            }
            if (aDiagonal && zmin == currentCellZmax) {
              decoder()->set(nID, m_thetaIndex, bin);
              cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
            }
          }
          // if the cID is in the negative-z side
          if (theta(cID) > M_PI / 2.) {
            if ((zmin >= currentCellZmin && zmin <= currentCellZmax) ||
                (zmax > currentCellZmin && zmax <= currentCellZmax) ||
                (currentCellZmin >= zmin && currentCellZmax <= zmax)) {
              decoder()->set(nID, m_thetaIndex, bin);
              cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
            }
            if (aDiagonal && zmax == currentCellZmin) {
              decoder()->set(nID, m_thetaIndex, bin);
              cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
            }
          }
        }
      }
      // if this is not the last layer then look for neighbours in the next layer
      if (currentLayerId < maxLayerId) {
        CellID nID = cID;
        int nextLayerId = currentLayerId + 1;
        const LayerInfo& nextLI = getLayerInfo(nextLayerId);
        decoder()->set(nID, m_layerIndex, nextLayerId);
        // find the ones that share at least part of a border with the current cell
        for (auto bin : nextLI.thetaBins) {
          const auto& edgesn = nextLI.cellEdges.find(bin)->second;
          double zmin = edgesn.first;
          double zmax = edgesn.second;
          // if the cID is in the positive-z side
          if (theta(cID) < M_PI / 2.) {
            if ((zmin >= currentCellZmin && zmin <= currentCellZmax) ||
                (zmax > currentCellZmin && zmax <= currentCellZmax) ||
                (currentCellZmin >= zmin && currentCellZmax <= zmax)) {
              decoder()->set(nID, m_thetaIndex, bin);
              cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
            }
            if (aDiagonal && zmax == currentCellZmin) {
              decoder()->set(nID, m_thetaIndex, bin);
              cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
            }
          }
          // if the cID is in the negative-z side
          if (theta(cID) > M_PI / 2.) {
            if ((zmin >= currentCellZmin && zmin < currentCellZmax) ||
                (zmax >= currentCellZmin && zmax <= currentCellZmax) ||
                (currentCellZmin >= zmin && currentCellZmax <= zmax)) {
              decoder()->set(nID, m_thetaIndex, bin);
              cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
            }
            if (aDiagonal && zmin == currentCellZmax) {
              decoder()->set(nID, m_thetaIndex, bin);
              cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
            }
          }
        }
      }

      // if the Endcap consists of more than 1 part/section then look for neighbours in different parts as well
      if (m_offsetZ.size() > 1) {
        //
        double currentLayerRmin = li.radius - li.halfDepth;
        double currentLayerRmax = li.radius + li.halfDepth;

        // if the cell is in negative-z part, then swap min and max theta bins
        if (theta(cID) > M_PI / 2.) {
          minCellThetaBin = li.thetaBins.back();
          maxCellThetaBin = li.thetaBins[li.thetaBins.size() / 2];
        }

        // if it is the last cell in the part1
        if (EndcapPart == 0 && currentCellThetaBin == minCellThetaBin) {
          // find the layers in the part2 that share a border with the current layer
          for (int part2layerId = minLayerIdEndcap[1]; part2layerId <= maxLayerIdEndcap[1]; part2layerId++) {
            const LayerInfo& part2li = getLayerInfo(part2layerId);
            double Rmin = part2li.radius - part2li.halfDepth;
            double Rmax = part2li.radius + part2li.halfDepth;

            if ((Rmin >= currentLayerRmin && Rmin <= currentLayerRmax) ||
                (Rmax > currentLayerRmin && Rmax <= currentLayerRmax) ||
                (currentLayerRmin >= Rmin && currentLayerRmin < Rmax) ||
                (currentLayerRmax >= Rmin && currentLayerRmax <= Rmax)) {
              CellID nID = cID;
              decoder()->set(nID, m_layerIndex, part2layerId);
              decoder()->set(nID, m_thetaIndex, maxCellThetaBin);
              cellNeighbours.push_back(nID); // add the last theta bin cell from part2 layer
            }
            if (aDiagonal && Rmax == currentLayerRmin) {
              CellID nID = cID;
              decoder()->set(nID, m_layerIndex, part2layerId);
              decoder()->set(nID, m_thetaIndex, maxCellThetaBin);
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
              const LayerInfo& prevPartli = getLayerInfo(prevPartLayerId);
              double Rmin = prevPartli.radius - prevPartli.halfDepth;
              double Rmax = prevPartli.radius + prevPartli.halfDepth;

              if ((Rmin >= currentLayerRmin && Rmin < currentLayerRmax) ||
                  (Rmax >= currentLayerRmin && Rmax <= currentLayerRmax) ||
                  (currentLayerRmin >= Rmin && currentLayerRmin <= Rmax) ||
                  (currentLayerRmax > Rmin && currentLayerRmax <= Rmax)) {
                CellID nID = cID;
                decoder()->set(nID, m_layerIndex, prevPartLayerId);
                decoder()->set(nID, m_thetaIndex, minCellThetaBin);
                cellNeighbours.push_back(nID); // add the first theta bin cell from the part1 layer
              }
              if (aDiagonal && Rmin == currentLayerRmax) {
                CellID nID = cID;
                decoder()->set(nID, m_layerIndex, prevPartLayerId);
                decoder()->set(nID, m_thetaIndex, minCellThetaBin);
                cellNeighbours.push_back(nID); // add the first theta bin cell from the part1 layer
              }
            }
          }

          // if it is the first theta bin cell then look for neighbours in the next part
          if (currentCellThetaBin == minCellThetaBin) {
            // find the layers in the next part that share a border with the current layer
            for (int nextPartLayerId = minLayerIdEndcap[i_section + 1];
                 nextPartLayerId <= maxLayerIdEndcap[i_section + 1]; nextPartLayerId++) {
              const LayerInfo& nextPartli = getLayerInfo(nextPartLayerId);
              double Rmin = nextPartli.radius - nextPartli.halfDepth;
              double Rmax = nextPartli.radius + nextPartli.halfDepth;

              if ((Rmin >= currentLayerRmin && Rmin <= currentLayerRmax) ||
                  (Rmax > currentLayerRmin && Rmax <= currentLayerRmax) ||
                  (currentLayerRmin >= Rmin && currentLayerRmin < Rmax) ||
                  (currentLayerRmax >= Rmin && currentLayerRmax <= Rmax)) {
                CellID nID = cID;
                decoder()->set(nID, m_layerIndex, nextPartLayerId);
                decoder()->set(nID, m_thetaIndex, maxCellThetaBin);
                cellNeighbours.push_back(nID); // add the first cell from the part3 layer
              }
              if (aDiagonal && Rmax == currentLayerRmin) {
                CellID nID = cID;
                decoder()->set(nID, m_layerIndex, nextPartLayerId);
                decoder()->set(nID, m_thetaIndex, maxCellThetaBin);
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
            const LayerInfo& prevPartli = getLayerInfo(prevPartLayerId);
            double Rmin = prevPartli.radius - prevPartli.halfDepth;
            double Rmax = prevPartli.radius + prevPartli.halfDepth;

            if ((Rmin >= currentLayerRmin && Rmin < currentLayerRmax) ||
                (Rmax >= currentLayerRmin && Rmax <= currentLayerRmax) ||
                (currentLayerRmin >= Rmin && currentLayerRmin <= Rmax) ||
                (currentLayerRmax > Rmin && currentLayerRmax <= Rmax)) {
              CellID nID = cID;
              decoder()->set(nID, m_layerIndex, prevPartLayerId);
              decoder()->set(nID, m_thetaIndex, minCellThetaBin);
              cellNeighbours.push_back(nID); // add the first theta bin cell from the part2 layer
            }
            if (aDiagonal && Rmin == currentLayerRmax) {
              CellID nID = cID;
              decoder()->set(nID, m_layerIndex, prevPartLayerId);
              decoder()->set(nID, m_thetaIndex, minCellThetaBin);
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
  void FCCSWHCalPhiTheta_k4geo::neighbours(const CellID& cellID, std::set<CellID>& neighbours) const {
    std::vector<uint64_t> neigh = this->neighbours(cellID, false);
    neighbours.clear();
    neighbours.insert(neigh.begin(), neigh.end());
  }

  /// Determine minimum and maximum polar angle of the cell
  std::array<double, 2> FCCSWHCalPhiTheta_k4geo::cellTheta(const CellID cID) const {
    std::array<double, 2> cTheta = {M_PI, M_PI};

    // get the cell index
    int idx = decoder()->get(cID, m_thetaIndex);
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

} // namespace DDSegmentation
} // namespace dd4hep
