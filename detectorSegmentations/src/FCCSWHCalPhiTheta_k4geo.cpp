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
  registerParameter("offset_z", "Offset in z-axis of the layer center", m_offsetZ, std::vector<double>());
  registerParameter("width_z", "Width in z of the layer", m_widthZ, std::vector<double>());
  registerParameter("offset_r", "Offset in radius of the layer (Rmin)", m_offsetR, std::vector<double>());
  registerParameter("numLayers", "Number of layers", m_numLayers, std::vector<int>());
  registerParameter("dRlayer", "dR of the layer", m_dRlayer, std::vector<double>());
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
}

FCCSWHCalPhiTheta_k4geo::FCCSWHCalPhiTheta_k4geo(const BitFieldCoder* decoder) : GridTheta_k4geo(decoder) {
  // define type and description
  _type = "FCCSWHCalPhiTheta_k4geo";
  _description = "Phi-theta segmentation in the global coordinates";

  // register all necessary parameters (additional to those registered in GridTheta_k4geo)
  registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
  registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
  registerParameter("offset_z", "Offset in z-axis of the layer center", m_offsetZ, std::vector<double>());
  registerParameter("width_z", "Width in z of the layer", m_widthZ, std::vector<double>());
  registerParameter("offset_r", "Offset in radius of the layer (Rmin)", m_offsetR, std::vector<double>());
  registerParameter("numLayers", "Number of layers", m_numLayers, std::vector<int>());
  registerParameter("dRlayer", "dR of the layer", m_dRlayer, std::vector<double>());
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
}


/// determine the global position based on the cell ID
Vector3D FCCSWHCalPhiTheta_k4geo::position(const CellID& cID) const {
  uint layer = _decoder->get(cID,"layer");
  double radius = 1.0;

  if(m_radii.empty()) defineCellsInRZplan();
  if(!m_radii.empty()) radius = m_radii[layer];

  return positionFromRThetaPhi(radius, theta(cID), phi(cID));
}


/// determine the global position based on the cell ID
/// returns the geometric center of the cell
Vector3D FCCSWHCalPhiTheta_k4geo::centerPosition(const CellID& cID) const {
  uint layer = _decoder->get(cID,"layer");
  int thetaID = _decoder->get(cID,"theta");
  double zpos = 0.;
  double radius = 1.0;

  if(m_radii.empty()) defineCellsInRZplan();
  if(!m_radii.empty()) radius = m_radii[layer];
  if(!m_cellEdges.empty()) zpos = m_cellEdges[layer][thetaID].first + (m_cellEdges[layer][thetaID].second - m_cellEdges[layer][thetaID].first) * 0.5;

  auto pos = positionFromRThetaPhi(radius, theta(cID), phi(cID));

  // retrun the position with corrected z corrdinate to match to the geometric center
  return Vector3D(pos.x(),pos.y(),zpos);
}


void FCCSWHCalPhiTheta_k4geo::defineCellsInRZplan() const {
  if(m_radii.empty())
  {
    // check if all necessary variables are available
    if(m_offsetZ.empty() || m_widthZ.empty() || m_offsetR.empty() || m_numLayers.empty() || m_dRlayer.empty())
    {
       dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo","Please check the readout description in the XML file!",
                                       "One of the variables is missing: offset_z | width_z | offset_r | numLayers | dRlayer");
       return;
    }

    // calculate the radius for each layer
    if(m_offsetZ.size() == 1) // Barrel
    {
      dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiTheta_k4geo","Barrel configuration found!");
      uint N_dR = m_numLayers.size();
      double moduleDepth = 0.;
      for(uint i_dR = 0; i_dR < N_dR; i_dR++)
      {
       	for(int i_row = 1; i_row <= m_numLayers[i_dR]; i_row++)
        {
          moduleDepth+=m_dRlayer[i_dR];
          m_radii.push_back(m_offsetR[0] + moduleDepth - m_dRlayer[i_dR]*0.5);
          // layer lower and upper edges in z-axis
          m_layerEdges.push_back( std::make_pair(m_offsetZ[0] - 0.5*m_widthZ[0], m_offsetZ[0] + 0.5*m_widthZ[0]) );
        }
      }
    } // Barrel

    if(m_offsetZ.size() > 1) // ThreePartsEndCap
    {
      dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiTheta_k4geo","ThreePartsEndCap configuration found!");
      uint N_dR = m_numLayers.size()/3;
      double moduleDepth1 = 0.;
      double moduleDepth2 = 0.;
      double moduleDepth3 = 0.;
      // part1
      for(uint i_dR = 0; i_dR < N_dR; i_dR++)
      {
        for(int i_row = 1; i_row <= m_numLayers[i_dR]; i_row++)
        {
          moduleDepth1+=m_dRlayer[i_dR];
          m_radii.push_back(m_offsetR[0] + moduleDepth1 - m_dRlayer[i_dR]*0.5);
          // layer lower and upper edges in z-axis
          m_layerEdges.push_back( std::make_pair(m_offsetZ[0] - 0.5*m_widthZ[0], m_offsetZ[0] + 0.5*m_widthZ[0]) );
        }
      }
      // part2
      for(uint i_dR = 0; i_dR < N_dR; i_dR++)
      {
       	for(int i_row = 1; i_row <= m_numLayers[i_dR + N_dR]; i_row++)
        {
          moduleDepth2+=m_dRlayer[i_dR];
          m_radii.push_back(m_offsetR[1] + moduleDepth2 - m_dRlayer[i_dR]*0.5);
          // layer lower and upper edges in z-axis
          m_layerEdges.push_back( std::make_pair(m_offsetZ[1] - 0.5*m_widthZ[1], m_offsetZ[1] + 0.5*m_widthZ[1]) );
        }
      }
      // part3
      for(uint i_dR = 0; i_dR < N_dR; i_dR++)
      {
       	for(int i_row = 1; i_row <= m_numLayers[i_dR + 2*N_dR]; i_row++)
        {
          moduleDepth3+=m_dRlayer[i_dR];
          m_radii.push_back(m_offsetR[2] + moduleDepth3 - m_dRlayer[i_dR]*0.5);
          // layer lower and upper edges in z-axis
          m_layerEdges.push_back( std::make_pair(m_offsetZ[2] - 0.5*m_widthZ[2], m_offsetZ[2] + 0.5*m_widthZ[2]) );
        }
      }
    } // ThreePartsEndcap

    // print info of calculated radii and edges
    for(uint i_layer = 0; i_layer < m_radii.size(); i_layer++){
       dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiTheta_k4geo","layer %d radius: %.2f, z range: %.2f - %.2f cm", 
                                      i_layer, m_radii[i_layer], m_layerEdges[i_layer].first, m_layerEdges[i_layer].second);
    }

    // allocate thetaBins vector for each layer
    m_thetaBins.resize(m_radii.size());
    // allocate cellEdges vector for each layer
    m_cellEdges.resize(m_radii.size());

    // determine theta bins and cell edges for each layer
    for(uint i_layer = 0; i_layer < m_radii.size(); i_layer++) defineCellEdges(i_layer);
  }
}


void FCCSWHCalPhiTheta_k4geo::defineCellEdges(const uint layer) const {
  if(m_thetaBins[layer].size()==0 && m_radii.size() > 0)
  {
    //find theta bins that fit within the given layer
    int ibin = positionToBin(0.02, m_gridSizeTheta, m_offsetTheta); // <--- start from theta bin outside the HCal theta range
    while(m_radii[layer]*std::cos(m_offsetTheta+ibin*m_gridSizeTheta)/std::sin(m_offsetTheta+ibin*m_gridSizeTheta) > m_layerEdges[layer].first)
    {
      if(m_radii[layer]*std::cos(m_offsetTheta+ibin*m_gridSizeTheta)/std::sin(m_offsetTheta+ibin*m_gridSizeTheta) < m_layerEdges[layer].second)
      {
        m_thetaBins[layer].push_back(ibin);
      }
      ibin++;
    }

    // find edges of each cell (theta bin) in the given layer
    auto prevBin = m_thetaBins[layer][0];
    // set the upper edge of the first cell in the given layer (starting from positive z part)
    m_cellEdges[layer][prevBin] = std::make_pair(0.,m_layerEdges[layer].second);
    for(auto bin : m_thetaBins[layer])
    {
      if(bin!=prevBin)
      {
        double z1 = m_radii[layer]*std::cos(m_offsetTheta+bin*m_gridSizeTheta)/std::sin(m_offsetTheta+bin*m_gridSizeTheta);
        double z2 = m_radii[layer]*std::cos(m_offsetTheta+prevBin*m_gridSizeTheta)/std::sin(m_offsetTheta+prevBin*m_gridSizeTheta);
        // set the lower edge of the prevBin cell
        m_cellEdges[layer][prevBin].first = z1 + 0.5*(z2 - z1);
        // set the upper edge of current bin cell
        m_cellEdges[layer][bin] = std::make_pair(0.,m_cellEdges[layer][prevBin].first);
        prevBin = bin;
      }
    }
    // set the lower edge of the last cell in the given layer
    m_cellEdges[layer][prevBin].first = m_layerEdges[layer].first;

    // for the EndCap, do it again but for negative z part
    if(m_offsetZ.size() > 1)
    {
      while(m_radii[layer]*std::cos(m_offsetTheta+ibin*m_gridSizeTheta)/std::sin(m_offsetTheta+ibin*m_gridSizeTheta) > (-m_layerEdges[layer].second))
      {
       	if(m_radii[layer]*std::cos(m_offsetTheta+ibin*m_gridSizeTheta)/std::sin(m_offsetTheta+ibin*m_gridSizeTheta) < (-m_layerEdges[layer].first))
        {
          m_thetaBins[layer].push_back(ibin);
        }
	ibin++;
      }

      // Create a span view over the theta bins corresponding to the Endcap in negative z part
      std::span<int> thetaBins(m_thetaBins[layer].begin() + m_thetaBins[layer].size()/2, m_thetaBins[layer].end());
      prevBin = thetaBins[0];

      // set the upper edge of the first cell in the given layer at negative z part
      m_cellEdges[layer][prevBin] = std::make_pair(0.,-m_layerEdges[layer].first);
      for(auto bin : thetaBins)
      {
        if(bin!=prevBin)
        {
          double z1 = m_radii[layer]*std::cos(m_offsetTheta+bin*m_gridSizeTheta)/std::sin(m_offsetTheta+bin*m_gridSizeTheta);
          double z2 = m_radii[layer]*std::cos(m_offsetTheta+prevBin*m_gridSizeTheta)/std::sin(m_offsetTheta+prevBin*m_gridSizeTheta);
          // set the lower edge of the prevBin cell
          m_cellEdges[layer][prevBin].first = z1 + 0.5*(z2 - z1);
          // set the upper edge of current bin cell
          m_cellEdges[layer][bin] = std::make_pair(0.,m_cellEdges[layer][prevBin].first);
          prevBin = bin;
        }
      }
      // set the lower edge of the last cell in the given layer
      m_cellEdges[layer][prevBin].first = (-m_layerEdges[layer].second);
    }// negative-z endcap

    dd4hep::printout(dd4hep::DEBUG, "FCCSWHCalPhiTheta_k4geo","Number of cells in layer %d: %d", layer, m_thetaBins[layer].size());
    for(auto bin : m_thetaBins[layer]) dd4hep::printout(dd4hep::DEBUG, "FCCSWHCalPhiTheta_k4geo","Layer %d cell theta bin: %d, edges: %.2f - %.2f cm", layer, bin, m_cellEdges[layer][bin]);
  }
}

/// create the cell ID based on the position
CellID FCCSWHCalPhiTheta_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                          const VolumeID& vID) const {

  CellID cID = vID;

  // The volume ID comes with "row" field information (number of sequences) that would screw up the topo-clustering using cell
  // neighbours map produced with RecFCCeeCalorimeter/src/components/CreateFCCeeCaloNeighbours.cpp,
  // therefore, lets set it to zero, as it is for the cell IDs in the neighbours map.
  _decoder->set(cID, "row", 0);

  // For endcap, the volume ID comes with "type" field information which would screw up the topo-clustering as the "row" field,
  // therefore, lets set it to zero, as it is for the cell IDs in the neighbours map.
  if(m_offsetZ.size() > 1) _decoder->set(cID, "type", 0);

  double lTheta = thetaFromXYZ(globalPosition);
  double lPhi = phiFromXYZ(globalPosition);
  uint layer = _decoder->get(vID,"layer");

  // define cell boundaries in R-z plan
  if(m_radii.empty()) defineCellsInRZplan();

  // check if the cells are defined for the given layer
  if(m_thetaBins[layer].empty()) dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiTheta_k4geo","No cells are defined for layer %d", layer);

  // find the cell (theta bin) corresponding to the hit and return the cellID
  for(auto bin : m_thetaBins[layer])
  {
    double posz =  globalPosition.z();
    if(posz > m_cellEdges[layer][bin].first && posz < m_cellEdges[layer][bin].second)
    {
      _decoder->set(cID, m_thetaID, bin);
      _decoder->set(cID, m_phiID, positionToBin(lPhi, 2 * M_PI / (double)m_phiBins, m_offsetPhi));
      return cID;
    }
  }

  dd4hep::printout(dd4hep::WARNING, "FCCSWHCalPhiTheta_k4geo","The hit is outside the defined range of the layer %d", layer);

  _decoder->set(cID, m_thetaID, positionToBin(lTheta, m_gridSizeTheta, m_offsetTheta));
  _decoder->set(cID, m_phiID, positionToBin(lPhi, 2 * M_PI / (double)m_phiBins, m_offsetPhi));
  return cID;
}

/// determine the azimuthal angle phi based on the cell ID
double FCCSWHCalPhiTheta_k4geo::phi(const CellID& cID) const {
  CellID phiValue = _decoder->get(cID, m_phiID);
  return binToPosition(phiValue, 2. * M_PI / (double)m_phiBins, m_offsetPhi);
}
}
}
