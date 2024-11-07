#include "detectorSegmentations/FCCSWHCalPhiRow_k4geo.h"
#include "DD4hep/Printout.h"

namespace dd4hep {
namespace DDSegmentation {

using std::runtime_error;


/// default constructor using an encoding string
FCCSWHCalPhiRow_k4geo::FCCSWHCalPhiRow_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) {
  // define type and description
  _type = "FCCSWHCalPhiRow_k4geo";
  _description = "Phi-theta segmentation in the global coordinates";

  // register all necessary parameters
  registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
  registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
  registerParameter("grid_size_row", "Cell size in row", m_gridSizeRow, std::vector<int>());
  registerParameter("dz_row", "dz of row", m_dz_row, 0.);
  registerParameter("offset_z", "Offset in z-axis of the layer center", m_offsetZ, std::vector<double>());
  registerParameter("width_z", "Width in z of the layer", m_widthZ, std::vector<double>());
  registerParameter("offset_r", "Offset in radius of the layer (Rmin)", m_offsetR, std::vector<double>());
  registerParameter("numLayers", "Number of layers", m_numLayers, std::vector<int>());
  registerParameter("dRlayer", "dR of the layer", m_dRlayer, std::vector<double>());
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
  registerIdentifier("identifier_row", "Cell ID identifier for row", m_rowID, "row");
}

FCCSWHCalPhiRow_k4geo::FCCSWHCalPhiRow_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) {
  // define type and description
  _type = "FCCSWHCalPhiRow_k4geo";
  _description = "Phi-theta segmentation in the global coordinates";

  // register all necessary parameters
  registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
  registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
  registerParameter("grid_size_row", "Cell size in row", m_gridSizeRow, std::vector<int>());
  registerParameter("dz_row", "dz of row", m_dz_row, 0.);
  registerParameter("offset_z", "Offset in z-axis of the layer center", m_offsetZ, std::vector<double>());
  registerParameter("width_z", "Width in z of the layer", m_widthZ, std::vector<double>());
  registerParameter("offset_r", "Offset in radius of the layer (Rmin)", m_offsetR, std::vector<double>());
  registerParameter("numLayers", "Number of layers", m_numLayers, std::vector<int>());
  registerParameter("dRlayer", "dR of the layer", m_dRlayer, std::vector<double>());
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
  registerIdentifier("identifier_row", "Cell ID identifier for row", m_rowID, "row");
}


/// determine the global position based on the cell ID
Vector3D FCCSWHCalPhiRow_k4geo::position(const CellID& cID) const {
  uint layer = _decoder->get(cID,"layer");
  double phi = _decoder->get(cID,m_phiID);

  if(m_radii.empty()) calculateLayerRadii();
  if(m_radii.empty() || m_layerEdges.empty())
  {
    dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiRow_k4geo","Could not calculate layer radii!");
    return Vector3D(0.,0.,0.);
  }

  double radius = m_radii[layer];
  double minLayerZ = m_layerEdges[layer].first;

  // get index of the cell in the layer (index starts from 1!)
  int idx = _decoder->get(cID, "row");
  // calculate z-coordinate of the cell center
  double zpos = minLayerZ + (idx-1) * m_dz_row * m_gridSizeRow[layer] + 0.5 * m_dz_row * m_gridSizeRow[layer];

  // for negative-z Endcap, the index is negetive (starts from -1!)
  if(idx < 0) zpos = -minLayerZ + (idx+1) * m_dz_row * m_gridSizeRow[layer] - 0.5 * m_dz_row * m_gridSizeRow[layer];

  return Vector3D( radius * std::cos(phi), radius * std::sin(phi), zpos );
}


void FCCSWHCalPhiRow_k4geo::calculateLayerRadii() const {
  if(m_radii.empty())
  {
    // check if all necessary variables are available
    if(m_offsetZ.empty() || m_widthZ.empty() || m_offsetR.empty() || m_numLayers.empty() || m_dRlayer.empty())
    {
       dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiRow_k4geo","Please check the readout description in the XML file!",
                                       "One of the variables is missing: offset_z | width_z | offset_r | numLayers | dRlayer");
       return;
    }

    // calculate the radius for each layer
    if(m_offsetZ.size() == 1) // Barrel
    {
      dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiRow_k4geo","Barrel configuration found!");
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
          m_layerDepth.push_back(m_dRlayer[i_dR]);
        }
      }
    } // Barrel

    if(m_offsetZ.size() > 1) // ThreePartsEndCap
    {
      dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiRow_k4geo","ThreePartsEndCap configuration found!");
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
          m_layerDepth.push_back(m_dRlayer[i_dR]);
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
          m_layerDepth.push_back(m_dRlayer[i_dR]);
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
          m_layerDepth.push_back(m_dRlayer[i_dR]);
        }
      }
    } // ThreePartsEndcap

    // print info of calculated radii and edges
    for(uint i_layer = 0; i_layer < m_radii.size(); i_layer++){
       dd4hep::printout(dd4hep::INFO, "FCCSWHCalPhiRow_k4geo","layer %d radius: %.2f, z range: %.2f - %.2f cm", 
                                      i_layer, m_radii[i_layer], m_layerEdges[i_layer].first, m_layerEdges[i_layer].second);
    }

    // allocate cellIndexes vector for each layer
    m_cellIndexes.resize(m_radii.size());
    // allocate cellEdges vector for each layer
    m_cellEdges.resize(m_radii.size());

    // determine row bins for each layer
    for(uint i_layer = 0; i_layer < m_radii.size(); i_layer++) defineCellIndexes(i_layer);
  }
}

/*
*  This function fills the m_cellIndexes vector per layer with the cell indexes.
*  The cell index is encoded in CellID with "row" field.
*  In case of a cell with single row/sequence, the index is directly the number of row in the layer.
*  In case of a cell with several rows/sequences merged, the index is the number of cell in the layer.
*  For the layers of negative-z Endcap, indexes of cells are negative (-1, -2, ..., -N).
*
*  m_cellIndexes vector can be used to define the neighbours using CreateFCCeeCaloNeighbours tool.
*/
void FCCSWHCalPhiRow_k4geo::defineCellIndexes(const uint layer) const {
  if(m_cellIndexes[layer].size()==0 && m_radii.size() > 0)
  {
    double minLayerZ = m_layerEdges[layer].first;
    double maxLayerZ = m_layerEdges[layer].second;

    // find rows/sequences that fit within the given layer range along z
    int irow = 0;
    while( (minLayerZ + (irow+1)* m_dz_row) < (maxLayerZ+0.0001) )
    {
      // define the cell index
      int idx = floor(irow/m_gridSizeRow[layer]) + 1;
      // add the index if it is not already there
      if(std::find(m_cellIndexes[layer].begin(),m_cellIndexes[layer].end(),idx) == m_cellIndexes[layer].end()) m_cellIndexes[layer].push_back(idx);
      irow++;
    }

    // for the EndCap, do it again but for negative-z part
    if(m_offsetZ.size() > 1)
    {
      irow = 0;
      while( (minLayerZ + (irow+1) * m_dz_row) < (maxLayerZ+0.0001) )
      {
        // define the cell index with negative sign
        int idx = - ( floor(irow/m_gridSizeRow[layer]) + 1 );
        // add the index if it is not already there
        if(std::find(m_cellIndexes[layer].begin(),m_cellIndexes[layer].end(),idx) == m_cellIndexes[layer].end()) m_cellIndexes[layer].push_back(idx);
        irow++;
      }
    }

    // find edges of each cell in the given layer along z axis
    for(auto idx : m_cellIndexes[layer])
    {
      // calculate z-coordinates of the cell edges
      double z1 = minLayerZ + (idx-1) * m_dz_row * m_gridSizeRow[layer]; // lower edge
      double z2 = minLayerZ + (idx) * m_dz_row * m_gridSizeRow[layer]; // upper edge

      // for negative-z Endcap, the index is negetive (starts from -1!)
      if(idx < 0)
      { 
        z1 = -minLayerZ + (idx) * m_dz_row * m_gridSizeRow[layer]; // lower edge
        z2 = -minLayerZ + (idx+1) * m_dz_row * m_gridSizeRow[layer]; // upper edge
      }

      m_cellEdges[layer][idx] = std::make_pair(z1,z2);
    }

    dd4hep::printout(dd4hep::DEBUG, "FCCSWHCalPhiRow_k4geo","Number of cells in layer %d: %d", layer, m_cellIndexes[layer].size());
  }
}

/// create the cell ID based on the position
CellID FCCSWHCalPhiRow_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                          const VolumeID& vID) const {

  // get the row number from volumeID (starts from 0!)
  int nrow = _decoder->get(vID, "row");
  // get the layer number from volumeID
  uint layer = _decoder->get(vID,"layer");

  CellID cID = vID;

  // get the cell index (start from 1!)
  int idx = floor(nrow/m_gridSizeRow[layer]) + 1;

  // if the hit is in the negative-z part of the Endcap then assign negative index
  if(m_offsetZ.size() > 1 && globalPosition.z() < 0) idx *= -1;

  _decoder->set(cID, m_rowID, idx);
  _decoder->set(cID, m_phiID, positionToBin(dd4hep::DDSegmentation::Util::phiFromXYZ(globalPosition), 2 * M_PI / (double)m_phiBins, m_offsetPhi));

  // For endcap, the volume ID comes with "type" field information which would screw up the topo-clustering,
  // therefore, lets set it to zero, as it is for the cell IDs in the neighbours map.
  if(m_offsetZ.size() > 1) _decoder->set(cID, "type", 0);

  return cID;
}

/// determine the azimuthal angle phi based on the cell ID
double FCCSWHCalPhiRow_k4geo::phi(const CellID& cID) const {
  CellID phiValue = _decoder->get(cID, m_phiID);
  return binToPosition(phiValue, 2. * M_PI / (double)m_phiBins, m_offsetPhi);
}

/// Get the min and max layer indexes for each part of the HCal
std::vector<std::pair<uint,uint> > FCCSWHCalPhiRow_k4geo::getMinMaxLayerId() const {

   std::vector<std::pair<uint,uint> > minMaxLayerId;

   if(m_radii.empty()) calculateLayerRadii();
   if(m_radii.empty()) return minMaxLayerId;


   std::vector<int> minLayerIdEndcap = {0, 0, 0};
   std::vector<int> maxLayerIdEndcap = {-1, 0, 0};
   if(m_offsetZ.size() > 1)
   {
     uint Nl = m_numLayers.size()/3;

     // count the number of layers in the first part of the Endcap
     for(uint i=0; i < Nl; i++) maxLayerIdEndcap[0] += m_numLayers[i];

     minLayerIdEndcap[1] = maxLayerIdEndcap[0]+1;
     maxLayerIdEndcap[1] = maxLayerIdEndcap[0];
     // count the number of layers in the second part of the Endcap
     for(uint i=0; i < Nl; i++) maxLayerIdEndcap[1] += m_numLayers[i+Nl];

     minLayerIdEndcap[2] = maxLayerIdEndcap[1]+1;
     maxLayerIdEndcap[2] = maxLayerIdEndcap[1];
     // count the number of layers in the third part of the Endcap
     for(uint i=0; i < Nl; i++) maxLayerIdEndcap[2] += m_numLayers[i+2*Nl];

     minMaxLayerId.push_back(std::make_pair(minLayerIdEndcap[0],maxLayerIdEndcap[0]));
     minMaxLayerId.push_back(std::make_pair(minLayerIdEndcap[1],maxLayerIdEndcap[1]));
     minMaxLayerId.push_back(std::make_pair(minLayerIdEndcap[2],maxLayerIdEndcap[2]));
   }
   else // for Barrel
   {
     minMaxLayerId.push_back(std::make_pair(0,m_radii.size()-1));
   }

   return minMaxLayerId;
}


/// Calculates the neighbours of the given cell ID and adds them to the list of neighbours
std::vector<uint64_t> FCCSWHCalPhiRow_k4geo::neighbours(const CellID& cID) const {
   std::vector<uint64_t> cellNeighbours;

   if(m_radii.empty()) calculateLayerRadii();
   if(m_cellIndexes.empty()) return cellNeighbours;

   bool EndcapPart1 = false;
   bool EndcapPart2 = false;
   bool EndcapPart3 = false;

   int minLayerId = -1;
   int maxLayerId = -1;

   int currentLayerId = _decoder->get(cID,"layer");
   int currentCellId = _decoder->get(cID,m_rowID);

   int minCellId = m_cellIndexes[currentLayerId].front();
   int maxCellId = m_cellIndexes[currentLayerId].back();

   //--------------------------------
   // Determine min and max layer Id
   //--------------------------------
   // if this is the segmentation of three parts Endcap
   /*
    *  hardcoded numbers would be the following:
    *  std::vector<int> minLayerIdEndcap = {0, 6, 15};
    *  std::vector<int> maxLayerIdEndcap = {5, 14, 36};
    *  but lets try to avoid hardcoding:
    */
   std::vector<int> minLayerIdEndcap = {0, 0, 0};
   std::vector<int> maxLayerIdEndcap = {0, 0, 0};
   if(m_offsetZ.size() > 1)
   {
     std::vector<std::pair<uint,uint> > minMaxLayerId(getMinMaxLayerId());
     if(minMaxLayerId.empty())
     {
       dd4hep::printout(dd4hep::ERROR, "FCCSWHCalPhiRow_k4geo","Can not get ThreePartsEndcap min and max layer indexes! --> returning empty neighbours");
       return cellNeighbours;
     }

     // part1 min and max layer index
     minLayerIdEndcap[0] = minMaxLayerId[0].first;
     maxLayerIdEndcap[0] = minMaxLayerId[0].second;
     // part2 min and max layer index
     minLayerIdEndcap[1] = minMaxLayerId[1].first;
     maxLayerIdEndcap[1] = minMaxLayerId[1].second;
     // part3 min and max layer index
     minLayerIdEndcap[2] = minMaxLayerId[2].first;
     maxLayerIdEndcap[2] = minMaxLayerId[2].second;

     // Part 1
     if(currentLayerId >= minLayerIdEndcap[0] && currentLayerId <= maxLayerIdEndcap[0])
     {
       minLayerId = minLayerIdEndcap[0];
       maxLayerId = maxLayerIdEndcap[0];
       EndcapPart1 = true;
     }
     // Part 2
     if(currentLayerId >= minLayerIdEndcap[1] && currentLayerId <= maxLayerIdEndcap[1])
     {
       minLayerId = minLayerIdEndcap[1];
       maxLayerId = maxLayerIdEndcap[1];
       EndcapPart2 = true;
     }
     // Part 3
     if(currentLayerId >= minLayerIdEndcap[2] && currentLayerId <= maxLayerIdEndcap[2])
     {
       minLayerId = minLayerIdEndcap[2];
       maxLayerId = maxLayerIdEndcap[2];
       EndcapPart3 = true;
     }

     // correct the min and max CellId for endcap
     if(currentCellId < 0) // negative-z part (cell index is negative (-1, -2, ... -N))
     {
       // second half of elements in m_cellIndexes[currentLayerId] vector corresponds to the negative-z layer cells
       minCellId = m_cellIndexes[currentLayerId].back();
       maxCellId = m_cellIndexes[currentLayerId][m_cellIndexes[currentLayerId].size()/2];
     }
     else // positive-z part (cell index is positive (1, 2, ... N))
     {
       // first half of elements in m_cellIndexes[currentLayerId] vector corresponds to the positive-z layer cells
       minCellId = m_cellIndexes[currentLayerId].front();
       maxCellId = m_cellIndexes[currentLayerId][m_cellIndexes[currentLayerId].size()/2 - 1];
     }
   }
   else // for Barrel
   {
     minLayerId = 0;
     maxLayerId = m_radii.size()-1;
   }
   //--------------------------------


   //--------------------------------
   // Find neighbours
   //--------------------------------

   // if this is not the first layer then look for neighbours in the previous layer
   if(currentLayerId > minLayerId)
   {
     CellID nID = cID ;
     int prevLayerId = currentLayerId - 1;
     _decoder->set(nID,"layer",prevLayerId);

     // if the granularity is the same for the previous layer then take the cells with currentCellId, currentCellId - 1, and currentCellId + 1
     if(m_gridSizeRow[prevLayerId] == m_gridSizeRow[currentLayerId])
     {
       _decoder->set(nID,m_rowID,currentCellId);
       cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
       if(currentCellId > minCellId)
       {
         _decoder->set(nID,m_rowID,currentCellId - 1);
         cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
       }
       if(currentCellId < maxCellId)
       {
         _decoder->set(nID,m_rowID,currentCellId + 1);
         cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
       }
     }
     // if the granularity is different
     else
     {
       // determine the cell index in the previous layer that is below of the current cell
       int idx = (currentCellId > 0) ? ((currentCellId-1)/m_gridSizeRow[prevLayerId] + 1) : ((currentCellId+1)/m_gridSizeRow[prevLayerId] - 1);
       _decoder->set(nID,m_rowID,idx);
       cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module

       //
       if((m_gridSizeRow[prevLayerId] - abs(currentCellId)%m_gridSizeRow[prevLayerId]) == (m_gridSizeRow[prevLayerId]-1) && currentCellId > minCellId)
       {
         _decoder->set(nID,m_rowID, (idx > 0) ? (idx - 1) : (idx + 1));
         cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
       }

       //
       if(abs(currentCellId)%m_gridSizeRow[prevLayerId] == 0 && currentCellId < maxCellId)
       {
         _decoder->set(nID,m_rowID, (idx > 0) ? (idx + 1) : (idx - 1));
         cellNeighbours.push_back(nID); // add the cell from the previous layer of the same phi module
       }
     }
   }

   // if this is not the last layer then look for neighbours in the next layer
   if(currentLayerId < maxLayerId)
   {
     CellID nID = cID ;
     int nextLayerId = currentLayerId + 1;
     _decoder->set(nID,"layer",nextLayerId);

     // if the granularity is the same for the next layer then take the cells with currentCellId, currentCellId - 1, and currentCellId + 1
     if(m_gridSizeRow[nextLayerId] == m_gridSizeRow[currentLayerId])
     {
       _decoder->set(nID,m_rowID,currentCellId);
       cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
       if(currentCellId > minCellId)
       {
         _decoder->set(nID,m_rowID,currentCellId - 1);
         cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
       }
       if(currentCellId < maxCellId)
       {
         _decoder->set(nID,m_rowID,currentCellId + 1);
         cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
       }
     }
     // if the granularity is different
     else
     {
       // determine the cell index in the next layer that is below of the current cell
       int idx = (currentCellId > 0) ? ((currentCellId-1)/m_gridSizeRow[nextLayerId] + 1) : ((currentCellId+1)/m_gridSizeRow[nextLayerId] - 1);
       _decoder->set(nID,m_rowID,idx);
       cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module

       //
       if((m_gridSizeRow[nextLayerId] - abs(currentCellId)%m_gridSizeRow[nextLayerId]) == (m_gridSizeRow[nextLayerId]-1) && currentCellId > minCellId)
       {
         _decoder->set(nID,m_rowID, (idx > 0) ? (idx - 1) : (idx + 1));
         cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
       }

       //
       if(abs(currentCellId)%m_gridSizeRow[nextLayerId] == 0 && currentCellId < maxCellId)
       {
         _decoder->set(nID,m_rowID, (idx > 0) ? (idx + 1) : (idx - 1));
         cellNeighbours.push_back(nID); // add the cell from the next layer of the same phi module
       }
     }
   }

   // if this is not the first cell in the given layer then add the previous cell
   if(currentCellId > minCellId)
   {
     CellID nID = cID ;
     _decoder->set(nID,m_rowID,currentCellId - 1);
     cellNeighbours.push_back(nID); // add the previous cell from current layer of the same phi module
   }
   // if this is not the last cell in the given layer then add the next cell
   if(currentCellId < maxCellId)
   {
     CellID nID = cID ;
     _decoder->set(nID,m_rowID,currentCellId + 1);
     cellNeighbours.push_back(nID); // add the next cell from current layer of the same phi module
   }

   // if this is the Endcap then look for neighbours in different parts as well
   if(m_offsetZ.size() > 1)
   {
     double currentLayerRmin = m_radii[currentLayerId] - 0.5*m_layerDepth[currentLayerId];
     double currentLayerRmax = m_radii[currentLayerId] + 0.5*m_layerDepth[currentLayerId];

     // if the cell is in negative-z part, then swap min and max cell indexes
     if(currentCellId < 0)
     {
       minCellId = m_cellIndexes[currentLayerId][m_cellIndexes[currentLayerId].size()/2]; // this should be -1
       maxCellId = m_cellIndexes[currentLayerId].back(); // this should be -N
     }

     // if it is the last cell in the part1
     if(EndcapPart1 && currentCellId == maxCellId )
     {
       // find the layers in the part2 that share a border with the current layer
       for(int part2layerId = minLayerIdEndcap[1]; part2layerId <= maxLayerIdEndcap[1]; part2layerId++)
       {
         double Rmin = m_radii[part2layerId] - 0.5*m_layerDepth[part2layerId];
         double Rmax = m_radii[part2layerId] + 0.5*m_layerDepth[part2layerId];

         if( (Rmin >= currentLayerRmin && Rmin <= currentLayerRmax)
          || (Rmax >= currentLayerRmin && Rmax <= currentLayerRmax)
          || (currentLayerRmin >= Rmin && currentLayerRmin <= Rmax)
          || (currentLayerRmax >= Rmin && currentLayerRmax <= Rmax)
         )
         {
           CellID nID = cID ;
           _decoder->set(nID,"layer",part2layerId);
           _decoder->set(nID,m_rowID,minCellId);
           cellNeighbours.push_back(nID); // add the first cell from part2 layer
         }
       }
     }

     // if it is the first cell in the part2
     if(EndcapPart2 && currentCellId == minCellId)
     {
       // find the layers in part1 that share a border with the current layer
       for(int part1layerId = minLayerIdEndcap[0]; part1layerId <= maxLayerIdEndcap[0]; part1layerId++)
       {
         double Rmin = m_radii[part1layerId] - 0.5*m_layerDepth[part1layerId];
         double Rmax = m_radii[part1layerId] + 0.5*m_layerDepth[part1layerId];

         if( (Rmin >= currentLayerRmin && Rmin <= currentLayerRmax)
          || (Rmax >= currentLayerRmin && Rmax <= currentLayerRmax)
          || (currentLayerRmin >= Rmin && currentLayerRmin <= Rmax)
          || (currentLayerRmax >= Rmin && currentLayerRmax <= Rmax)
         )
         {
           CellID nID = cID ;
           _decoder->set(nID,"layer",part1layerId);
           _decoder->set(nID,m_rowID,maxCellId);
           cellNeighbours.push_back(nID); // add the last cell from the part1 layer
         }
       }
     }

     // if it is the last cell in the part2
     if(EndcapPart2 && currentCellId == maxCellId)
     {
       // find the layers in the part3 that share a border with the current layer
       for(int part3layerId = minLayerIdEndcap[2]; part3layerId <= maxLayerIdEndcap[2]; part3layerId++)
       {
         double Rmin = m_radii[part3layerId] - 0.5*m_layerDepth[part3layerId];
         double Rmax = m_radii[part3layerId] + 0.5*m_layerDepth[part3layerId];

         if( (Rmin >= currentLayerRmin && Rmin <= currentLayerRmax)
          || (Rmax >= currentLayerRmin && Rmax <= currentLayerRmax)
          || (currentLayerRmin >= Rmin && currentLayerRmin <= Rmax)
          || (currentLayerRmax >= Rmin && currentLayerRmax <= Rmax)
         )
         {
           CellID nID = cID ;
           _decoder->set(nID,"layer",part3layerId);
           _decoder->set(nID,m_rowID,minCellId);
           cellNeighbours.push_back(nID); // add the first cell from the part3 layer
         }
       }
     }

     // if it is the first cell in the part3
     if(EndcapPart3 && currentCellId == minCellId)
     {
       // find the layers in the part2 that share a border with the current layer
       for(int part2layerId = minLayerIdEndcap[1]; part2layerId <= maxLayerIdEndcap[1]; part2layerId++)
       {
         double Rmin = m_radii[part2layerId] - 0.5*m_layerDepth[part2layerId];
         double Rmax = m_radii[part2layerId] + 0.5*m_layerDepth[part2layerId];

         if( (Rmin >= currentLayerRmin && Rmin <= currentLayerRmax)
          || (Rmax >= currentLayerRmin && Rmax <= currentLayerRmax)
          || (currentLayerRmin >= Rmin && currentLayerRmin <= Rmax)
          || (currentLayerRmax >= Rmin && currentLayerRmax <= Rmax)
         )
         {
           CellID nID = cID ;
           _decoder->set(nID,"layer",part2layerId);
           _decoder->set(nID,m_rowID,maxCellId);
           cellNeighbours.push_back(nID); // add the last cell from the part2 layer
         }
       }
     }
   }

   // Now loop over the neighbours and add the cells from next/previous phi module
   for(auto nID : cellNeighbours)
   {
     CellID newID = nID;
     // previous: if the current is 0 then previous is the last bin (id = m_phiBins - 1) else current - 1
     _decoder->set(newID,m_phiID, (_decoder->get(nID,m_phiID) == 0) ? m_phiBins - 1 : _decoder->get(nID,m_phiID) - 1);
     cellNeighbours.push_back(newID);
     // next: if the current is the last bin (id = m_phiBins - 1) then the next is the first bin (id = 0) else current + 1
     _decoder->set(newID,m_phiID, (_decoder->get(nID,m_phiID) == (m_phiBins - 1)) ? 0 : _decoder->get(nID,m_phiID) + 1);
     cellNeighbours.push_back(newID);
   }

   // At the end, find neighbours with the same layer/row in next/previous phi module
   CellID nID = cID ;
   // previous: if the current is 0 then previous is the last bin (id = m_phiBins - 1) else current - 1
   _decoder->set(nID,m_phiID, (_decoder->get(cID,m_phiID) == 0) ? m_phiBins - 1 : _decoder->get(cID,m_phiID) - 1);
   cellNeighbours.push_back(nID);
   // next: if the current is the last bin (id = m_phiBins - 1) then the next is the first bin (id = 0) else current + 1
   _decoder->set(nID,m_phiID, (_decoder->get(cID,m_phiID) == (m_phiBins - 1)) ? 0 : _decoder->get(cID,m_phiID) + 1);
   cellNeighbours.push_back(nID);

   return cellNeighbours;
}

/// Determine minimum and maximum polar angle of the cell
std::array<double, 2> FCCSWHCalPhiRow_k4geo::cellTheta(const CellID& cID) const {
  std::array<double, 2> cTheta = {M_PI,M_PI};

  // get the cell index
  int idx = _decoder->get(cID, m_rowID);
  // get the layer index
  uint layer = _decoder->get(cID,"layer");

  if(m_radii.empty()) calculateLayerRadii();
  if(m_cellEdges.empty()) return cTheta;

  double zlow = m_cellEdges[layer][idx].first;
  double zhigh = m_cellEdges[layer][idx].second;

  double Rmin = m_radii[layer] - 0.5*m_layerDepth[layer];
  double Rmax = m_radii[layer] + 0.5*m_layerDepth[layer];

  if( theta(cID) < M_PI/2. )
  {
    cTheta[0] = std::atan2(Rmin,zhigh); // theta min
    cTheta[1] = std::atan2(Rmax,zlow); // theta max
  }
  else
  {
    cTheta[0] = std::atan2(Rmax,zhigh); // theta min
    cTheta[1] = std::atan2(Rmin,zlow); // theta max
  }

  return cTheta;
}

}
}
