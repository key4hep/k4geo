#include "detectorSegmentations/FCCSWGridModuleThetaMerged.h"

#include <iostream>
#include "DD4hep/Detector.h"

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
FCCSWGridModuleThetaMerged::FCCSWGridModuleThetaMerged(const std::string& cellEncoding) : GridTheta(cellEncoding) {
  // define type and description
  _type = "FCCSWGridModuleThetaMerged";
  _description = "Module-theta segmentation with per-layer merging along theta and/or module";

  // register all necessary parameters (additional to those registered in GridTheta)
  registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");
  registerIdentifier("identifier_module", "Cell ID identifier for readout module", m_moduleID, "module");
  registerParameter("mergedCells_Theta", "Numbers of merged cells in theta per layer", m_mergedCellsTheta, std::vector<int>());
  registerParameter("mergedModules", "Numbers of merged modules per layer", m_mergedModules, std::vector<int>());
  GetNModulesFromGeom();
  GetNLayersFromGeom();
}

FCCSWGridModuleThetaMerged::FCCSWGridModuleThetaMerged(const BitFieldCoder* decoder) : GridTheta(decoder) {
  // define type and description
  _type = "FCCSWGridModuleThetaMerged";
  _description = "Module-theta segmentation with per-layer merging along theta and/or module";

  // register all necessary parameters (additional to those registered in GridTheta)
  registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");
  registerIdentifier("identifier_module", "Cell ID identifier for module", m_moduleID, "module");
  registerParameter("mergedCells_Theta", "Numbers of merged cells in theta per layer", m_mergedCellsTheta, std::vector<int>());
  registerParameter("mergedModules", "Numbers of merged modules per layer", m_mergedModules, std::vector<int>());
  GetNModulesFromGeom();
  GetNLayersFromGeom();
}

void FCCSWGridModuleThetaMerged::GetNModulesFromGeom() {
  dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
  try {
    m_nModules = dd4hepgeo->constant<int>("ECalBarrelNumPlanes");
  }
  catch(...) {
    std::cout << "Number of modules not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  std::cout << "Number of modules read from detector metadata and used in readout class: " << m_nModules << std::endl;
}

void FCCSWGridModuleThetaMerged::GetNLayersFromGeom() {
  dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
  try {
    m_nLayers = dd4hepgeo->constant<int>("ECalBarrelNumLayers");
  }
  catch(...) {
    std::cout << "Number of layers not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  std::cout << "Number of layers read from detector metadata and used in readout class: " << m_nLayers << std::endl;
}

/// determine the local position based on the cell ID
Vector3D FCCSWGridModuleThetaMerged::position(const CellID& cID) const {

  // debug
  // std::cout << "cellID: " << cID << std::endl;

  // cannot return for R=0 otherwise will lose phi info, return for R=1
  return positionFromRThetaPhi(1.0, theta(cID), phi(cID));
}

/// determine the cell ID based on the global position
CellID FCCSWGridModuleThetaMerged::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                          const VolumeID& vID) const {
  CellID cID = vID;

  // retrieve layer (since merging depends on layer)
  int layer = _decoder->get(vID, m_layerID);

  // retrieve theta
  double lTheta = thetaFromXYZ(globalPosition);

  // calculate theta bin with original segmentation
  int thetaBin = positionToBin(lTheta, m_gridSizeTheta, m_offsetTheta);

  // adjust theta bin if cells are merged along theta in this layer
  // assume that m_mergedCellsTheta[layer]>=1
  thetaBin -= (thetaBin % m_mergedCellsTheta[layer]);

  // set theta field of cellID
  _decoder->set(cID, m_thetaID, thetaBin);

  // retrieve module number
  int module = _decoder->get(vID, m_moduleID);

  // adjust module number if modules are merged in this layer
  // assume that m_mergedModules[layer]>=1
  module -= (module % m_mergedModules[layer]);

  // set module field of cellID
  _decoder->set(cID, m_moduleID, module);

  return cID;
}

/// determine the azimuth based on the cell ID
/// the value returned is the relative shift in phi
/// with respect to the first module in the group of
/// merged ones - which will be then added on top of
/// the phi of the volume containing the first cell
/// by the positioning tool
double FCCSWGridModuleThetaMerged::phi(const CellID& cID) const {

  // retrieve layer
  int layer = _decoder->get(cID, m_layerID);

  // calculate phi offset due to merging
  // assume that m_mergedModules[layer]>=1
  double phi = (m_mergedModules[layer]-1) * M_PI / m_nModules ;

  // debug
  // std::cout << "layer: " << layer << std::endl;
  // std::cout << "merged modules: " << m_mergedModules[layer] << std::endl;
  // std::cout << "phi: " << phi << std::endl;

  return phi;
}

/// determine the polar angle based on the cell ID and the
/// number of merged theta cells
double FCCSWGridModuleThetaMerged::theta(const CellID& cID) const {

  // retrieve layer
  int layer = _decoder->get(cID, m_layerID);

  // retrieve theta bin from cellID and determine theta position
  CellID thetaValue = _decoder->get(cID, m_thetaID);
  double _theta = binToPosition(thetaValue, m_gridSizeTheta, m_offsetTheta);

  // adjust return value if cells are merged along theta in this layer
  // shift by (N-1)*half theta grid size
  // assume that m_mergedCellsTheta[layer]>=1
  _theta += (m_mergedCellsTheta[layer]-1) * m_gridSizeTheta / 2.0 ;

  // debug
  // std::cout << "layer: " << layer << std::endl;
  // std::cout << "theta bin: " << thetaValue << std::endl;
  // std::cout << "gridSizeTheta, offsetTheta: " << m_gridSizeTheta << " , " << m_offsetTheta << std::endl;
  // std::cout << "merged cells: " << m_mergedCellsTheta[layer] << std::endl;
  // std::cout << "theta: " << _theta << std::endl;

  return _theta;
}

}
}
