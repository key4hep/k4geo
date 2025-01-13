#include "detectorSegmentations/FCCSWEndcapTurbine_k4geo.h"
#include "DD4hep/Detector.h"

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
  FCCSWEndcapTurbine_k4geo::FCCSWEndcapTurbine_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) {

    commonSetup();
}

  FCCSWEndcapTurbine_k4geo::FCCSWEndcapTurbine_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) {
  // define type and description

    commonSetup();

}

/// initialize variables, etc (needed for either version of the ctor)
void FCCSWEndcapTurbine_k4geo::commonSetup() {
  _type = "FCCSWEndcapTurbine_k4geo";
  _description = "Turbine-specific segmentation in the global coordinates";

  // register all necessary parameters
  m_offsetRho.resize(3);
  m_gridSizeRho.resize(3);
  m_gridSizeZ.resize(3);
  m_offsetZ.resize(3);
  
  registerParameter("offset_rho1", "Offset in rho1", m_offsetRho[0], 0.);
  registerParameter("offset_rho2", "Offset in rho2", m_offsetRho[1], 0.);
  registerParameter("offset_rho3", "Offset in rho3", m_offsetRho[2], 0.);
  
  registerIdentifier("identifier_rho", "Cell ID identifier for rho", m_rhoID, "rho");
  
  registerParameter("grid_size_rho1", "Grid size in rho1", m_gridSizeRho[0], 0.);
  registerParameter("grid_size_rho2", "Grid size in rho2", m_gridSizeRho[1], 0.);
  registerParameter("grid_size_rho3", "Grid size in rho3", m_gridSizeRho[2], 0.);
  
  registerParameter("grid_size_z1", "Grid size in z1", m_gridSizeZ[0], 0.);
  registerParameter("grid_size_z2", "Grid size in z2", m_gridSizeZ[1], 0.);
  registerParameter("grid_size_z3", "Grid size in z3", m_gridSizeZ[2], 0.);  
  registerParameter("offset_z1", "Offset in z1", m_offsetZ[0], 0.);
  registerParameter("offset_z2", "Offset in z2", m_offsetZ[1], 0.);
  registerParameter("offset_z3", "Offset in z3", m_offsetZ[2], 0.);
  
  registerParameter("offset_theta", "Angular offset in theta", m_offsetTheta, 0., SegmentationParameter::AngleUnit, true);
  registerIdentifier("identifier_z", "Cell ID identifier for z", m_zID, "z");
  registerIdentifier("identifier_side", "Cell ID identifier for side", m_sideID, "side");
  registerIdentifier("identifier_wheel", "Cell ID identifier for wheel", m_wheelID, "wheel");
  registerIdentifier("identifier_module", "Cell ID identifier for module", m_moduleID, "module");
  registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");  
  dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
  
  m_bladeAngle.clear();
  
  try {
    m_bladeAngle.push_back(dd4hepgeo->constant<double>("BladeAngle1"));
  }
  catch(...) {
    std::cout << "BladeAngle1 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  try {
    m_bladeAngle.push_back(dd4hepgeo->constant<double>("BladeAngle2"));
  }
  catch(...) {
    std::cout << "BladeAngle2 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }  
  try {
    m_bladeAngle.push_back(dd4hepgeo->constant<double>("BladeAngle3"));
  }
  catch(...) {
    std::cout << "BladeAngle3 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }

  
  m_nUnitCells.clear();
  try {
    m_nUnitCells.push_back(dd4hepgeo->constant<int>("nUnitCells1"));
  }
  catch(...) {
    std::cout << "nUnitCells1 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  try {
    m_nUnitCells.push_back(dd4hepgeo->constant<int>("nUnitCells2"));
  }
  catch(...) {
    std::cout << "nUnitCells2 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  try {
    m_nUnitCells.push_back(dd4hepgeo->constant<int>("nUnitCells3"));
  }
  catch(...) {
    std::cout << "nUnitCells3 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }

  m_numReadoutRhoLayers.clear();
  try {
    m_numReadoutRhoLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumReadoutRhoLayersWheel1"));
  }
  catch(...) {
    std::cout << "ECalEndcapNumReadoutRhoLayersWheel1 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  try {
    m_numReadoutRhoLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumReadoutRhoLayersWheel2"));
  }
  catch(...) {
    std::cout << "ECalEndcapNumReadoutRhoLayersWheel2 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  try {
    m_numReadoutRhoLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumReadoutRhoLayersWheel3"));
  }
  catch(...) {
    std::cout << "ECalEndcapNumReadoutRhoLayersWheel3 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }

  m_numReadoutZLayers.clear();
  try {
    m_numReadoutZLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumReadoutZLayersWheel1"));
  }
  catch(...) {
    std::cout << "ECalEndcapNumReadoutZLayersWheel1 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  try {
    m_numReadoutZLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumReadoutZLayersWheel2"));
  }
  catch(...) {
    std::cout << "ECalEndcapNumReadoutZLayersWheel2 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  try {
    m_numReadoutZLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumReadoutZLayersWheel3"));
  }
  catch(...) {
    std::cout << "ECalEndcapNumReadoutZLayersWheel3 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }  

 m_numCalibRhoLayers.clear();
  try {
    m_numCalibRhoLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumCalibRhoLayersWheel1"));
  }
  catch(...) {
    std::cout << "ECalEndcapNumCalibRhoLayersWheel1 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  try {
    m_numCalibRhoLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumCalibRhoLayersWheel2"));
  }
  catch(...) {
    std::cout << "ECalEndcapNumCalibRhoLayersWheel2 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  try {
    m_numCalibRhoLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumCalibRhoLayersWheel3"));
  }
  catch(...) {
    std::cout << "ECalEndcapNumCalibRhoLayersWheel3 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }

  m_numCalibZLayers.clear();
  try {
    m_numCalibZLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumCalibZLayersWheel1"));
  }
  catch(...) {
    std::cout << "ECalEndcapNumCalibZLayersWheel1 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  try {
    m_numCalibZLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumCalibZLayersWheel2"));
  }
  catch(...) {
    std::cout << "ECalEndcapNumCalibZLayersWheel2 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  try {
    m_numCalibZLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumCalibZLayersWheel3"));
  }
  catch(...) {
    std::cout << "ECalEndcapNumCalibZLayersWheel3 not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
}
  
/// determine the local position based on the cell ID
Vector3D FCCSWEndcapTurbine_k4geo::position(const CellID& cID) const {

  double rhoVal = rho(cID);
  double zVal = z(cID);
  double phiVal = phi(cID);
  Vector3D pos = PositionRhoZPhi(rhoVal, zVal, phiVal);
  // account for the fact that the -z endcap is mirrored wrt to the +z one
  if (pos.Z < 0.) pos.Y = -pos.Y;

  return pos;
}

/// determine the cell ID based on the position
CellID FCCSWEndcapTurbine_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                          const VolumeID& vID) const {
  CellID cID = vID;
  CellID iWheel = _decoder->get(cID, m_wheelID);
  CellID iLayer = _decoder->get(cID, m_layerID);
  double lRho = rhoFromXYZ(globalPosition);
  int iRho = positionToBin(lRho, m_gridSizeRho[iWheel], m_offsetRho[iWheel]);
  if (iRho < 0) iRho = 0;
  if (iRho >= m_numReadoutRhoLayers[iWheel]) iRho = m_numReadoutRhoLayers[iWheel]-1;

  _decoder->set(cID, m_rhoID, iRho);

  double lZ = TMath::Abs(globalPosition.Z);
  int iZ = positionToBin(lZ, m_gridSizeZ[iWheel], m_offsetZ[iWheel]);
  if (iZ < 0) iZ = 0;
  if (iZ >= m_numReadoutZLayers[iWheel]) iZ = m_numReadoutZLayers[iWheel]-1;
  _decoder->set(cID, m_zID, iZ);

  if (expLayer(iWheel, iRho, iZ) != iLayer) {
    _decoder->set(cID, m_layerID, expLayer(iWheel, iRho, iZ));
  }
  
  return cID;
}

  /// determine rho based on the cell ID
double FCCSWEndcapTurbine_k4geo::rho(const CellID& cID) const {
  CellID rhoValue = _decoder->get(cID, m_rhoID);
  CellID iWheel = _decoder->get(cID, m_wheelID);

  return binToPosition(rhoValue,m_gridSizeRho[iWheel], m_offsetRho[iWheel]);
}  
  
/// determine the azimuthal angle phi based on the cell ID
double FCCSWEndcapTurbine_k4geo::phi(const CellID& cID) const {
  CellID iModule = _decoder->get(cID, m_moduleID);
  CellID iWheel = _decoder->get(cID, m_wheelID);
  
  double phiCent = twopi*(iModule+0.5)/m_nUnitCells[iWheel];
  double rhoLoc = rho(cID);
  double zLoc = TMath::Abs(z(cID))-m_offsetZ[iWheel] - 45/2.;  // hard-code midpoint in z for now

  double zCotBladeAngle = zLoc/TMath::Tan(m_bladeAngle[iWheel]);

  double x = zCotBladeAngle;
  double y = TMath::Sqrt(rhoLoc*rhoLoc - x*x);
  // rotate about z axis by phiCent
  double xprime = x*TMath::Cos(phiCent) +y*TMath::Sin(phiCent);
  double yprime = y*TMath::Cos(phiCent) -x*TMath::Sin(phiCent);
  
  return TMath::ATan2(xprime,yprime);
 }

  
/// determine local x in plane of blade based on the cell ID
double FCCSWEndcapTurbine_k4geo::z(const CellID& cID) const {
  CellID zValue = _decoder->get(cID, m_zID);
  CellID sideValue = _decoder->get(cID, m_sideID);
  CellID iWheel = _decoder->get(cID, m_wheelID);
  return ((long long int)sideValue)*binToPosition(zValue,m_gridSizeZ[iWheel],m_offsetZ[iWheel]);
}

  /// determine expected layer value based on wheel, rho, and z indices
  unsigned FCCSWEndcapTurbine_k4geo::expLayer(unsigned iWheel, unsigned iRho, unsigned iZ) const {
    unsigned layerOffset = 0;
    if (iWheel == 1) {
      layerOffset = m_numCalibZLayers[0]*m_numCalibRhoLayers[0];
    }
    else if (iWheel == 2) {
      layerOffset = m_numCalibZLayers[0]*m_numCalibRhoLayers[0]+layerOffset + m_numCalibZLayers[1]*m_numCalibRhoLayers[1];
    }
    return layerOffset + iZ/(m_numReadoutZLayers[iWheel]/m_numCalibZLayers[iWheel]) + iRho*m_numCalibZLayers[iWheel]/(m_numReadoutRhoLayers[iWheel]/m_numCalibRhoLayers[iWheel]);
  }
}
}

