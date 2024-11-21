#include "detectorSegmentations/FCCSWEndcapTurbine_k4geo.h"
#include "DD4hep/Detector.h"

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
  FCCSWEndcapTurbine_k4geo::FCCSWEndcapTurbine_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) {
  // define type and description
  _type = "FCCSWEndcapTurbine_k4geo";
  _description = "Turbine-specific segmentation in the global coordinates";

  // register all necessary parameters
  //  registerParameter("grid_size_rho", "Grid size in rho", m_gridSizeRho, std::vector<float>());
  // registerParameter("offset_rho", "Offset in rho", m_offsetRho, std::vector<float>());
  // registerIdentifier("identifier_rho", "Cell ID identifier for rho", m_rhoID, "rho"); // might want to have separate concepts for rho and layer...
  registerParameter("grid_size_z", "Grid size in z", m_gridSizeZ, 0.);
  registerParameter("offset_z", "Offset in z", m_offsetZ, 0.);
  registerParameter("offset_theta", "Angular offset in theta", m_offsetTheta, 0., SegmentationParameter::AngleUnit, true);
  registerIdentifier("identifier_z", "Cell ID identifier for z", m_zID, "z");
  registerIdentifier("identifier_side", "Cell ID identifier for side", m_sideID, "side");

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
    
}

  FCCSWEndcapTurbine_k4geo::FCCSWEndcapTurbine_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) {
  // define type and description
  _type = "FCCSWEndcapTurbine_k4geo";
  _description = "Turbine-specific segmentation in the global coordinates";

  // register all necessary parameters
  //  registerParameter("grid_size_rho", "Grid size in rho", m_gridSizeRho, std::vector<float>());
  // registerParameter("offset_rho", "Offset in rho", m_offsetRho, std::vector<float>());
  // registerIdentifier("identifier_rho", "Cell ID identifier for rho", m_rhoID, "rho");
  registerParameter("grid_size_z", "Grid size in z", m_gridSizeZ, 0.);
  registerParameter("offset_z", "Offset in z", m_offsetZ, 0.);
  registerParameter("offset_theta", "Angular offset in theta", m_offsetTheta, 0., SegmentationParameter::AngleUnit, true);
  registerIdentifier("identifier_z", "Cell ID identifier for z", m_zID, "z");
  registerIdentifier("identifier_side", "Cell ID identifier for side", m_sideID, "side");
  registerIdentifier("identifier_wheel", "Cell ID identifier for wheel", m_wheelID, "wheel");
  registerIdentifier("identifier_module", "Cell ID identifier for module", m_moduleID, "module");
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
}

/// determine the local position based on the cell ID
Vector3D FCCSWEndcapTurbine_k4geo::position(const CellID& cID) const {

  //  double rhoVal = rho(cID);  // just the center of the blade
  double zVal = z(cID);
  double phiVal = phi(cID);
  double xVal = x(cID);
  
  //Vector3D pos = PositionRhoZPhi(rhoVal, zVal, phiVal);

  // make a dummy "position" vector that has phi, x, and z
  Vector3D pos(phiVal, 0., 0);
  
  //  // account for the fact that the -z endcap is mirrored wrt to the +z one
  // if (pos.Z < 0.) pos.Y = -pos.Y;
  return pos;
}

/// determine the cell ID based on the position
CellID FCCSWEndcapTurbine_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                          const VolumeID& vID) const {
  CellID cID = vID;

  //  _decoder->set(cID, m_rhoID, positionToBin(lRho, m_gridSizeRho[iWheel], m_offsetRho[iWheel]));
  //_decoder->set(cID, m_rhoID, 0);
  double lZ = TMath::Abs(globalPosition.Z);
  _decoder->set(cID, m_zID, positionToBin(lZ, m_gridSizeZ, m_offsetZ));

  return cID;
}

/// determine the azimuthal angle phi based on the cell ID
double FCCSWEndcapTurbine_k4geo::phi(const CellID& cID) const {

  // just return the phi of the module center
  CellID iModule = _decoder->get(cID, m_moduleID);
  CellID iWheel = _decoder->get(cID, m_wheelID);
  CellID iSide = _decoder->get(cID, m_sideID);
  
  double phiCent = ((long long int)iSide)*twopi*(iModule)/m_nUnitCells[iWheel];

  return phiCent;
  
  //double rhoLoc = rho(cID);
  //double zLoc = TMath::Abs(z(cID))-m_offsetZ - 45;  // hard-code midpoint in z for now

  //double zCotBladeAngle = zLoc/TMath::Tan(m_bladeAngle[iWheel]);
  // double x = zCotBladeAngle;
  //double y = TMath::Sqrt(rhoLoc*rhoLoc - x*x);
  //// rotate about z axis by phiCent
  //double xprime = x*TMath::Cos(phiCent) +y*TMath::Sin(phiCent);
  //double yprime = y*TMath::Cos(phiCent) -x*TMath::Sin(phiCent);
  
  //return TMath::ATan2(xprime,yprime);
}
/// determine the value of the x coordinate relative to the center of the blade
  double FCCSWEndcapTurbine_k4geo::x(const CellID& cID) const {

    CellID iWheel = _decoder->get(cID, m_wheelID);
    double zLoc = TMath::Abs(z(cID))-m_offsetZ - 45;  // hard-code midpoint in z for now

    double zCotBladeAngle = zLoc/TMath::Tan(m_bladeAngle[iWheel]);
    double x = zCotBladeAngle;

    return x;
}

  
/// determine local x in plane of blade based on the cell ID
double FCCSWEndcapTurbine_k4geo::z(const CellID& cID) const {
  CellID zValue = _decoder->get(cID, m_zID);
  CellID sideValue = _decoder->get(cID, m_sideID);
  return ((long long int)sideValue)*binToPosition(zValue,m_gridSizeZ,m_offsetZ);
}

}
}

