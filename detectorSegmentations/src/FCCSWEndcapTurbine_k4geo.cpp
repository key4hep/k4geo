#include "detectorSegmentations/FCCSWEndcapTurbine_k4geo.h"
#include "DD4hep/Detector.h"

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
  FCCSWEndcapTurbine_k4geo::FCCSWEndcapTurbine_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) {
  std::cout << "In string ctor" << std::endl;
  // define type and description
  _type = "FCCSWEndcapTurbine_k4geo";
  _description = "Turbine-specific segmentation in the global coordinates";

  // register all necessary parameters (additional to those registered in GridTheta_k4geo)
  //  registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
  //registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
  //registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
  //registerParameter("rho_bins", "Number of bins in rho", m_rhoBins, 1);
  registerParameter("grid_size_rho", "Grid size in rho", m_gridSizeRho, 0.);
  registerParameter("offset_rho", "Offset in rho", m_offsetRho, 0.);
  registerIdentifier("identifier_rho", "Cell ID identifier for rho", m_rhoID, "rho");
  //registerParameter("z_bins", "Number of bins in z", m_zBins, 1);
  registerParameter("grid_size_z", "Grid size in z", m_gridSizeZ, 0.);
  registerParameter("offset_z", "Offset in z", m_offsetZ, 0.);
  registerIdentifier("identifier_z", "Cell ID identifier for z", m_zID, "z");
  registerIdentifier("identifier_side", "Cell ID identifier for side", m_sideID, "side");

  dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
  try {
    m_bladeAngle = dd4hepgeo->constant<double>("BladeAngle");
    std::cout << "Blade angle is " << m_bladeAngle << std::endl;
  }
  catch(...) {
    std::cout << "Blade angle not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }

}

  FCCSWEndcapTurbine_k4geo::FCCSWEndcapTurbine_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) {
    std::cout << "In BitFieldCoder ctor" << std::endl;
  // define type and description
  _type = "FCCSWEndcapTurbine_k4geo";
  _description = "Turbine-specific segmentation in the global coordinates";

  // register all necessary parameters (additional to those registered in GridTheta_k4geo)
  //  registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
  // registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
  // registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
  //  registerParameter("rho_bins", "Number of bins in rho", m_rhoBins, 1);
  registerParameter("grid_size_rho", "Grid size in rho", m_gridSizeRho, 0.);
  registerParameter("offset_rho", "Offset in rho", m_offsetRho, 0.);
  registerIdentifier("identifier_rho", "Cell ID identifier for rho", m_rhoID, "rho");
  //  registerParameter("z_bins", "Number of bins in z", m_zBins, 1);
  registerParameter("grid_size_z", "Grid size in z", m_gridSizeZ, 0.);
  registerParameter("offset_z", "Offset in z", m_offsetZ, 0.);
  registerIdentifier("identifier_z", "Cell ID identifier for z", m_zID, "z");
  registerIdentifier("identifier_side", "Cell ID identifier for side", m_sideID, "side");
  dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
  try {
    m_bladeAngle = dd4hepgeo->constant<double>("BladeAngle");
    std::cout << "Blade angle is " << m_bladeAngle << std::endl;
  }
  catch(...) {
    std::cout << "Blade angle not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }


}

/// determine the local position based on the cell ID
Vector3D FCCSWEndcapTurbine_k4geo::position(const CellID& cID) const {
  //  Vector3D localPos = PositionRhoZPhi(rho(cID), z(cID), phi(cID));
  //  std::cout << "Local position is " << localPos.x() << " " << localPos.y() << " " << localPos.z() << std::endl;
  //Vector3D dummyPos(z(cID),0,rho(cID));
  // Vector3D dummyPos(x(cID),0,z(cID));  /* this one seemed close when G4 elements were entire blades */
  Vector3D dummyPos(0, 0, 0);
  return dummyPos;
  // return PositionRhoZPhi(rho(cID), z(cID), phi(cID));
  // return PositionRhoZPhi(0, rho(cID), z(cID)/rho(cID));
}

/// determine the cell ID based on the position
CellID FCCSWEndcapTurbine_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                          const VolumeID& vID) const {
  CellID cID = vID;
  double lRho = rhoFromXYZ(globalPosition);
  _decoder->set(cID, m_rhoID, positionToBin(lRho, m_gridSizeRho, m_offsetRho));
  std::cout << "vID = " << std::hex <<vID << std::dec << std::endl;
  std::cout << "lRho = " << lRho << std::endl;
  std::cout << "bin = " << positionToBin(lRho, m_gridSizeRho, m_offsetRho) << std::endl;

  double lZ = TMath::Abs(globalPosition.Z);
  std::cout << "Actual Z position is " << lZ << std::endl;
  std::cout << "Thus Z bin is " << positionToBin(lZ, m_gridSizeZ, m_offsetZ) << std::endl;
  _decoder->set(cID, m_zID, positionToBin(lZ, m_gridSizeZ, m_offsetZ));

  std::cout << "cID = " << std::hex << cID << std::dec << std::endl;  
  std::cout << "Global position (x,y,z) and cellID: " << globalPosition.X << " "  << globalPosition.Y << " " << globalPosition.Z << " " << cID << std::endl;

  //  double lTheta = thetaFromXYZ(globalPosition);
  // double lPhi = phiFromXYZ(globalPosition);
  //_decoder->set(cID, m_thetaID, positionToBin(lTheta, m_gridSizeTheta, m_offsetTheta));
  //_decoder->set(cID, m_phiID, positionToBin(lPhi, 2 * M_PI / (double)m_phiBins, m_offsetPhi));
  return cID;
}

/// determine the azimuthal angle phi based on the current cell ID
//double FCCSWEndcapTurbine_k4geo::phi() const {
//  CellID phiValue = (*_decoder)[m_phiID].value();
//  return binToPosition(phiValue, 2. * M_PI / (double)m_phiBins, m_offsetPhi);
//}

/// determine the azimuthal angle phi based on the cell ID
double FCCSWEndcapTurbine_k4geo::phi(const CellID& cID) const {
  // CellID phiValue = _decoder->get(cID, m_phiID);
  //return binToPosition(phiValue, 2. * M_PI / (double)m_phiBins, m_offsetPhi);
  double rhoLoc = rho(cID)-m_offsetRho;
  double zLoc = TMath::Abs(x(cID))-m_offsetZ;
  double zCotBladeAngle = zLoc/TMath::Tan(m_bladeAngle);

  std::cout << "Calculated phi is " << TMath::ATan(zCotBladeAngle/TMath::Sqrt(rhoLoc*rhoLoc-zCotBladeAngle*zCotBladeAngle)) << std::endl;
    
  //  return 0;

  return TMath::ATan(zCotBladeAngle/TMath::Sqrt(rhoLoc*rhoLoc-zCotBladeAngle*zCotBladeAngle));
}
/// determine the transverse distance from the beamline r based on the cell ID
double FCCSWEndcapTurbine_k4geo::rho(const CellID& cID) const {
  std::cout << "Calculating rho... " << std::endl;
  CellID rhoValue = _decoder->get(cID, m_rhoID);
  std::cout << "Calculated bin is " << rhoValue << std::endl;
  std::cout << "Thus, calculated rho position is " << binToPosition(rhoValue,m_gridSizeRho) << std::endl;
  return binToPosition(rhoValue,m_gridSizeRho);
}
/// determine local x in plane of blade based on the cell ID
double FCCSWEndcapTurbine_k4geo::x(const CellID& cID) const {
  std::cout << "Calculating x for cID " << std::hex << cID << std::dec <<std::endl;
  CellID zValue = _decoder->get(cID, m_zID);
  CellID sideValue = _decoder->get(cID, m_sideID);
  std::cout << "Calculated bin is " << (long long int)zValue << std::endl;
  std::cout << "Side value is " << (long long int)sideValue << std::endl;
  std::cout << "Thus, calculated x position is " << binToPosition(zValue,m_gridSizeZ)/TMath::Cos(m_bladeAngle) << std::endl;
  return -((long long int)sideValue)*binToPosition(zValue,m_gridSizeZ)/TMath::Cos(m_bladeAngle);
}

/// determine local y in plane of blade based on cellID
double FCCSWEndcapTurbine_k4geo::z(const CellID& cID) const {

  std::cout << "Calculating z for cID " << std::hex << cID << std::dec <<std::endl;
  CellID zValue = _decoder->get(cID, m_zID);
  CellID rhoValue = _decoder->get(cID, m_rhoID);

  double zGlob=  binToPosition(zValue,m_gridSizeZ);
  double rho= binToPosition(rhoValue,m_gridSizeRho)+m_offsetRho;
  std::cout << "Rho, zGlob are " << rho << ", " << zGlob << std::endl; 
  std::cout << "Calculated z position is " << -TMath::Sqrt(rho*rho-zGlob*zGlob) - m_offsetRho << std::endl;
  return -TMath::Sqrt(rho*rho-zGlob*zGlob) - m_offsetRho;

}
}
}
