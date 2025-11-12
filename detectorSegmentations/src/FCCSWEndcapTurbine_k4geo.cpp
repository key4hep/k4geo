#include "detectorSegmentations/FCCSWEndcapTurbine_k4geo.h"
#include "DD4hep/Detector.h"

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  FCCSWEndcapTurbine_k4geo::FCCSWEndcapTurbine_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) {
    commonSetup();
  }

  FCCSWEndcapTurbine_k4geo::FCCSWEndcapTurbine_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) {
    commonSetup();
  }

  /// initialize variables, etc (needed for either version of the ctor)
  void FCCSWEndcapTurbine_k4geo::commonSetup() {
    // define type and description
    _type = "FCCSWEndcapTurbine_k4geo";
    _description = "Turbine-specific segmentation in the global coordinates";

    // register all necessary parameters
    registerParameter("offset_rho", "Offset in rho", m_offsetRho, std::vector<double>());

    registerIdentifier("identifier_rho", "Cell ID identifier for rho", m_rhoID, "rho");

    registerParameter("grid_size_rho", "Grid size in rho", m_gridSizeRho, std::vector<double>());
    registerParameter("grid_size_z", "Grid size in z", m_gridSizeZ, std::vector<double>());
    registerParameter("offset_z", "Offset in z1", m_offsetZ, std::vector<double>());
    registerParameter("offset_theta", "Angular offset in theta", m_offsetTheta, 0., SegmentationParameter::AngleUnit,
                      true);
    registerParameter("mergedModules", "Number of merged modules per wheel", m_mergedModules, std::vector<int>());
    registerIdentifier("identifier_z", "Cell ID identifier for z", m_zID, "z");
    registerIdentifier("identifier_side", "Cell ID identifier for side", m_sideID, "side");
    registerIdentifier("identifier_wheel", "Cell ID identifier for wheel", m_wheelID, "wheel");
    registerIdentifier("identifier_module", "Cell ID identifier for module", m_moduleID, "module");
    registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");
    dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());

    m_bladeAngle.clear();

    try {
      m_bladeAngle.push_back(dd4hepgeo->constant<double>("BladeAngle1"));
    } catch (...) {
      std::cout << "BladeAngle1 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_bladeAngle.push_back(dd4hepgeo->constant<double>("BladeAngle2"));
    } catch (...) {
      std::cout << "BladeAngle2 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_bladeAngle.push_back(dd4hepgeo->constant<double>("BladeAngle3"));
    } catch (...) {
      std::cout << "BladeAngle3 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }

    m_nUnitCells.clear();
    try {
      m_nUnitCells.push_back(dd4hepgeo->constant<int>("nUnitCells1"));
    } catch (...) {
      std::cout << "nUnitCells1 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_nUnitCells.push_back(dd4hepgeo->constant<int>("nUnitCells2"));
    } catch (...) {
      std::cout << "nUnitCells2 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_nUnitCells.push_back(dd4hepgeo->constant<int>("nUnitCells3"));
    } catch (...) {
      std::cout << "nUnitCells3 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }

    m_numReadoutRhoLayers.clear();
    try {
      m_numReadoutRhoLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumReadoutRhoLayersWheel1"));
    } catch (...) {
      std::cout << "ECalEndcapNumReadoutRhoLayersWheel1 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numReadoutRhoLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumReadoutRhoLayersWheel2"));
    } catch (...) {
      std::cout << "ECalEndcapNumReadoutRhoLayersWheel2 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numReadoutRhoLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumReadoutRhoLayersWheel3"));
    } catch (...) {
      std::cout << "ECalEndcapNumReadoutRhoLayersWheel3 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }

    m_numReadoutZLayers.clear();
    try {
      m_numReadoutZLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumReadoutZLayersWheel1"));
    } catch (...) {
      std::cout << "ECalEndcapNumReadoutZLayersWheel1 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numReadoutZLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumReadoutZLayersWheel2"));
    } catch (...) {
      std::cout << "ECalEndcapNumReadoutZLayersWheel2 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numReadoutZLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumReadoutZLayersWheel3"));
    } catch (...) {
      std::cout << "ECalEndcapNumReadoutZLayersWheel3 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }

    m_numCalibRhoLayers.clear();
    try {
      m_numCalibRhoLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumCalibRhoLayersWheel1"));
    } catch (...) {
      std::cout << "ECalEndcapNumCalibRhoLayersWheel1 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numCalibRhoLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumCalibRhoLayersWheel2"));
    } catch (...) {
      std::cout << "ECalEndcapNumCalibRhoLayersWheel2 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numCalibRhoLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumCalibRhoLayersWheel3"));
    } catch (...) {
      std::cout << "ECalEndcapNumCalibRhoLayersWheel3 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }

    m_numCalibZLayers.clear();
    try {
      m_numCalibZLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumCalibZLayersWheel1"));
    } catch (...) {
      std::cout << "ECalEndcapNumCalibZLayersWheel1 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numCalibZLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumCalibZLayersWheel2"));
    } catch (...) {
      std::cout << "ECalEndcapNumCalibZLayersWheel2 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numCalibZLayers.push_back(dd4hepgeo->constant<int>("ECalEndcapNumCalibZLayersWheel3"));
    } catch (...) {
      std::cout << "ECalEndcapNumCalibZLayersWheel3 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }

    m_rhoIndex = decoder()->index(m_rhoID);
    m_wheelIndex = decoder()->index(m_wheelID);
    m_moduleIndex = decoder()->index(m_moduleID);
    m_zIndex = decoder()->index(m_zID);
    m_sideIndex = decoder()->index(m_sideID);
    m_layerIndex = decoder()->index(m_layerID);
  }

  /// determine the local position based on the cell ID
  Vector3D FCCSWEndcapTurbine_k4geo::position(const CellID& cID) const {

    double rhoVal = rho(cID);
    double zVal = z(cID);
    double phiVal = phi(cID);
    Vector3D pos = PositionRhoZPhi(rhoVal, zVal, phiVal);
    // account for the fact that the -z endcap is mirrored wrt to the +z one
    if (pos.Z < 0.)
      pos.Y = -pos.Y;

    return pos;
  }

  /// determine the cell ID based on the position
  CellID FCCSWEndcapTurbine_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                                          const VolumeID& vID) const {
    CellID cID = vID;
    CellID iWheel = decoder()->get(cID, m_wheelIndex);
    CellID iLayer = decoder()->get(cID, m_layerIndex);
    CellID iModule = decoder()->get(cID, m_moduleIndex);

    double lRho = rhoFromXYZ(globalPosition);
    int iRho = positionToBin(lRho, m_gridSizeRho[iWheel], m_offsetRho[iWheel] + m_gridSizeRho[iWheel] / 2.);
    if (iRho < 0) {
      iRho = 0;
    }
    if (iRho >= m_numReadoutRhoLayers[iWheel]) {
      iRho = m_numReadoutRhoLayers[iWheel] - 1;
    }
    decoder()->set(cID, m_rhoIndex, iRho);

    double lZ = TMath::Abs(globalPosition.Z);
    int iZ = positionToBin(lZ, m_gridSizeZ[iWheel], m_offsetZ[iWheel] + m_gridSizeZ[iWheel] / 2.);
    if (iZ < 0) {
      iZ = 0;
    }
    if (iZ >= m_numReadoutZLayers[iWheel]) {
      iZ = m_numReadoutZLayers[iWheel] - 1;
    }
    decoder()->set(cID, m_zIndex, iZ);

    if (expLayer(iWheel, iRho, iZ) != iLayer) {
      decoder()->set(cID, m_layerIndex, expLayer(iWheel, iRho, iZ));
    }

    // adjust module number to account for merging
    iModule = iModule / m_mergedModules[iWheel];

    decoder()->set(cID, m_moduleIndex, iModule);

    return cID;
  }

  /// determine rho based on the cell ID
  double FCCSWEndcapTurbine_k4geo::rho(const CellID cID) const {
    CellID rhoValue = decoder()->get(cID, m_rhoIndex);
    CellID iWheel = decoder()->get(cID, m_wheelIndex);

    return binToPosition(rhoValue, m_gridSizeRho[iWheel], m_offsetRho[iWheel]) + m_gridSizeRho[iWheel] / 2.;
  }

  /// determine the azimuthal angle phi based on the cell ID
  double FCCSWEndcapTurbine_k4geo::phi(const CellID cID) const {
    CellID iModule = decoder()->get(cID, m_moduleIndex);
    CellID iWheel = decoder()->get(cID, m_wheelIndex);

    double phiCent = twopi * (iModule + 0.5) / (m_nUnitCells[iWheel] / m_mergedModules[iWheel]);
    double rhoLoc = rho(cID);

    double zdepth = m_numReadoutZLayers[iWheel] * m_gridSizeZ[iWheel];

    double zLoc = TMath::Abs(z(cID)) - m_offsetZ[iWheel] - zdepth / 2;
    double x = zLoc / TMath::Tan(m_bladeAngle[iWheel]);
    double y = TMath::Sqrt(rhoLoc * rhoLoc - x * x);
    // rotate about z axis by phiCent
    double xprime = x * TMath::Cos(phiCent) + y * TMath::Sin(phiCent);
    double yprime = y * TMath::Cos(phiCent) - x * TMath::Sin(phiCent);

    return TMath::ATan2(xprime, yprime);
  }

  /// determine the absolute z position based on the CellID
  double FCCSWEndcapTurbine_k4geo::z(const CellID cID) const {
    CellID zValue = decoder()->get(cID, m_zIndex);
    CellID sideValue = decoder()->get(cID, m_sideIndex);
    CellID iWheel = decoder()->get(cID, m_wheelIndex);
    return ((long long int)sideValue) *
           (binToPosition(zValue, m_gridSizeZ[iWheel], m_offsetZ[iWheel]) + m_gridSizeZ[iWheel] / 2.);
  }

  /// determine expected layer value based on wheel, rho, and z indices
  unsigned FCCSWEndcapTurbine_k4geo::expLayer(unsigned iWheel, unsigned iRho, unsigned iZ) const {
    unsigned layerOffset = 0;
    if (iWheel == 1) {
      layerOffset = m_numCalibZLayers[0] * m_numCalibRhoLayers[0];
    } else if (iWheel == 2) {
      layerOffset =
          m_numCalibZLayers[0] * m_numCalibRhoLayers[0] + layerOffset + m_numCalibZLayers[1] * m_numCalibRhoLayers[1];
    }
    return layerOffset + iZ / (m_numReadoutZLayers[iWheel] / m_numCalibZLayers[iWheel]) +
           m_numCalibZLayers[iWheel] * (iRho / (m_numReadoutRhoLayers[iWheel] / m_numCalibRhoLayers[iWheel]));
  }
} // namespace DDSegmentation
} // namespace dd4hep
