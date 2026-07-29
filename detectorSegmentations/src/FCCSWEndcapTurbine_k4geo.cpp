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
    registerParameter("offset_z", "Offset in z", m_offsetZ, std::vector<double>());
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
      m_numWheels = dd4hepgeo->constant<int>("EMECnWheels");
    } catch (...) {
      std::cout << "EMECnWheels not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_bladeAngle.push_back(dd4hepgeo->constant<double>("EMECBladeAngle1"));
    } catch (...) {
      std::cout << "EMECBladeAngle1 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_bladeAngle.push_back(dd4hepgeo->constant<double>("EMECBladeAngle2"));
    } catch (...) {
      std::cout << "EMECBladeAngle2 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_bladeAngle.push_back(dd4hepgeo->constant<double>("EMECBladeAngle3"));
    } catch (...) {
      std::cout << "EMECBladeAngle3 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }

    m_nUnitCells.clear();
    try {
      m_nUnitCells.push_back(dd4hepgeo->constant<int>("EMECnUnitCells1"));
    } catch (...) {
      std::cout << "EMECnUnitCells1 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_nUnitCells.push_back(dd4hepgeo->constant<int>("EMECnUnitCells2"));
    } catch (...) {
      std::cout << "EMECnUnitCells2 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_nUnitCells.push_back(dd4hepgeo->constant<int>("EMECnUnitCells3"));
    } catch (...) {
      std::cout << "EMECnUnitCells3 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }

    m_numReadoutRhoLayers.clear();
    try {
      m_numReadoutRhoLayers.push_back(dd4hepgeo->constant<int>("EMECNumReadoutRhoLayersWheel1"));
    } catch (...) {
      std::cout << "EMECNumReadoutRhoLayersWheel1 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numReadoutRhoLayers.push_back(dd4hepgeo->constant<int>("EMECNumReadoutRhoLayersWheel2"));
    } catch (...) {
      std::cout << "EMECNumReadoutRhoLayersWheel2 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numReadoutRhoLayers.push_back(dd4hepgeo->constant<int>("EMECNumReadoutRhoLayersWheel3"));
    } catch (...) {
      std::cout << "EMECNumReadoutRhoLayersWheel3 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }

    m_numReadoutZLayers.clear();
    try {
      m_numReadoutZLayers.push_back(dd4hepgeo->constant<int>("EMECNumReadoutZLayersWheel1"));
    } catch (...) {
      std::cout << "EMECNumReadoutZLayersWheel1 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numReadoutZLayers.push_back(dd4hepgeo->constant<int>("EMECNumReadoutZLayersWheel2"));
    } catch (...) {
      std::cout << "EMECNumReadoutZLayersWheel2 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numReadoutZLayers.push_back(dd4hepgeo->constant<int>("EMECNumReadoutZLayersWheel3"));
    } catch (...) {
      std::cout << "EMECNumReadoutZLayersWheel3 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }

    m_numCalibRhoLayers.clear();
    try {
      m_numCalibRhoLayers.push_back(dd4hepgeo->constant<int>("EMECNumCalibRhoLayersWheel1"));
    } catch (...) {
      std::cout << "EMECNumCalibRhoLayersWheel1 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numCalibRhoLayers.push_back(dd4hepgeo->constant<int>("EMECNumCalibRhoLayersWheel2"));
    } catch (...) {
      std::cout << "EMECNumCalibRhoLayersWheel2 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numCalibRhoLayers.push_back(dd4hepgeo->constant<int>("EMECNumCalibRhoLayersWheel3"));
    } catch (...) {
      std::cout << "EMECNumCalibRhoLayersWheel3 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }

    m_numCalibZLayers.clear();
    try {
      m_numCalibZLayers.push_back(dd4hepgeo->constant<int>("EMECNumCalibZLayersWheel1"));
    } catch (...) {
      std::cout << "EMECNumCalibZLayersWheel1 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numCalibZLayers.push_back(dd4hepgeo->constant<int>("EMECNumCalibZLayersWheel2"));
    } catch (...) {
      std::cout << "EMECNumCalibZLayersWheel2 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    try {
      m_numCalibZLayers.push_back(dd4hepgeo->constant<int>("EMECNumCalibZLayersWheel3"));
    } catch (...) {
      std::cout << "EMECNumCalibZLayersWheel3 not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }

    m_rhoIndex = decoder()->index(m_rhoID);
    m_wheelIndex = decoder()->index(m_wheelID);
    m_moduleIndex = decoder()->index(m_moduleID);
    m_zIndex = decoder()->index(m_zID);
    m_sideIndex = decoder()->index(m_sideID);
    m_layerIndex = decoder()->index(m_layerID);

    // precalculate some quantities that will be used to determine positions
    unsigned nWheels = m_bladeAngle.size();
    for (unsigned iWheel = 0; iWheel < nWheels; iWheel++) {
      m_cscBladeAngle.push_back(1. / std::sin(m_bladeAngle[iWheel]));
      m_gangedRhoLayers.push_back(m_numReadoutRhoLayers[iWheel] / m_numCalibRhoLayers[iWheel]);
      m_gangedZLayers.push_back(m_numReadoutZLayers[iWheel] / m_numCalibZLayers[iWheel]);
    }
  }

  /// determine the local position of the readout cell with respect to its
  /// parent calibration cell based on the cell ID
  Vector3D FCCSWEndcapTurbine_k4geo::position(const CellID& cID) const {

    double xVal = 0.;             // this is the direction perpendicular to the blade face
    double zVal = getLocalZ(cID); // coincides with the rho direction in the center
                                  // of the blade
    double yVal = getLocalY(cID);

    Vector3D pos = Position(xVal, yVal, zVal);

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

    double lZ = std::abs(globalPosition.Z);
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

  /// determine the radial position in transverse plane (rho) based on the
  /// cell ID (this is in the global coordinate frame of the detector)
  double FCCSWEndcapTurbine_k4geo::getGlobalRho(const CellID cID) const {
    CellID rhoValue = decoder()->get(cID, m_rhoIndex);
    CellID iWheel = decoder()->get(cID, m_wheelIndex);

    return binToPosition(rhoValue, m_gridSizeRho[iWheel], m_offsetRho[iWheel]) + m_gridSizeRho[iWheel] / 2.;
  }

  /// determine the longitudinal position (z) based on the cell ID, in the
  /// global coordinate frame of the detector
  double FCCSWEndcapTurbine_k4geo::getGlobalZ(const CellID cID) const {
    CellID zValue = decoder()->get(cID, m_zIndex);
    CellID sideValue = decoder()->get(cID, m_sideIndex);
    CellID iWheel = decoder()->get(cID, m_wheelIndex);
    return ((long long int)sideValue) *
           (binToPosition(zValue, m_gridSizeZ[iWheel], m_offsetZ[iWheel]) + m_gridSizeZ[iWheel] / 2.);
  }

  /// determine the phi position based on the cell ID, in the
  /// global coordinate frame of the detector
  double FCCSWEndcapTurbine_k4geo::getGlobalPhi(const CellID cID) const {
    CellID iModule = decoder()->get(cID, m_moduleIndex);
    CellID iWheel = decoder()->get(cID, m_wheelIndex);
    CellID iSide = decoder()->get(cID, m_sideIndex);

    double phiCent = dd4hep::twopi * (iModule + 0.5) / (m_nUnitCells[iWheel] / m_mergedModules[iWheel]);
    double rho = getGlobalRho(cID);

    double zdepth = m_numReadoutZLayers[iWheel] * m_gridSizeZ[iWheel];

    double zFromCent = std::abs(getGlobalZ(cID)) - m_offsetZ[iWheel] - zdepth / 2;

    // calculation position in frame with unit cell at phi = 0
    // line below is equivalent to zFromCent/tan(bladeAngle)
    double y = zFromCent * std::sqrt(m_cscBladeAngle[iWheel] * m_cscBladeAngle[iWheel] - 1.);
    double x = std::sqrt(rho * rho - y * y);
    double locPhi = std::atan2(y, x);

    // now rotate by phi position of the unit cell
    double phi = locPhi + phiCent;
    if (phi > std::numbers::pi)
      phi = phi - dd4hep::twopi;

    if (iSide != 1)
      phi = -phi;

    return phi;
  }

  /// determine the local z position of a readout cell with respect to the
  /// parent calibration cell.  In the local coordinates, z is the direction
  /// pointing upward from the beamline through the center of a turbine blade

  double FCCSWEndcapTurbine_k4geo::getLocalZ(const CellID cID) const {

    /// initialize vector of local Z positions if necessary
    const std::vector<EndcapTurbineWheelLocalZ>* wlzs = m_localZPositions.load();
    if (!wlzs) {
      std::vector<EndcapTurbineWheelLocalZ>* wlzs_new = new std::vector<EndcapTurbineWheelLocalZ>;
      int nWheels = m_numReadoutRhoLayers.size();
      for (int jWheel = 0; jWheel < nWheels; jWheel++) {
        EndcapTurbineWheelLocalZ wlz(m_numReadoutRhoLayers[jWheel], m_numReadoutZLayers[jWheel]);
        wlzs_new->emplace_back(wlz);
        for (int jRho = 0; jRho < m_numReadoutRhoLayers[jWheel]; jRho++) {
          for (int jZ = 0; jZ < m_numReadoutZLayers[jWheel]; jZ++) {
            /// get center position of calibration cell.  The volume is formed from the intersection of a Trd2
            /// and a Tube.  The z position is the center of the Trd2.  So, redo the calculation of the Trd2 size
            /// from k4geo, and find its z center in global coordinates.  Then find the z center of the readout
            /// cell in global coordinates, and subtract the two to get the local z position of the readout cell.
            unsigned jCalibRho = jRho / m_gangedRhoLayers[jWheel];
            float calibRhoMin = m_offsetRho[jWheel] + jCalibRho * m_gangedRhoLayers[jWheel] * m_gridSizeRho[jWheel];
            float calibRhoMax = calibRhoMin + m_gangedRhoLayers[jWheel] * m_gridSizeRho[jWheel];
            float calibZmax = calibRhoMax;
            /// angFactor is equivalent to 1/tan(bladeAngle)^2
            float angFactor = m_cscBladeAngle[jWheel] * m_cscBladeAngle[jWheel] - 1.;
            float calibZmin =
                std::sqrt(calibRhoMin * calibRhoMin - (m_numReadoutZLayers[jWheel] * m_numReadoutZLayers[jWheel] *
                                                       m_gridSizeZ[jWheel] * m_gridSizeZ[jWheel] / 4) *
                                                          angFactor);
            float calibZcent = (calibZmax + calibZmin) / 2.;

            // get center position of readout cell
            float readoutRhoMin = m_offsetRho[jWheel] + jRho * m_gridSizeRho[jWheel];
            float readoutRhoMax = readoutRhoMin + m_gridSizeRho[jWheel];
            float readoutZPos[3];
            readoutZPos[0] = (jZ - m_numReadoutZLayers[jWheel] / 2) * m_gridSizeZ[jWheel] +
                             (1 - m_numReadoutZLayers[jWheel] % 2) * m_gridSizeZ[jWheel] / 2;
            readoutZPos[1] = readoutZPos[0] + m_gridSizeZ[jWheel] / 2;
            readoutZPos[2] = readoutZPos[1] - m_gridSizeZ[jWheel];

            std::vector<float> readoutZmaxArr, readoutZminArr;
            for (int jPos = 0; jPos < 3; jPos++) {
              if (readoutRhoMax > std::abs(readoutZPos[jPos])) {
                readoutZmaxArr.push_back(
                    std::sqrt(readoutRhoMax * readoutRhoMax - readoutZPos[jPos] * readoutZPos[jPos]));
              } else {
                readoutZmaxArr.push_back(0);
              }

              if (readoutRhoMin > std::abs(readoutZPos[jPos])) {
                readoutZminArr.push_back(
                    std::sqrt(readoutRhoMin * readoutRhoMin - readoutZPos[jPos] * readoutZPos[jPos]));
              } else {
                readoutZminArr.push_back(1000000);
              }
            }
            float readoutZmax = *(std::ranges::max_element(readoutZmaxArr));
            float readoutZmin = *(std::ranges::min_element(readoutZminArr));
            float readoutZcent = (readoutZmax + readoutZmin) / 2.;

            wlzs_new->at(jWheel).setLocalZ(jRho, jZ, readoutZcent - calibZcent);
          }
        }
      }
      if (m_localZPositions.compare_exchange_strong(wlzs, wlzs_new)) {
        wlzs = wlzs_new;
      } else {
        delete wlzs_new;
      }
    }

    CellID iWheel = decoder()->get(cID, m_wheelIndex);
    CellID iRho = decoder()->get(cID, m_rhoIndex);
    CellID iZ = decoder()->get(cID, m_rhoIndex);

    return wlzs->at(iWheel).getLocalZ(iRho, iZ);
  }

  /// determine the local y position of a readout cell with respect to the
  /// parent calibration cell.  In the local coordinates, y is the direction
  /// pointing across the turbine blade, from one straight edge to the other

  double FCCSWEndcapTurbine_k4geo::getLocalY(const CellID cID) const {

    CellID iWheel = decoder()->get(cID, m_wheelIndex);
    CellID iZ = decoder()->get(cID, m_zIndex);

    return ((iZ % m_gangedZLayers[iWheel] - m_gangedZLayers[iWheel] / 2.) * m_gridSizeZ[iWheel] +
            (1 - m_gangedZLayers[iWheel] % 2) * m_gridSizeZ[iWheel] / 2) *
           m_cscBladeAngle[iWheel];
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
