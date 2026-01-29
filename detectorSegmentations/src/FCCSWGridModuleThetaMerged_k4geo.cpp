#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

#include "DD4hep/Detector.h"
#include "DD4hep/Printout.h"
#include "DD4hep/VolumeManager.h"
#include <iostream>

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  FCCSWGridModuleThetaMerged_k4geo::FCCSWGridModuleThetaMerged_k4geo(const std::string& cellEncoding)
      : GridTheta_k4geo(cellEncoding) {
    commonSetup();
  }

  FCCSWGridModuleThetaMerged_k4geo::FCCSWGridModuleThetaMerged_k4geo(const BitFieldCoder* decoder)
      : GridTheta_k4geo(decoder) {
    commonSetup();
  }

  /// Initialization common to all ctors.
  void FCCSWGridModuleThetaMerged_k4geo::commonSetup() {
    // define type and description
    _type = "FCCSWGridModuleThetaMerged_k4geo";
    _description = "Module-theta segmentation with per-layer merging along theta and/or module";

    // register all necessary parameters (additional to those registered in GridTheta_k4geo)
    registerIdentifier("identifier_cryo", "Cell ID identifier for the cryostat", m_cryoID, "cryo");
    registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");
    registerIdentifier("identifier_module", "Cell ID identifier for readout module", m_moduleID, "module");
    registerParameter("mergedCells_Theta", "Numbers of merged cells in theta per layer", m_mergedCellsTheta,
                      std::vector<int>());
    registerParameter("mergedModules", "Numbers of merged modules per layer", m_mergedModules, std::vector<int>());
    GetNModulesFromGeom();
    GetNLayersFromGeom();

    m_cryoIndex = decoder()->index(m_cryoID);
    m_layerIndex = decoder()->index(m_layerID);
    m_thetaIndex = decoder()->index(fieldNameTheta());
    m_moduleIndex = decoder()->index(m_moduleID);
  }

  FCCSWGridModuleThetaMerged_k4geo::~FCCSWGridModuleThetaMerged_k4geo() { delete m_layerInfo; }

  void FCCSWGridModuleThetaMerged_k4geo::GetNModulesFromGeom() {
    dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
    try {
      m_nModules = dd4hepgeo->constant<int>("ECalBarrelNumPlanes");
    } catch (...) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWGridModuleThetaMerged_k4geo",
                       "Number of modules not found in detector metadata, exiting...");
      exit(1);
    }
    dd4hep::printout(dd4hep::INFO, "FCCSWGridModuleThetaMerged_k4geo",
                     "Number of modules read from detector metadata and used in readout class = %d", m_nModules);
  }

  void FCCSWGridModuleThetaMerged_k4geo::GetNLayersFromGeom() {
    dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
    try {
      m_nLayers = dd4hepgeo->constant<int>("ECalBarrelNumLayers");
    } catch (...) {
      dd4hep::printout(dd4hep::ERROR, "FCCSWGridModuleThetaMerged_k4geo",
                       "Number of layers not found in detector metadata, exiting...");
      exit(1);
    }
    dd4hep::printout(dd4hep::INFO, "FCCSWGridModuleThetaMerged_k4geo",
                     "Number of layers read from detector metadata and used in readout class = %d", m_nLayers);
  }

  /// Tabulate the cylindrical radii of all layers, as well as the
  /// local x and z components needed for the proper phi offset.
  std::vector<FCCSWGridModuleThetaMerged_k4geo::LayerInfo>
  FCCSWGridModuleThetaMerged_k4geo::initLayerInfo(const CellID& cID) const {

    dd4hep::printout(dd4hep::INFO, "FCCSWGridModuleThetaMerged_k4geo", "Precalculating position info of radial layers");
    dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
    VolumeManager vman = VolumeManager::getVolumeManager(*dd4hepgeo);

    std::vector<LayerInfo> out;
    out.reserve(m_nLayers);
    VolumeID vID = cID;
    decoder()->set(vID, m_thetaIndex, 0);
    for (int l = 0; l < m_nLayers; l++) {
      dd4hep::printout(dd4hep::INFO, "FCCSWGridModuleThetaMerged_k4geo", "Layer = %d", l);
      // Look up a volume in layer l in the volume manager, and find its radius
      // by transforming the origin in the local coordinate system to global
      // coordinates.
      decoder()->set(vID, m_layerIndex, l);
      VolumeManagerContext* vc = vman.lookupContext(vID);
      Position wpos = vc->localToWorld({0, 0, 0});
      double rho = wpos.Rho();

      // If different modules are ganged together, we want to put hits
      // in the center of all the ganged modules.  This represents
      // a rotation in phi (in global coordinates).  Convert this rotation
      // to displacement in local coordinates.  This will be in the x-z plane,
      // and will be constant for all modules in a layer, but will be different
      // for different layers (even with identical ganging).
      double xloc = 0;
      double zloc = 0;
      double phioff = phi(vID);
      if (phioff > 0) {
        // We need to apply a phi offset.  Calculate it by rotating
        // the global position in phi and converting back to local
        // coordinates.
        //
        // For the Allegro calorimeter, this can also be calculated as:
        //
        //   d = rho * 2 * sin(phioff/2)
        //   beta = pi/2 - phioff/2 - asin(sin(InclinationAngle)*EMBarrel_rmin/rho)
        //   xloc = - d * sin(beta)
        //   yloc =   d * cos(beta)
        //
        // but we prefer to do the calculation via explicit rotations because
        // it's easier to see that that was is correct, and it also avoids
        // the explicit dependencies on the geometry parameters.
        Position wpos2 = RotateZ(wpos, phioff);
        Position lpos2 = vc->worldToLocal(wpos2);
        xloc = lpos2.X();
        zloc = lpos2.Z();
      }
      dd4hep::printout(dd4hep::INFO, "FCCSWGridModuleThetaMerged_k4geo", "rho, xloc, zloc = %lf %lf %lf (cm)",
                       rho / dd4hep::cm, xloc / dd4hep::cm, zloc / dd4hep::cm);
      out.emplace_back(rho, xloc, zloc);
    }
    return out;
  }

  /// determine the local position based on the cell ID
  Vector3D FCCSWGridModuleThetaMerged_k4geo::position(const CellID& cID) const {

    // Get the vector of layer info.  If it hasn't been made yet,
    // calculate it now.
    // skip cells in the cryostat (if the cryostat is active, there are cells in there,
    // but it has no longitudinal layers)
    if (decoder()->get(cID, m_cryoIndex) != 0)
      return Vector3D(0, 0, 0);
    const std::vector<LayerInfo>* liv = m_layerInfo.load();
    if (!liv) {
      auto liv_new = new std::vector<LayerInfo>(initLayerInfo(cID));
      if (m_layerInfo.compare_exchange_strong(liv, liv_new)) {
        liv = liv_new;
      } else {
        delete liv_new;
      }
    }

    VolumeID vID = cID;
    decoder()->set(vID, m_thetaIndex, 0);
    int layer = this->layer(vID);

    // debug (run ddsim with --printLevel 1 option to see these messages)
    // dd4hep::printout(dd4hep::VERBOSE, "FCCSWGridModuleThetaMerged_k4geo", "cellID = %lu", cID);

    // Calculate the position in local coordinates.
    // The volume here has the cross-section of a cell in the x-z plane;
    // it extends the length of the calorimeter along the y-axis.
    // We set the y-coordinate based on the theta bin, and x and z
    // based on the phi offset required for this layer.
    const LayerInfo& li = liv->at(layer);
    return Vector3D(li.xloc, -li.rho / tan(theta(cID)), li.zloc);
  }

  /// determine the cell ID based on the global position
  CellID FCCSWGridModuleThetaMerged_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                                                  const VolumeID& vID) const {
    CellID cID = vID;

    // retrieve layer (since merging depends on layer)
    int layer = this->layer(vID);

    // retrieve theta
    double lTheta = thetaFromXYZ(globalPosition);

    // calculate theta bin with original segmentation
    int thetaBin = positionToBin(lTheta, gridSizeTheta(), offsetTheta());

    // adjust theta bin if cells are merged along theta in this layer
    // assume that m_mergedCellsTheta[layer]>=1
    thetaBin -= (thetaBin % m_mergedCellsTheta[layer]);

    // set theta field of cellID
    decoder()->set(cID, m_thetaIndex, thetaBin);

    // retrieve module number
    int module = decoder()->get(vID, m_moduleIndex);

    // adjust module number if modules are merged in this layer
    // assume that m_mergedModules[layer]>=1
    module -= (module % m_mergedModules[layer]);

    // set module field of cellID
    decoder()->set(cID, m_moduleIndex, module);

    return cID;
  }

  /// determine the azimuth based on the cell ID
  /// the value returned is the relative shift in phi
  /// with respect to the first module in the group of
  /// merged ones - which will be then added on top of
  /// the phi of the volume containing the first cell
  /// by the positioning tool
  double FCCSWGridModuleThetaMerged_k4geo::phi(const CellID cID) const {

    // retrieve layer
    int layer = this->layer(cID);

    // calculate phi offset due to merging
    // assume that m_mergedModules[layer]>=1
    double phi = (m_mergedModules[layer] - 1) * M_PI / m_nModules;

    // debug
    // dd4hep::printout(dd4hep::VERBOSE, "FCCSWGridModuleThetaMerged_k4geo", "layer = %d, merged modules = %d, phi =
    // %lf",
    //                  layer, m_mergedModules[layer], phi);

    return phi;
  }

  /// determine the polar angle based on the cell ID and the
  /// number of merged theta cells
  double FCCSWGridModuleThetaMerged_k4geo::theta(const CellID cID) const {

    // retrieve layer
    int layer = this->layer(cID);

    // retrieve theta bin from cellID and determine theta position
    CellID thetaValue = decoder()->get(cID, m_thetaIndex);
    double _theta = binToPosition(thetaValue, gridSizeTheta(), offsetTheta());

    // adjust return value if cells are merged along theta in this layer
    // shift by (N-1)*half theta grid size
    // assume that m_mergedCellsTheta[layer]>=1
    _theta += (m_mergedCellsTheta[layer] - 1) * gridSizeTheta() / 2.0;

    // debug
    // dd4hep::printout(dd4hep::VERBOSE, "FCCSWGridModuleThetaMerged_k4geo",
    //                  "layer = %d, theta bin = %d, merged cells in theta = %d, theta = %lf", layer, thetaValue,
    //                  m_mergedCellsTheta[layer], _theta);
    // std::cout << "gridSizeTheta, offsetTheta: " << m_gridSizeTheta << " , " << m_offsetTheta << std::endl;
    return _theta;
  }

  /// Extract the layer index fom a cell ID.
  int FCCSWGridModuleThetaMerged_k4geo::layer(const CellID cID) const { return decoder()->get(cID, m_layerIndex); }

  /// Determine the volume ID from the full cell ID by removing all local fields
  VolumeID FCCSWGridModuleThetaMerged_k4geo::volumeID(const CellID& cID) const {
    VolumeID vID = cID;
    decoder()->set(vID, m_thetaIndex, 0);
    return vID;
  }

} // namespace DDSegmentation
} // namespace dd4hep
