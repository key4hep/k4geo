#include "detectorSegmentations/GridDRcalo_k4geo.h"
#include "DD4hep/Printout.h"

#include <climits>
#include <cmath>
#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  GridDRcalo_k4geo::GridDRcalo_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) { commonSetup(); }

  GridDRcalo_k4geo::GridDRcalo_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) { commonSetup(); }

  GridDRcalo_k4geo::~GridDRcalo_k4geo() {
    delete m_paramBarrel;
    delete m_paramEndcap;
  }

  /// Initialization common to all ctors.
  void GridDRcalo_k4geo::commonSetup() {
    // define type and description
    _type = "GridDRcalo_k4geo";
    _description = "DRcalo segmentation based on the tower / (Cerenkov or Scintillation) fiber / SiPM hierarchy";

    // register all necessary parameters
    registerIdentifier("identifier_assembly", "Cell ID identifier for assembly", m_assemblyID, "assembly");
    registerIdentifier("identifier_eta", "Cell ID identifier for numEta", m_numEtaID, "eta");
    registerIdentifier("identifier_phi", "Cell ID identifier for numPhi", m_numPhiID, "phi");
    registerIdentifier("identifier_x", "Cell ID identifier for x", m_xID, "x");
    registerIdentifier("identifier_y", "Cell ID identifier for y", m_yID, "y");
    registerIdentifier("identifier_IsCerenkov", "Cell ID identifier for IsCerenkov", m_isCerenkovID, "c");
    registerIdentifier("identifier_module", "Cell ID identifier for module", m_moduleID, "module");

    m_paramBarrel = new DRparamBarrel_k4geo();
    m_paramEndcap = new DRparamEndcap_k4geo();

    m_systemIndex = decoder()->index("system");
    m_numEtaIndex = decoder()->index(m_numEtaID);
    m_numPhiIndex = decoder()->index(m_numPhiID);
    m_moduleIndex = decoder()->index(m_moduleID);
    m_xIndex = decoder()->index(m_xID);
    m_yIndex = decoder()->index(m_yID);
    m_isCerenkovIndex = decoder()->index(m_isCerenkovID);
    m_assemblyIndex = decoder()->index(m_assemblyID);
  }

  // front end position (default calo hit position)
  Vector3D GridDRcalo_k4geo::position(const CellID& cID) const {
    int noEta = numEta(cID); // always positive since we replicated the RHS
    int noPhi = numPhi(cID);
    bool isRHS = IsRHS(cID);

    // first get a vector to the center of the wafer (per tower)
    DRparamBase_k4geo* paramBase = setParamBase(noEta);
    auto waferPos = paramBase->GetSipmLayerPos(noPhi);

    dd4hep::Position localPos = dd4hep::Position(0., 0., 0.);

    // get local coordinate
    // if IsSiPM is false, the local position is (0, 0, 0)
    // and the function returns the tower's position
    if (IsSiPM(cID))
      localPos = dd4hep::Position(localPosition(cID));

    // rotate the local coordinate on par to the tower
    auto rot = paramBase->GetRotationZYX(noPhi);
    auto translation = rot * localPos;

    // we want the cell position to be the front end of the tower
    // retrieve the fiber length
    int row = y(cID);
    int col = x(cID);
    // retrieving fiber lengths - if the fiber is not within the wing region of the tower
    // then the method will return the full length
    double fiberLen = paramBase->GetShortFibers(noEta).retrieveFiberLength(row, col);
    // and vector to the front end
    auto diff = paramBase->GetTowerPos(noPhi) - waferPos;
    auto unitVec = diff / std::sqrt(diff.Mag2());

    // total vector is sum of the vector to the wafer center + rotated local coordinate
    // + the translation to the front end of the tower
    auto total = translation + waferPos + unitVec * fiberLen;

    // if LHS rotate by 180 deg w.r.t. X axis (on par to the DRconstructor)
    if (!isRHS)
      total = dd4hep::RotationX(M_PI) * total;

    dd4hep::printout(dd4hep::VERBOSE, "GridDRcalo_k4geo::position",
                     "Hit position of (isRHS, noEta, noPhi, x, y) = (%d, %d, %d, %d, %d) is (x, y, z) = (%f, %f, %f)",
                     isRHS, noEta, noPhi, col, row, total.x(), total.y(), total.z());

    return Vector3D(total.x(), total.y(), total.z());
  }

  // rear end position (used for the digitization)
  Vector3D GridDRcalo_k4geo::sipmPosition(const CellID& cID) const {
    int noEta = numEta(cID); // always positive since we replicated the RHS
    int noPhi = numPhi(cID);
    bool isRHS = IsRHS(cID);

    // first get a vector to the center of the wafer (per tower)
    DRparamBase_k4geo* paramBase = setParamBase(noEta);
    auto waferPos = paramBase->GetSipmLayerPos(noPhi);

    dd4hep::Position localPos = dd4hep::Position(0., 0., 0.);

    // get local coordinate
    // if IsSiPM is false, the local position is (0, 0, 0)
    // and the function returns the wafer's position
    if (IsSiPM(cID))
      localPos = dd4hep::Position(localPosition(cID));

    // rotate the local coordinate on par to the tower
    auto rot = paramBase->GetRotationZYX(noPhi);
    auto translation = rot * localPos;

    // total vector is sum of the vector to the wafer center + rotated local coordinate
    auto total = translation + waferPos;

    // if LHS rotate by 180 deg w.r.t. X axis (on par to the DRconstructor)
    if (!isRHS)
      total = dd4hep::RotationX(M_PI) * total;

    int row = y(cID);
    int col = x(cID);

    dd4hep::printout(dd4hep::VERBOSE, "GridDRcalo_k4geo::sipmPosition",
                     "SiPM position of (isRHS, noEta, noPhi, x, y) = (%d, %d, %d, %d, %d) is (x, y, z) = (%f, %f, %f)",
                     isRHS, noEta, noPhi, col, row, total.x(), total.y(), total.z());

    return Vector3D(total.x(), total.y(), total.z());
  }

  Vector3D GridDRcalo_k4geo::localPosition(const CellID cID) const {
    int numx = numX(cID);
    int numy = numY(cID);
    int x_ = x(cID);
    int y_ = y(cID);

    return localPosition(numx, numy, x_, y_);
  }

  Vector3D GridDRcalo_k4geo::localPosition(int numx, int numy, int x_, int y_) const {
    // Here x & y range from 0 to numx & numy (used for the geometry construction)
    float ptX = -m_gridSize * static_cast<float>(numx / 2) + static_cast<float>(x_) * m_gridSize +
                (numx % 2 == 0 ? m_gridSize / 2. : 0.);
    float ptY = -m_gridSize * static_cast<float>(numy / 2) + static_cast<float>(y_) * m_gridSize +
                (numy % 2 == 0 ? m_gridSize / 2. : 0.);

    return Vector3D(ptX, ptY, 0.);
  }

  /// determine the cell ID based on the position
  CellID GridDRcalo_k4geo::cellID(const Vector3D& localPosition, const Vector3D& globalPosition,
                                  const VolumeID& vID) const {
    // vID is assume to be the tower's
    // so we use z position instead to determine isRHS
    int systemId = static_cast<int>(decoder()->get(vID, m_systemIndex));
    int numx = numX(vID);
    int numy = numY(vID);
    bool isRHS = globalPosition.z() > 0.;

    auto localX = localPosition.x();
    auto localY = localPosition.y();

    int x = std::floor((localX + (numx % 2 == 0 ? 0. : m_gridSize / 2.)) / m_gridSize) + numx / 2;
    int y = std::floor((localY + (numy % 2 == 0 ? 0. : m_gridSize / 2.)) / m_gridSize) + numy / 2;

    return setCellID(isRHS, systemId, numEta(vID), numPhi(vID), x, y);
  }

  VolumeID GridDRcalo_k4geo::setVolumeID(int aSystem, int numEta, int numPhi) const {
    VolumeID systemId = static_cast<VolumeID>(aSystem);
    VolumeID numEtaId = static_cast<VolumeID>(numEta);
    VolumeID numPhiId = static_cast<VolumeID>(numPhi);
    VolumeID vID = 0;
    decoder()->set(vID, m_systemIndex, systemId);
    decoder()->set(vID, m_numEtaIndex, numEtaId);
    decoder()->set(vID, m_numPhiIndex, numPhiId);

    VolumeID module = 0; // Tower, SiPM layer attached to the tower, etc.
    decoder()->set(vID, m_moduleIndex, module);

    return vID;
  }

  CellID GridDRcalo_k4geo::setCellID(bool isRHS, int aSystem, int numEta, int numPhi, int x, int y) const {
    VolumeID systemId = static_cast<VolumeID>(aSystem);
    VolumeID numEtaId = static_cast<VolumeID>(numEta);
    VolumeID numPhiId = static_cast<VolumeID>(numPhi);
    VolumeID xId = static_cast<VolumeID>(x);
    VolumeID yId = static_cast<VolumeID>(y);
    VolumeID vID = 0;
    decoder()->set(vID, m_systemIndex, systemId);
    decoder()->set(vID, m_numEtaIndex, numEtaId);
    decoder()->set(vID, m_numPhiIndex, numPhiId);
    decoder()->set(vID, m_xIndex, xId);
    decoder()->set(vID, m_yIndex, yId);

    VolumeID module = 1; // Fiber, SiPM, etc.
    decoder()->set(vID, m_moduleIndex, module);

    VolumeID isCeren = IsCerenkov(x, y) ? 1 : 0;
    decoder()->set(vID, m_isCerenkovIndex, isCeren);

    VolumeID assemblyId = isRHS ? 0 : 1;
    decoder()->set(vID, m_assemblyIndex, assemblyId);

    return vID;
  }

  // neighbor finding algorithm regardless of the type of the channel (mixing scintillation and Cherenkov)
  void GridDRcalo_k4geo::neighbours(const CellID& cID, std::set<CellID>& neighbours) const {
    int systemId = static_cast<int>(decoder()->get(cID, m_systemIndex));
    int noEta = numEta(cID);
    int noPhi = numPhi(cID);
    int nX = x(cID); // col
    int nY = y(cID); // row
    int totX = numX(cID);
    int totY = numY(cID);
    bool isRHS = IsRHS(cID);
    bool isCeren = IsCerenkov(cID);

    DRparamBase_k4geo* paramBase = setParamBase(noEta);
    auto fl = paramBase->GetFullLengthFibers(noEta);
    int numZRot = paramBase->GetNumZRot();

    std::set<CellID> nb;
    float neighbourSize = paramBase->GetNeighborSize();
    int radiusFloor = static_cast<int>(std::floor(neighbourSize));

    // loop over (2n+1) x (2n+1) cells (n = radius)
    // and define the minimal neighborhood with the cells that dist <= radius
    for (int ix = 0; ix <= radiusFloor; ix++) {
      for (int iy = 0; iy <= radiusFloor; iy++) {
        // skip the cell of interest
        if (ix == 0 && iy == 0)
          continue;

        if (static_cast<float>(ix * ix + iy * iy) <= neighbourSize * neighbourSize) {
          // need to cover four quadrants of (+,+), (-,-), (+,-), (-,+)
          nb.insert(setCellID(isRHS, systemId, noEta, noPhi, nX + ix, nY + iy));
          nb.insert(setCellID(isRHS, systemId, noEta, noPhi, nX - ix, nY - iy));
          // duplicated cellID when either ix or iy = 0 but std::set will handle it
          nb.insert(setCellID(isRHS, systemId, noEta, noPhi, nX + ix, nY - iy));
          nb.insert(setCellID(isRHS, systemId, noEta, noPhi, nX - ix, nY + iy));
        }
      }
    }

    // function to remove different channel in the set
    auto removeDifferentChannel = [this](bool isC, std::set<CellID>& input) {
      for (auto it = input.begin(); it != input.end();) {
        if (isC != IsCerenkov(*it))
          it = input.erase(it);
        else
          ++it;
      }
    };

    // now handle the (tower) edge cases
    // margin to switch the definition of the neighborhood at the edge of the tower
    int margin = paramBase->GetMargin();

    // if the seed fiber is not on the edge of the tower, return the minimal neighborhood
    if (nX > fl.cmin + margin && nX < fl.cmax - margin && nY > fl.rmin + margin && nY < fl.rmax - margin) {
      if (m_removeDifferentCh)
        removeDifferentChannel(isCeren, nb);

      neighbours = nb;
      return;
    }

    // WARNING the following edge stitching part has very strict assumption
    // on the way how geometry is constructed
    // e.g. how do we replicate the RHS to LHS
    // if someone wants to change the geometry
    // then always check that the neighboring algorithm is not broken

    // now the seed fiber is on the edge of the tower
    // then stitch the short length fibers
    // up to the closest full length fiber in the neighboring tower
    // first stitch phi direction first (it's easier)
    auto modulo = [](int in, int n) -> int { return (in % n + n) % n; };

    if (nX >= fl.cmax - margin) { // to east
      // same tower
      for (int idx = nX + 1; idx < totX; idx++)
        nb.insert(setCellID(isRHS, systemId, noEta, noPhi, idx, nY));

      int nextPhi = modulo(noPhi - 1, numZRot);

      // next tower
      for (int idx = 0; idx <= fl.cmin + margin; idx++)
        nb.insert(setCellID(isRHS, systemId, noEta, nextPhi, idx, nY));
    }

    if (nX <= fl.cmin + margin) { // to west
      // same tower
      for (int idx = nX - 1; idx >= 0; idx--)
        nb.insert(setCellID(isRHS, systemId, noEta, noPhi, idx, nY));

      int nextPhi = modulo(noPhi + 1, numZRot);

      // next tower
      for (int idx = totX - 1; idx >= fl.cmax - margin; idx--)
        nb.insert(setCellID(isRHS, systemId, noEta, nextPhi, idx, nY));
    }

    // now stitch in eta direction
    // but first we need an explicit treatment at eta=0
    // because we use rotation instead of reflection
    // hence cellID is not consistent
    bool isBarrelCenter = (noEta == 0 && nY >= fl.rmax - margin);

    if (isBarrelCenter) {
      // same tower
      for (int idx = nY + 1; idx < totY; idx++)
        nb.insert(setCellID(isRHS, systemId, noEta, noPhi, nX, idx));

      // since we rotate by 180 deg, nextPhi is numZRot - noPhi
      int nextPhi = modulo(numZRot - noPhi, numZRot);

      // same applies for nX
      int nextX = totX - nX;

      // next tower
      for (int idx = totY - 1; idx >= fl.rmax - margin; idx--)
        nb.insert(setCellID(!isRHS, systemId, noEta, nextPhi, nextX, idx));

      if (m_removeDifferentCh)
        removeDifferentChannel(isCeren, nb);

      neighbours = nb;
      return;
    }

    // now we forget about the eta=0 case
    // and treat the north & south bound case
    if (nY >= fl.rmax) { // to north
      // same tower
      for (int idx = nY + 1; idx < totY; idx++)
        nb.insert(setCellID(isRHS, systemId, noEta, noPhi, nX, idx));

      // for different noEta rmin and rmax can be different
      // also protect from map::at exception at the barrel-endcap boundary
      auto flNext = paramBase->unsignedTowerNo(noEta) == m_paramBarrel->GetTotTowerNum()
                        ? m_paramBarrel->GetFullLengthFibers(noEta - 1)
                        : paramBase->GetFullLengthFibers(noEta - 1);

      // next tower
      for (int idx = 0; idx <= flNext.rmin + margin; idx++)
        nb.insert(setCellID(isRHS, systemId, noEta - 1, noPhi, nX, idx));
    }

    if (nY <= fl.rmin) { // to south
      // same tower
      for (int idx = nY - 1; idx >= 0; idx--)
        nb.insert(setCellID(isRHS, systemId, noEta, noPhi, nX, idx));

      // next tower
      // in principle totY is also different
      // but to southbound the next totY is always smaller
      // than the current one
      // also protect from looking for a tower with numEta greater than the total # of towers
      if (noEta + 1 < m_paramEndcap->GetTotTowerNum() + m_paramBarrel->GetTotTowerNum()) {
        // protect from map::at exception at the barrel-endcap boundary
        auto flNext = paramBase->unsignedTowerNo(noEta) + 1 == m_paramBarrel->GetTotTowerNum()
                          ? m_paramEndcap->GetFullLengthFibers(noEta + 1)
                          : paramBase->GetFullLengthFibers(noEta + 1);

        for (int idx = totY - 1; idx >= flNext.rmax - margin; idx--)
          nb.insert(setCellID(isRHS, systemId, noEta + 1, noPhi, nX, idx));
      }
    }

    // finalize
    if (m_removeDifferentCh)
      removeDifferentChannel(isCeren, nb);

    neighbours = nb;
    return;
  }

  // Get the identifier number of a mother tower in eta or phi direction
  int GridDRcalo_k4geo::numEta(const CellID aCellID) const {
    VolumeID numEta = static_cast<VolumeID>(decoder()->get(aCellID, m_numEtaIndex));
    return static_cast<int>(numEta);
  }

  int GridDRcalo_k4geo::numPhi(const CellID aCellID) const {
    VolumeID numPhi = static_cast<VolumeID>(decoder()->get(aCellID, m_numPhiIndex));
    return static_cast<int>(numPhi);
  }

  // Get the total number of SiPMs of the mother tower in x or y direction (local coordinate)
  int GridDRcalo_k4geo::numX(const CellID aCellID) const {
    int noEta = numEta(aCellID);

    DRparamBase_k4geo* paramBase = setParamBase(noEta);

    int noX =
        static_cast<int>(std::floor((paramBase->GetTl2() * 2. - m_sipmSize / 2.) / m_gridSize)) + 1; // in phi direction

    return noX;
  }

  int GridDRcalo_k4geo::numY(const CellID aCellID) const {
    int noEta = numEta(aCellID);

    DRparamBase_k4geo* paramBase = setParamBase(noEta);

    int noY =
        static_cast<int>(std::floor((paramBase->GetH2() * 2. - m_sipmSize / 2.) / m_gridSize)) + 1; // in eta direction

    return noY;
  }

  // Get the identifier number of a SiPM in x or y direction (local coordinate)
  int GridDRcalo_k4geo::x(const CellID aCellID) const { // approx phi direction
    VolumeID x = static_cast<VolumeID>(decoder()->get(aCellID, m_xIndex));
    return static_cast<int>(x);
  }
  int GridDRcalo_k4geo::y(const CellID aCellID) const { // approx eta direction
    VolumeID y = static_cast<VolumeID>(decoder()->get(aCellID, m_yIndex));
    return static_cast<int>(y);
  }

  bool GridDRcalo_k4geo::IsCerenkov(const CellID aCellID) const {
    VolumeID isCeren = static_cast<VolumeID>(decoder()->get(aCellID, m_isCerenkovIndex));
    return static_cast<bool>(isCeren);
  }
  // Identify if the fiber is Cerenkov or scintillation by its column and row number
  bool GridDRcalo_k4geo::IsCerenkov(int col, int row) const {
    bool isCeren = false;
    if (col % 2 == 1) {
      isCeren = !isCeren;
    }
    if (row % 2 == 1) {
      isCeren = !isCeren;
    }
    return isCeren;
  }

  bool GridDRcalo_k4geo::IsTower(const CellID aCellID) const {
    VolumeID module = static_cast<VolumeID>(decoder()->get(aCellID, m_moduleIndex));
    return module == 0;
  }

  bool GridDRcalo_k4geo::IsSiPM(const CellID aCellID) const {
    VolumeID module = static_cast<VolumeID>(decoder()->get(aCellID, m_moduleIndex));
    return module == 1;
  }

  bool GridDRcalo_k4geo::IsRHS(const CellID aCellID) const {
    VolumeID assembly = static_cast<VolumeID>(decoder()->get(aCellID, m_assemblyIndex));
    return assembly == 0;
  }

  int GridDRcalo_k4geo::getLast32bits(const CellID aCellID) const {
    CellID aId64 = aCellID >> sizeof(int) * CHAR_BIT;
    int aId32 = (int)aId64;

    return aId32;
  }

  CellID GridDRcalo_k4geo::convertLast32to64(const int aId32) const {
    CellID aId64 = (CellID)aId32;
    aId64 <<= sizeof(int) * CHAR_BIT;

    return aId64;
  }

  DRparamBase_k4geo* GridDRcalo_k4geo::setParamBase(int noEta) const {
    DRparamBase_k4geo* paramBase = nullptr;

    if (m_paramEndcap->unsignedTowerNo(noEta) >= m_paramBarrel->GetTotTowerNum())
      paramBase = static_cast<DRparamBase_k4geo*>(m_paramEndcap);
    else
      paramBase = static_cast<DRparamBase_k4geo*>(m_paramBarrel);

    if (paramBase->GetCurrentTowerNum() == noEta)
      return paramBase;

    // This should not be called while building detector geometry
    if (!paramBase->IsFinalized())
      throw std::runtime_error("GridDRcalo_k4geo::position should not be called while building detector geometry!");

    paramBase->SetDeltaThetaByTowerNo(noEta, m_paramBarrel->GetTotTowerNum());
    paramBase->SetThetaOfCenterByTowerNo(noEta, m_paramBarrel->GetTotTowerNum());
    paramBase->SetIsRHSByTowerNo(noEta);
    paramBase->SetCurrentTowerNum(noEta);
    paramBase->init();

    return paramBase;
  }

} // namespace DDSegmentation
} // namespace dd4hep
