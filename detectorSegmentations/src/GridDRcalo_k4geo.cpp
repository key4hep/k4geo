#include "detectorSegmentations/GridDRcalo_k4geo.h"

#include <climits>
#include <cmath>
#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  GridDRcalo_k4geo::GridDRcalo_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) {
    // define type and description
    _type = "GridDRcalo_k4geo";
    _description = "DRcalo segmentation based on the tower / (Cerenkov or Scintillation) fiber / SiPM hierarchy";

    // register all necessary parameters
    registerIdentifier("identifier_assembly", "Cell ID identifier for assembly", fAssemblyId, "assembly");
    registerIdentifier("identifier_eta", "Cell ID identifier for numEta", fNumEtaId, "eta");
    registerIdentifier("identifier_phi", "Cell ID identifier for numPhi", fNumPhiId, "phi");
    registerIdentifier("identifier_x", "Cell ID identifier for x", fXId, "x");
    registerIdentifier("identifier_y", "Cell ID identifier for y", fYId, "y");
    registerIdentifier("identifier_IsCerenkov", "Cell ID identifier for IsCerenkov", fIsCerenkovId, "c");
    registerIdentifier("identifier_module", "Cell ID identifier for module", fModule, "module");

    fParamBarrel = new DRparamBarrel_k4geo();
    fParamEndcap = new DRparamEndcap_k4geo();
  }

  GridDRcalo_k4geo::GridDRcalo_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) {
    // define type and description
    _type = "GridDRcalo_k4geo";
    _description = "DRcalo segmentation based on the tower / (Cerenkov or Scintillation) fiber / SiPM hierarchy";

    // register all necessary parameters
    registerIdentifier("identifier_assembly", "Cell ID identifier for assembly", fAssemblyId, "assembly");
    registerIdentifier("identifier_eta", "Cell ID identifier for numEta", fNumEtaId, "eta");
    registerIdentifier("identifier_phi", "Cell ID identifier for numPhi", fNumPhiId, "phi");
    registerIdentifier("identifier_x", "Cell ID identifier for x", fXId, "x");
    registerIdentifier("identifier_y", "Cell ID identifier for y", fYId, "y");
    registerIdentifier("identifier_IsCerenkov", "Cell ID identifier for IsCerenkov", fIsCerenkovId, "c");
    registerIdentifier("identifier_module", "Cell ID identifier for module", fModule, "module");

    fParamBarrel = new DRparamBarrel_k4geo();
    fParamEndcap = new DRparamEndcap_k4geo();
  }

  GridDRcalo_k4geo::~GridDRcalo_k4geo() {
    delete fParamBarrel;
    delete fParamEndcap;
  }

  Vector3D GridDRcalo_k4geo::position(const CellID& cID) const {
    int noEta = numEta(cID); // always positive since we replicated the RHS
    int noPhi = numPhi(cID);
    bool isRHS = IsRHS(cID);

    // first get a vector to the center of the wafer (per tower)
    DRparamBase_k4geo* paramBase = setParamBase(noEta);
    auto waferPos = paramBase->GetSipmLayerPos(noPhi);

    dd4hep::Position localPos = dd4hep::Position(0., 0., 0.);

    // get local coordinate
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

    return Vector3D(total.x(), total.y(), total.z());
  }

  Vector3D GridDRcalo_k4geo::localPosition(const CellID& cID) const {
    int numx = numX(cID);
    int numy = numY(cID);
    int x_ = x(cID);
    int y_ = y(cID);

    return localPosition(numx, numy, x_, y_);
  }

  Vector3D GridDRcalo_k4geo::localPosition(int numx, int numy, int x_, int y_) const {
    // Here x & y range from 0 to numx & numy (used for the geometry construction)
    float ptX = -fGridSize * static_cast<float>(numx / 2) + static_cast<float>(x_) * fGridSize +
                (numx % 2 == 0 ? fGridSize / 2. : 0.);
    float ptY = -fGridSize * static_cast<float>(numy / 2) + static_cast<float>(y_) * fGridSize +
                (numy % 2 == 0 ? fGridSize / 2. : 0.);

    return Vector3D(ptX, ptY, 0.);
  }

  /// determine the cell ID based on the position
  CellID GridDRcalo_k4geo::cellID(const Vector3D& localPosition, const Vector3D& globalPosition,
                                  const VolumeID& vID) const {
    // vID is assume to be the tower's
    // so we use z position instead to determine isRHS
    int systemId = static_cast<int>(_decoder->get(vID, "system"));
    int numx = numX(vID);
    int numy = numY(vID);
    bool isRHS = globalPosition.z() > 0.;

    auto localX = localPosition.x();
    auto localY = localPosition.y();

    int x = std::floor((localX + (numx % 2 == 0 ? 0. : fGridSize / 2.)) / fGridSize) + numx / 2;
    int y = std::floor((localY + (numy % 2 == 0 ? 0. : fGridSize / 2.)) / fGridSize) + numy / 2;

    return setCellID(isRHS, systemId, numEta(vID), numPhi(vID), x, y);
  }

  VolumeID GridDRcalo_k4geo::setVolumeID(int aSystem, int numEta, int numPhi) const {
    VolumeID systemId = static_cast<VolumeID>(aSystem);
    VolumeID numEtaId = static_cast<VolumeID>(numEta);
    VolumeID numPhiId = static_cast<VolumeID>(numPhi);
    VolumeID vID = 0;
    _decoder->set(vID, "system", systemId);
    _decoder->set(vID, fNumEtaId, numEtaId);
    _decoder->set(vID, fNumPhiId, numPhiId);

    VolumeID module = 0; // Tower, SiPM layer attached to the tower, etc.
    _decoder->set(vID, fModule, module);

    return vID;
  }

  CellID GridDRcalo_k4geo::setCellID(bool isRHS, int aSystem, int numEta, int numPhi, int x, int y) const {
    VolumeID systemId = static_cast<VolumeID>(aSystem);
    VolumeID numEtaId = static_cast<VolumeID>(numEta);
    VolumeID numPhiId = static_cast<VolumeID>(numPhi);
    VolumeID xId = static_cast<VolumeID>(x);
    VolumeID yId = static_cast<VolumeID>(y);
    VolumeID vID = 0;
    _decoder->set(vID, "system", systemId);
    _decoder->set(vID, fNumEtaId, numEtaId);
    _decoder->set(vID, fNumPhiId, numPhiId);
    _decoder->set(vID, fXId, xId);
    _decoder->set(vID, fYId, yId);

    VolumeID module = 1; // Fiber, SiPM, etc.
    _decoder->set(vID, fModule, module);

    VolumeID isCeren = IsCerenkov(x, y) ? 1 : 0;
    _decoder->set(vID, fIsCerenkovId, isCeren);

    VolumeID assemblyId = isRHS ? 0 : 1;
    _decoder->set(vID, fAssemblyId, assemblyId);

    return vID;
  }

  void GridDRcalo_k4geo::neighbours(const CellID& cID, std::set<CellID>& neighbours) const {
    int systemId = static_cast<int>(_decoder->get(cID, "system"));
    int noEta = numEta(cID);
    int noPhi = numPhi(cID);
    int nX = x(cID); // col
    int nY = y(cID); // row
    int totX = numX(cID);
    int totY = numY(cID);
    bool isCeren = IsCerenkov(cID);
    bool isRHS = IsRHS(cID);

    DRparamBase_k4geo* paramBase = setParamBase(noEta);
    auto fl = paramBase->GetFullLengthFibers(noEta);
    int numZRot = paramBase->GetNumZRot();

    // First we look for the closest (dist=sqrt(2))
    // and the second-closest (dist=2) cells in the checkerboard
    // in the same tower
    // dist=sqrt(2) cells
    auto northEast = setCellID(isRHS, systemId, noEta, noPhi, nX + 1, nY + 1);
    auto southEast = setCellID(isRHS, systemId, noEta, noPhi, nX + 1, nY - 1);
    auto southWest = setCellID(isRHS, systemId, noEta, noPhi, nX - 1, nY - 1);
    auto northWest = setCellID(isRHS, systemId, noEta, noPhi, nX - 1, nY + 1);
    // dist=2 cells
    auto north = setCellID(isRHS, systemId, noEta, noPhi, nX, nY + 2);
    auto east = setCellID(isRHS, systemId, noEta, noPhi, nX + 2, nY);
    auto south = setCellID(isRHS, systemId, noEta, noPhi, nX, nY - 2);
    auto west = setCellID(isRHS, systemId, noEta, noPhi, nX - 2, nY);

    // the minimal neighborhood
    std::set<CellID> nb = {north, northEast, east, southEast, south, southWest, west, northWest};

    // function to remove different channel in the set
    auto removeDifferentChannel = [this](bool isC, std::set<CellID>& input) {
      for (auto it = input.begin(); it != input.end();) {
        if (isC != IsCerenkov(*it))
          it = input.erase(it);
        else
          ++it;
      }
    };

    // margin to switch the definiton of the neighborhood at the edge of the tower
    int margin = 2;

    // if the seed fiber is not on the edge of the tower, return the minimal neighborhood
    if (nX > fl.cmin + margin && nX < fl.cmax - margin && nY > fl.rmin + margin && nY < fl.rmax - margin) {
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
      auto flNext = paramBase->GetFullLengthFibers(noEta - 1);

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
      if (noEta + 1 < fParamEndcap->GetTotTowerNum() + fParamBarrel->GetTotTowerNum()) {
        auto flNext = paramBase->GetFullLengthFibers(noEta + 1);

        for (int idx = totY - 1; idx >= flNext.rmax - margin; idx--)
          nb.insert(setCellID(isRHS, systemId, noEta + 1, noPhi, nX, idx));
      }
    }

    // finalize
    removeDifferentChannel(isCeren, nb);
    neighbours = nb;
    return;
  }

  // Get the identifier number of a mother tower in eta or phi direction
  int GridDRcalo_k4geo::numEta(const CellID& aCellID) const {
    VolumeID numEta = static_cast<VolumeID>(_decoder->get(aCellID, fNumEtaId));
    return static_cast<int>(numEta);
  }

  int GridDRcalo_k4geo::numPhi(const CellID& aCellID) const {
    VolumeID numPhi = static_cast<VolumeID>(_decoder->get(aCellID, fNumPhiId));
    return static_cast<int>(numPhi);
  }

  // Get the total number of SiPMs of the mother tower in x or y direction (local coordinate)
  int GridDRcalo_k4geo::numX(const CellID& aCellID) const {
    int noEta = numEta(aCellID);

    DRparamBase_k4geo* paramBase = setParamBase(noEta);

    int noX =
        static_cast<int>(std::floor((paramBase->GetTl2() * 2. - fSipmSize / 2.) / fGridSize)) + 1; // in phi direction

    return noX;
  }

  int GridDRcalo_k4geo::numY(const CellID& aCellID) const {
    int noEta = numEta(aCellID);

    DRparamBase_k4geo* paramBase = setParamBase(noEta);

    int noY =
        static_cast<int>(std::floor((paramBase->GetH2() * 2. - fSipmSize / 2.) / fGridSize)) + 1; // in eta direction

    return noY;
  }

  // Get the identifier number of a SiPM in x or y direction (local coordinate)
  int GridDRcalo_k4geo::x(const CellID& aCellID) const { // approx phi direction
    VolumeID x = static_cast<VolumeID>(_decoder->get(aCellID, fXId));
    return static_cast<int>(x);
  }
  int GridDRcalo_k4geo::y(const CellID& aCellID) const { // approx eta direction
    VolumeID y = static_cast<VolumeID>(_decoder->get(aCellID, fYId));
    return static_cast<int>(y);
  }

  bool GridDRcalo_k4geo::IsCerenkov(const CellID& aCellID) const {
    VolumeID isCeren = static_cast<VolumeID>(_decoder->get(aCellID, fIsCerenkovId));
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

  bool GridDRcalo_k4geo::IsTower(const CellID& aCellID) const {
    VolumeID module = static_cast<VolumeID>(_decoder->get(aCellID, fModule));
    return module == 0;
  }

  bool GridDRcalo_k4geo::IsSiPM(const CellID& aCellID) const {
    VolumeID module = static_cast<VolumeID>(_decoder->get(aCellID, fModule));
    return module == 1;
  }

  bool GridDRcalo_k4geo::IsRHS(const CellID& aCellID) const {
    VolumeID assembly = static_cast<VolumeID>(_decoder->get(aCellID, fAssemblyId));
    return assembly == 0;
  }

  int GridDRcalo_k4geo::getLast32bits(const CellID& aCellID) const {
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

    if (fParamEndcap->unsignedTowerNo(noEta) >= fParamBarrel->GetTotTowerNum())
      paramBase = static_cast<DRparamBase_k4geo*>(fParamEndcap);
    else
      paramBase = static_cast<DRparamBase_k4geo*>(fParamBarrel);

    if (paramBase->GetCurrentTowerNum() == noEta)
      return paramBase;

    // This should not be called while building detector geometry
    if (!paramBase->IsFinalized())
      throw std::runtime_error("GridDRcalo_k4geo::position should not be called while building detector geometry!");

    paramBase->SetDeltaThetaByTowerNo(noEta, fParamBarrel->GetTotTowerNum());
    paramBase->SetThetaOfCenterByTowerNo(noEta, fParamBarrel->GetTotTowerNum());
    paramBase->SetIsRHSByTowerNo(noEta);
    paramBase->SetCurrentTowerNum(noEta);
    paramBase->init();

    return paramBase;
  }

} // namespace DDSegmentation
} // namespace dd4hep
