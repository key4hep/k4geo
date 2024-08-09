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
  int noEta = numEta(cID);
  int noPhi = numPhi(cID);

  DRparamBase_k4geo* paramBase = setParamBase(noEta);

  auto transformA = paramBase->GetSipmTransform3D(noPhi);
  dd4hep::Position localPos = dd4hep::Position(0.,0.,0.);
  if ( IsSiPM(cID) ) localPos = dd4hep::Position( localPosition(cID) );

  dd4hep::RotationZYX rot = dd4hep::RotationZYX(M_PI, 0., 0.); // AdHoc rotation, potentially bug
  dd4hep::Transform3D transformB = dd4hep::Transform3D(rot,localPos);
  auto total = transformA*transformB;

  dd4hep::Position translation = dd4hep::Position(0.,0.,0.);
  total.GetTranslation(translation);

  return Vector3D(translation.x(),translation.y(),translation.z());
}

Vector3D GridDRcalo_k4geo::localPosition(const CellID& cID) const {
  int numx = numX(cID);
  int numy = numY(cID);
  int x_ = x(cID);
  int y_ = y(cID);

  return localPosition(numx,numy,x_,y_);
}

Vector3D GridDRcalo_k4geo::localPosition(int numx, int numy, int x_, int y_) const {
  // Here x & y range from 0 to numx & numy (used for the geometry construction)
  float ptX = -fGridSize*static_cast<float>(numx/2) + static_cast<float>(x_)*fGridSize + ( numx%2==0 ? fGridSize/2. : 0. );
  float ptY = -fGridSize*static_cast<float>(numy/2) + static_cast<float>(y_)*fGridSize + ( numy%2==0 ? fGridSize/2. : 0. );

  return Vector3D(ptX,ptY,0.);
}

/// determine the cell ID based on the position
CellID GridDRcalo_k4geo::cellID(const Vector3D& localPosition, const Vector3D& globalPosition, const VolumeID& vID) const {
  int numx = numX(vID);
  int numy = numY(vID);

  auto localX = localPosition.x();
  auto localY = localPosition.y();

  int x = std::floor( ( localX + ( numx%2==0 ? 0. : fGridSize/2. ) ) / fGridSize ) + numx/2;
  int y = std::floor( ( localY + ( numy%2==0 ? 0. : fGridSize/2. ) ) / fGridSize ) + numy/2;

  int signedNumEta = (globalPosition.z() >= 0) ? numEta(vID) : (-numEta(vID) - 1);
  return setCellID( signedNumEta, numPhi(vID), x, y );
}

VolumeID GridDRcalo_k4geo::setVolumeID(int numEta, int numPhi) const {
  VolumeID numEtaId = static_cast<VolumeID>(numEta);
  VolumeID numPhiId = static_cast<VolumeID>(numPhi);
  VolumeID vID = 0;
  _decoder->set(vID, fNumEtaId, numEtaId);
  _decoder->set(vID, fNumPhiId, numPhiId);

  VolumeID module = 0; // Tower, SiPM layer attached to the tower, etc.
  _decoder->set(vID, fModule, module);

  return vID;
}

CellID GridDRcalo_k4geo::setCellID(int numEta, int numPhi, int x, int y) const {
  VolumeID numEtaId = static_cast<VolumeID>(numEta);
  VolumeID numPhiId = static_cast<VolumeID>(numPhi);
  VolumeID xId = static_cast<VolumeID>(x);
  VolumeID yId = static_cast<VolumeID>(y);
  VolumeID vID = 0;
  _decoder->set(vID, fNumEtaId, numEtaId);
  _decoder->set(vID, fNumPhiId, numPhiId);
  _decoder->set(vID, fXId, xId);
  _decoder->set(vID, fYId, yId);

  VolumeID module = 1; // Fiber, SiPM, etc.
  _decoder->set(vID, fModule, module);

  VolumeID isCeren = IsCerenkov(x,y) ? 1 : 0;
  _decoder->set(vID, fIsCerenkovId, isCeren);

  VolumeID assemblyId = (numEta >= 0) ? 0 : 1;
  _decoder->set(vID, fAssemblyId, assemblyId);

  return vID;
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

  int noX = static_cast<int>( std::floor( ( paramBase->GetTl2()*2. - fSipmSize/2. )/fGridSize ) ) + 1; // in phi direction

  return noX;
}

int GridDRcalo_k4geo::numY(const CellID& aCellID) const {
  int noEta = numEta(aCellID);

  DRparamBase_k4geo* paramBase = setParamBase(noEta);

  int noY = static_cast<int>( std::floor( ( paramBase->GetH2()*2. - fSipmSize/2. )/fGridSize ) ) + 1; // in eta direction

  return noY;
}

// Get the identifier number of a SiPM in x or y direction (local coordinate)
int GridDRcalo_k4geo::x(const CellID& aCellID) const { // approx eta direction
  VolumeID x = static_cast<VolumeID>(_decoder->get(aCellID, fXId));
  return static_cast<int>(x);
}
int GridDRcalo_k4geo::y(const CellID& aCellID) const { // approx phi direction
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
  if ( col%2 == 1 ) { isCeren = !isCeren; }
  if ( row%2 == 1 ) { isCeren = !isCeren; }
  return isCeren;
}

bool GridDRcalo_k4geo::IsTower(const CellID& aCellID) const {
  VolumeID module = static_cast<VolumeID>(_decoder->get(aCellID, fModule));
  return module==0;
}

bool GridDRcalo_k4geo::IsSiPM(const CellID& aCellID) const {
  VolumeID module = static_cast<VolumeID>(_decoder->get(aCellID, fModule));
  return module==1;
}

int GridDRcalo_k4geo::getLast32bits(const CellID& aCellID) const {
  CellID aId64 = aCellID >> sizeof(int)*CHAR_BIT;
  int aId32 = (int)aId64;

  return aId32;
}

CellID GridDRcalo_k4geo::convertLast32to64(const int aId32) const {
  CellID aId64 = (CellID)aId32;
  aId64 <<= sizeof(int)*CHAR_BIT;

  return aId64;
}

DRparamBase_k4geo* GridDRcalo_k4geo::setParamBase(int noEta) const {
  DRparamBase_k4geo* paramBase = nullptr;

  if ( fParamEndcap->unsignedTowerNo(noEta) >= fParamBarrel->GetTotTowerNum() ) paramBase = static_cast<DRparamBase_k4geo*>(fParamEndcap);
  else paramBase = static_cast<DRparamBase_k4geo*>(fParamBarrel);

  if ( paramBase->GetCurrentTowerNum()==noEta ) return paramBase;

  // This should not be called while building detector geometry
  if (!paramBase->IsFinalized()) throw std::runtime_error("GridDRcalo_k4geo::position should not be called while building detector geometry!");

  paramBase->SetDeltaThetaByTowerNo(noEta, fParamBarrel->GetTotTowerNum());
  paramBase->SetThetaOfCenterByTowerNo(noEta, fParamBarrel->GetTotTowerNum());
  paramBase->SetIsRHSByTowerNo(noEta);
  paramBase->SetCurrentTowerNum(noEta);
  paramBase->init();

  return paramBase;
}

}
}
