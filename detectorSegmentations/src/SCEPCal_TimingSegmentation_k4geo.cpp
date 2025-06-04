//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#include "detectorSegmentations/SCEPCal_TimingSegmentation_k4geo.h"
#include <climits>
#include <cmath>
#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

SCEPCal_TimingSegmentation_k4geo::SCEPCal_TimingSegmentation_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) {
    _type = "SCEPCal_TimingSegmentation_k4geo";
    _description = "SCEPCal timing layer segmentation";
    registerIdentifier("identifier_system", "Cell ID identifier for System", fSystemId, "system");
    registerIdentifier("identifier_phi",    "Cell ID identifier for Phi",    fPhiId,    "phi");
    registerIdentifier("identifier_theta",  "Cell ID identifier for Theta",  fThetaId,  "theta");
    registerIdentifier("identifier_gamma",  "Cell ID identifier for Gamma",  fGammaId,  "gamma");
}

SCEPCal_TimingSegmentation_k4geo::SCEPCal_TimingSegmentation_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) {
    _type = "SCEPCal_TimingSegmentation_k4geo";
    _description = "SCEPCal timing layer segmentation";
    registerIdentifier("identifier_system", "Cell ID identifier for System", fSystemId, "system");
    registerIdentifier("identifier_phi",    "Cell ID identifier for Phi",    fPhiId,    "phi");
    registerIdentifier("identifier_theta",  "Cell ID identifier for Theta",  fThetaId,  "theta");
    registerIdentifier("identifier_gamma",  "Cell ID identifier for Gamma",  fGammaId,  "gamma");
}

SCEPCal_TimingSegmentation_k4geo::~SCEPCal_TimingSegmentation_k4geo() {}

Vector3D SCEPCal_TimingSegmentation_k4geo::position(const CellID& cellId) const {
    if (fPositionOf.find(cellId) != fPositionOf.end()) {
        return fPositionOf.find(cellId)->second;
    }

    return Vector3D(0,0,0);
}
}
}
