//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#include "detectorSegmentations/SCEPCal_MainSegmentation_k4geo.h"
#include <climits>
#include <cmath>
#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

SCEPCal_MainSegmentation_k4geo::SCEPCal_MainSegmentation_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) {
    _type = "SCEPCal_MainSegmentation_k4geo";
    _description = "SCEPCal main layer segmentation";
    registerIdentifier("identifier_system", "Cell ID identifier for numSystem", fSystemId, "system");
    registerIdentifier("identifier_phi", "Cell ID identifier for numPhi", fPhiId, "phi");
    registerIdentifier("identifier_eta", "Cell ID identifier for numEta", fEtaId, "eta");
    registerIdentifier("identifier_gamma", "Cell ID identifier for numGamma", fGammaId, "gamma");
    registerIdentifier("identifier_epsilon", "Cell ID identifier for numEpsilon", fEpsilonId, "epsilon");
    registerIdentifier("identifier_depth", "Cell ID identifier for numDepth", fDepthId, "depth");
}

SCEPCal_MainSegmentation_k4geo::SCEPCal_MainSegmentation_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) {
    _type = "SCEPCal_MainSegmentation_k4geo";
    _description = "SCEPCal main layer segmentation";
    registerIdentifier("identifier_system", "Cell ID identifier for numSystem", fSystemId, "system");
    registerIdentifier("identifier_phi", "Cell ID identifier for Phi", fPhiId, "phi");
    registerIdentifier("identifier_eta", "Cell ID identifier for Eta", fEtaId, "eta");
    registerIdentifier("identifier_gamma", "Cell ID identifier for Gamma", fGammaId, "gamma");
    registerIdentifier("identifier_epsilon", "Cell ID identifier for Epsilon", fEpsilonId, "epsilon");
    registerIdentifier("identifier_depth", "Cell ID identifier for Depth", fDepthId, "depth");
}

SCEPCal_MainSegmentation_k4geo::~SCEPCal_MainSegmentation_k4geo() {}

Vector3D SCEPCal_MainSegmentation_k4geo::position(const CellID& cID) const {
    return myPosition((int)cID);
};

Vector3D SCEPCal_MainSegmentation_k4geo::myPosition(const int& cellId) const {

    if (fPositionOf.find(cellId) != fPositionOf.end()) {
        return fPositionOf.find(cellId)->second;
    }

    return Vector3D(0,0,0);
}
}
}
