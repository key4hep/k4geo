//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#include "detectorSegmentations/SCEPCalSegmentation_k4geo.h"
#include <climits>
#include <cmath>
#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

SCEPCalSegmentation_k4geo::SCEPCalSegmentation_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) {
    _type = "SCEPCalSegmentation_k4geo";
    _description = "SCEPCal segmentation based on side/eta/phi/depth/S/C";
    registerIdentifier("identifier_system", "Cell ID identifier for numSystem", fSystemId, "system");
    registerIdentifier("identifier_eta", "Cell ID identifier for numEta", fEtaId, "eta");
    registerIdentifier("identifier_phi", "Cell ID identifier for numPhi", fPhiId, "phi");
    registerIdentifier("identifier_depth", "Cell ID identifier for numDepth", fDepthId, "depth");
}

SCEPCalSegmentation_k4geo::SCEPCalSegmentation_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) {
    _type = "SCEPCalSegmentation_k4geo";
    _description = "SCEPCal segmentation based on side/eta/phi/depth/S/C";
    registerIdentifier("identifier_system", "Cell ID identifier for numSystem", fSystemId, "system");
    registerIdentifier("identifier_eta", "Cell ID identifier for Eta", fEtaId, "eta");
    registerIdentifier("identifier_phi", "Cell ID identifier for Phi", fPhiId, "phi");
    registerIdentifier("identifier_depth", "Cell ID identifier for Depth", fDepthId, "depth");
}

SCEPCalSegmentation_k4geo::~SCEPCalSegmentation_k4geo() {}

Vector3D SCEPCalSegmentation_k4geo::position(const CellID& cID) const {
    return myPosition(cID);
};

Vector3D SCEPCalSegmentation_k4geo::myPosition(const CellID& cID) const {

    int copyNum = (int)cID;

    double Fdz           =getFdz();
    double Rdz           =getRdz();
    double nomfw         =getnomfw();
    double nomth         =getnomth();
    double EBz           =getEBz();
    double Rin           =getRin();
    double sipmth        =getSipmth();
    int PHI_SEGMENTS     =getphiSegments();
    int N_PROJECTIVE_FILL=getNProjectiveFill();
    int system           =System(copyNum);
    int nEta_in          =Eta(copyNum);
    int nPhi_in          =Phi(copyNum);
    int nDepth_in        =Depth(copyNum);

    double  D_PHI_GLOBAL        =2*M_PI/PHI_SEGMENTS;
    double  PROJECTIVE_GAP      =(N_PROJECTIVE_FILL*nomfw)/2;
    // double  THETA_SIZE_BARREL   =atan(EBz/Rin);
    double  THETA_SIZE_ENDCAP   =atan(Rin/EBz);
    int     N_THETA_BARREL      =2*floor(EBz/nomfw);
    int     N_THETA_ENDCAP      =floor(Rin/nomfw);
    double  D_THETA_BARREL      =(M_PI-2*THETA_SIZE_ENDCAP)/(N_THETA_BARREL);
    double  D_THETA_ENDCAP      =THETA_SIZE_ENDCAP/N_THETA_ENDCAP;
    int     N_PHI_BARREL_CRYSTAL=floor(2*M_PI*Rin/(PHI_SEGMENTS*nomfw));
    double  D_PHI_BARREL_CRYSTAL=D_PHI_GLOBAL/N_PHI_BARREL_CRYSTAL;

    // timing
    if (system == 6) {
        double thC_end         =THETA_SIZE_ENDCAP+D_THETA_BARREL/2;
        double r0slice_end     =Rin/sin(thC_end);
        double y0slice_end     =r0slice_end*tan(D_THETA_BARREL/2.);
        double slice_front_jut =y0slice_end*sin(M_PI/2-thC_end);
        double z1slice         =Rin -slice_front_jut;
        // double z2slice         =Rin +Fdz +Rdz +slice_front_jut;
        double y1slice         =z1slice*tan(M_PI/2-THETA_SIZE_ENDCAP) +PROJECTIVE_GAP;
        double phiTiming       =nPhi_in*D_PHI_GLOBAL;
        double rT              =z1slice -2*nomth;
        double wT              =rT *tan(D_PHI_GLOBAL/2);
        int    nTiles          =ceil(y1slice/wT);
        double lT              =2*y1slice/nTiles;
        int    nCy             =floor(lT/nomth);
        double actY            =lT/nCy;
        double actX            =2*wT/nCy; 
        double rTimingAssembly =rT+nomth;
        int    nTile           =int(nEta_in/nCy);
        ROOT::Math::XYZVector dispTimingAssembly(rTimingAssembly*cos(0),rTimingAssembly*sin(0),0);
        ROOT::Math::XYZVector dispTileAssembly(0,0,-y1slice +nTile*lT +lT/2);
        int nC=nEta_in-nTile*nCy;
        int phiTimingSign=nPhi_in%2==0? 1:-1;
        int sign         =nTile%2==0? 1:-1;
        ROOT::Math::RotationZ rotZ(phiTiming);
        ROOT::Math::XYZVector dispLg(sign*phiTimingSign*(nomth/2),-wT +actX/2 + nC*actX,0);
        ROOT::Math::XYZVector dispTr(sign*phiTimingSign*(-nomth/2),0,-lT/2 +actY/2 +nC*actY);
        ROOT::Math::XYZVector dispSipmLg(0, 0, lT/2-sipmth/2);
        ROOT::Math::XYZVector dispSipmTr(0, wT-sipmth/2, 0);
        if (nDepth_in==3) {return rotZ*(dispTimingAssembly +dispTileAssembly +dispLg);}
        else if (nDepth_in==4) {return rotZ*(dispTimingAssembly +dispTileAssembly +dispLg +dispSipmLg);}
        else if (nDepth_in==5) {return rotZ*(dispTimingAssembly +dispTileAssembly +dispLg -dispSipmLg);}
        else if (nDepth_in==6) {return rotZ*(dispTimingAssembly +dispTileAssembly +dispTr);}
        else if (nDepth_in==7) {return rotZ*(dispTimingAssembly +dispTileAssembly +dispTr +dispSipmTr);}
        else if (nDepth_in==8) {return rotZ*(dispTimingAssembly +dispTileAssembly +dispTr -dispSipmTr);}
    }
    // endcap
    else if (system == 5) {
        int nTheta;
        if (nEta_in<N_THETA_ENDCAP) {nTheta=nEta_in;}
        else if (nEta_in>N_THETA_ENDCAP+N_THETA_BARREL) {
            nTheta = N_THETA_ENDCAP+N_THETA_BARREL+N_THETA_ENDCAP -nEta_in;
        }
        double thC               =D_THETA_ENDCAP/2+nTheta*D_THETA_ENDCAP;
        double RinEndcap         =EBz*tan(thC);
        int    nPhiEndcapCrystal =floor(2*M_PI*RinEndcap/(PHI_SEGMENTS*nomfw));
        double dPhiEndcapCrystal =D_PHI_GLOBAL/nPhiEndcapCrystal;
        double r0e               =RinEndcap/sin(thC);
        double r1e               =r0e+Fdz;
        int    nPhi              =int(nPhi_in/nPhiEndcapCrystal);
        int    nGamma            =nPhi_in%nPhiEndcapCrystal;
        double phi               =nPhi*D_PHI_GLOBAL;
        double gamma             =-D_PHI_GLOBAL/2+dPhiEndcapCrystal/2+dPhiEndcapCrystal*nGamma;
        int    mirror            =nEta_in<N_THETA_ENDCAP? 0:1;
        double rSlice            =RinEndcap+(Fdz+Rdz)/2;
        Position dispSlice(rSlice*cos(phi),rSlice*sin(phi),0);
        ROOT::Math::RotationZ rotZ(phi);
        ROOT::Math::RotationY rotY(M_PI*mirror);
        if (nDepth_in==1) {
            double rF=r0e+Fdz/2.;
            ROOT::Math::XYZVector dispF(rF*sin(thC)-rSlice,rF*sin(thC)*tan(gamma),rF*cos(thC)+PROJECTIVE_GAP);
            return rotY*(dispSlice+rotZ*dispF);
        }
        else if (nDepth_in==2) {
            double rR=r1e+Rdz/2.;
            ROOT::Math::XYZVector dispR(rR*sin(thC)-rSlice,rR*sin(thC)*tan(gamma),rR*cos(thC)+PROJECTIVE_GAP);
            return rotY*(dispSlice+rotZ*dispR);
        }
    }
    // barrel
    else if (system == 4) {
        int    nTheta         =nEta_in-N_THETA_ENDCAP;
        int    nPhi           =int(nPhi_in/N_PHI_BARREL_CRYSTAL);
        int    nGamma         =nPhi_in%N_PHI_BARREL_CRYSTAL;
        double phi            =nPhi*D_PHI_GLOBAL;
        double gamma          =-D_PHI_GLOBAL/2+D_PHI_BARREL_CRYSTAL/2+D_PHI_BARREL_CRYSTAL*nGamma;
        double thC            =THETA_SIZE_ENDCAP+D_THETA_BARREL/2+(nTheta*D_THETA_BARREL);
        int    projective_sign=cos(thC)>0? -1:1;
        double r0e            =Rin/sin(thC);
        double r1e            =r0e+Fdz;
        double rSlice         =Rin+(Fdz+Rdz)/2;
        Position dispSlice(rSlice*cos(phi),rSlice*sin(phi),0);
        ROOT::Math::RotationZ rotZ(phi);
        if (nDepth_in==1) {
            double rF=r0e+Fdz/2.;
            ROOT::Math::XYZVector dispF(rF*sin(thC)-rSlice,rF*sin(thC)*tan(gamma),rF*cos(thC)-projective_sign*PROJECTIVE_GAP);
            return (dispSlice+rotZ*dispF);
        }
        else if (nDepth_in==2) {
            double rR=r1e+Rdz/2.;
            ROOT::Math::XYZVector dispR(rR*sin(thC)-rSlice,rR*sin(thC)*tan(gamma),rR*cos(thC)-projective_sign*PROJECTIVE_GAP);
            return (dispSlice+rotZ*dispR);
        }
    }
    // projective fill
    else if (system == 7) {
        int nTheta   =nEta_in;
        int nPhi     =int(nPhi_in/N_PHI_BARREL_CRYSTAL);
        int nGamma   =nPhi_in%N_PHI_BARREL_CRYSTAL;
        double phi   =nPhi*D_PHI_GLOBAL;
        double gamma =-D_PHI_GLOBAL/2+D_PHI_BARREL_CRYSTAL/2+D_PHI_BARREL_CRYSTAL*nGamma;
        double thC   =M_PI/2;
        double r0e   =Rin/sin(thC);
        double r1e   =r0e+Fdz;
        double rSlice=Rin +(Fdz+Rdz)/2;
        Position dispSlice(rSlice*cos(phi),rSlice*sin(phi),0);
        ROOT::Math::RotationZ rotZ(phi);
        if (nDepth_in==1) {
            double rF=r0e+Fdz/2.;
            ROOT::Math::XYZVector dispF(rF*sin(thC)-rSlice,rF*sin(thC)*tan(gamma),rF*cos(thC)-nomfw*(N_PROJECTIVE_FILL-1)/2+nTheta*nomfw);
            return (dispSlice+rotZ*dispF);
        }
        else if (nDepth_in==2) {
            double rR=r1e+Rdz/2.;
            ROOT::Math::XYZVector dispR(rR*sin(thC)-rSlice,rR*sin(thC)*tan(gamma),rR*cos(thC)-nomfw*(N_PROJECTIVE_FILL-1)/2+nTheta*nomfw);
            return (dispSlice+rotZ*dispR);
        }
    }
    return Vector3D(0,0,0);
}

CellID SCEPCalSegmentation_k4geo::cellID(const Vector3D& /*localPosition*/, 
                                   const Vector3D& /*globalPosition*/, 
                                   const VolumeID& vID) const {
    return setCellID(System(vID), Eta(vID), Phi(vID), Depth(vID) );
}

VolumeID SCEPCalSegmentation_k4geo::setVolumeID(int System, int Eta, int Phi, int Depth) const {
    VolumeID SystemId = static_cast<VolumeID>(System);
    VolumeID EtaId = static_cast<VolumeID>(Eta);
    VolumeID PhiId = static_cast<VolumeID>(Phi);
    VolumeID DepthId = static_cast<VolumeID>(Depth);
    VolumeID vID = 0;
    _decoder->set(vID, fSystemId, SystemId);
    _decoder->set(vID, fEtaId, EtaId);
    _decoder->set(vID, fPhiId, PhiId);
    _decoder->set(vID, fDepthId, DepthId);
    return vID;
}

CellID SCEPCalSegmentation_k4geo::setCellID(int System, int Eta, int Phi, int Depth) const {
    VolumeID SystemId = static_cast<VolumeID>(System);
    VolumeID EtaId = static_cast<VolumeID>(Eta);
    VolumeID PhiId = static_cast<VolumeID>(Phi);
    VolumeID DepthId = static_cast<VolumeID>(Depth);
    VolumeID vID = 0;
    _decoder->set(vID, fSystemId, SystemId);
    _decoder->set(vID, fEtaId, EtaId);
    _decoder->set(vID, fPhiId, PhiId);
    _decoder->set(vID, fDepthId, DepthId);
    return vID;
}

int SCEPCalSegmentation_k4geo::System(const CellID& aCellID) const {
    VolumeID System = static_cast<VolumeID>(_decoder->get(aCellID, fSystemId));
    return static_cast<int>(System);
}

int SCEPCalSegmentation_k4geo::Eta(const CellID& aCellID) const {
    VolumeID Eta = static_cast<VolumeID>(_decoder->get(aCellID, fEtaId));
    return static_cast<int>(Eta);
}

int SCEPCalSegmentation_k4geo::Phi(const CellID& aCellID) const {
    VolumeID Phi = static_cast<VolumeID>(_decoder->get(aCellID, fPhiId));
    return static_cast<int>(Phi);
}

int SCEPCalSegmentation_k4geo::Depth(const CellID& aCellID) const {
    VolumeID Depth = static_cast<VolumeID>(_decoder->get(aCellID, fDepthId));
    return static_cast<int>(Depth);
}

int SCEPCalSegmentation_k4geo::getLast32bits(const CellID& aCellID) const {
    CellID aId64 = aCellID >> sizeof(int)*CHAR_BIT;
    int aId32 = (int)aId64;
    return aId32;
}

CellID SCEPCalSegmentation_k4geo::convertLast32to64(const int aId32) const {
    CellID aId64 = (CellID)aId32;
    aId64 <<= sizeof(int)*CHAR_BIT;
    return aId64;
}

}
}
