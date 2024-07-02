#ifndef DRconstructor_h
#define DRconstructor_h 1

#include "detectorSegmentations/DRparamBarrel_k4geo.h"
#include "detectorSegmentations/GridDRcaloHandle_k4geo.h"

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Detector.h"

namespace ddDRcalo {
  class DRconstructor {
  public:
    DRconstructor(xml_det_t& x_det);
    ~DRconstructor() {}

    void setExpHall(dd4hep::Assembly* experimentalHall) { fExperimentalHall = experimentalHall; }
    void setDRparamBarrel(dd4hep::DDSegmentation::DRparamBarrel_k4geo* paramBarrel) { fParamBarrel = paramBarrel; }
    void setDRparamEndcap(dd4hep::DDSegmentation::DRparamEndcap_k4geo* paramEndcap) { fParamEndcap = paramEndcap; }
    void setDescription(dd4hep::Detector* description) { fDescription = description; }
    void setDetElement(dd4hep::DetElement* drDet) { fDetElement = drDet; }
    void setSipmSurf(dd4hep::OpticalSurface* sipmSurf) { fSipmSurf = sipmSurf; }
    void setSensDet(dd4hep::SensitiveDetector* sensDet) {
      fSensDet = sensDet;
      fSegmentation = dynamic_cast<dd4hep::DDSegmentation::GridDRcalo_k4geo*>( sensDet->readout().segmentation().segmentation() );
    }

    void construct();

  private:
    void initiateFibers();
    void implementTowers(xml_comp_t& x_theta, dd4hep::DDSegmentation::DRparamBase_k4geo* param, dd4hep::Volume& AssemblyBoxVol);
    void placeAssembly(xml_comp_t& x_theta, xml_comp_t& x_wafer, dd4hep::DDSegmentation::DRparamBase_k4geo* param,
                       dd4hep::Volume& AssemblyBoxVol, dd4hep::Volume& towerVol, dd4hep::Volume& sipmWaferVol,
                       int towerNo, int nPhi, bool isRHS=true);
    void implementFibers(xml_comp_t& x_theta, dd4hep::Volume& towerVol, dd4hep::Trap& trap, dd4hep::DDSegmentation::DRparamBase_k4geo* param, int towerNo);
    void implementFiber(dd4hep::Volume& towerVol, dd4hep::Trap& trap, dd4hep::Position pos, int col, int row, float fiberLen = 200.);
    double calculateDistAtZ(TGeoTrap* rootTrap, dd4hep::Position& pos, double* norm, double z);
    float calculateFiberLen(TGeoTrap* rootTrap, dd4hep::Position& pos, double* norm, double z1, double diff, double towerHeight);
    dd4hep::Box calculateFullBox(TGeoTrap* rootTrap, int& rmin, int& rmax, int& cmin, int& cmax, double dz);
    bool checkContained(TGeoTrap* rootTrap, dd4hep::Position& pos, double z, bool throwExcept=false);
    void getNormals(TGeoTrap* rootTrap, int numxBl2, double z, double* norm1, double* norm2, double* norm3, double* norm4);
    void placeUnitBox(dd4hep::Volume& fullBox, dd4hep::Volume& unitBox, int rmin, int rmax, int cmin, int cmax, bool& isEvenRow, bool& isEvenCol);

    xml_det_t fX_det;
    xml_comp_t fX_barrel;
    xml_comp_t fX_endcap;
    xml_comp_t fX_sipmDim;
    xml_comp_t fX_struct;
    xml_comp_t fX_dim;
    xml_comp_t fX_cladC;
    xml_comp_t fX_coreC;
    xml_comp_t fX_coreS;
    xml_comp_t fX_worldTube;
    xml_comp_t fX_barrelTube;
    dd4hep::Assembly* fExperimentalHall;
    dd4hep::Detector* fDescription;
    dd4hep::DDSegmentation::DRparamBarrel_k4geo* fParamBarrel;
    dd4hep::DDSegmentation::DRparamEndcap_k4geo* fParamEndcap;
    dd4hep::DetElement* fDetElement;
    dd4hep::SensitiveDetector* fSensDet;
    dd4hep::OpticalSurface* fSipmSurf;
    dd4hep::DDSegmentation::GridDRcalo_k4geo* fSegmentation;

    bool fVis;
    int fNumx, fNumy;
    std::vector< std::pair<int,int> > fFiberCoords;

    std::vector< dd4hep::Tube > fFiberEnvVec;
    std::vector< dd4hep::Tube > fFiberCoreCVec;
    std::vector< dd4hep::Tube > fFiberCoreSVec;
  };
}

#endif
