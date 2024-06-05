#include "DRconstructor.h"

ddDRcalo::DRconstructor::DRconstructor(xml_det_t& x_det)
: fX_det(x_det),
  // no default initializer for xml_comp_t
  fX_barrel( x_det.child( _Unicode(barrel) ) ),
  fX_endcap( x_det.child( _Unicode(endcap) ) ),
  fX_sipmDim( x_det.child( _Unicode(sipmDim) ) ),
  fX_struct( x_det.child( _Unicode(structure) ) ),
  fX_dim( fX_struct.child( _Unicode(dim) ) ),
  fX_cladC( fX_struct.child( _Unicode(cladC) ) ),
  fX_coreC( fX_struct.child( _Unicode(coreC) ) ),
  fX_coreS( fX_struct.child( _Unicode(coreS) ) ),
  fX_hole( fX_struct.child( _Unicode(hole) ) ),
  fX_dark( fX_struct.child( _Unicode(dark) ) ),
  fX_mirror( fX_struct.child( _Unicode(mirror) ) ) {
  fExperimentalHall = nullptr;
  fParamBarrel = nullptr;
  fDescription = nullptr;
  fDetElement = nullptr;
  fSensDet = nullptr;
  fSipmSurf = nullptr;
  fMirrorSurf = nullptr;
  fSegmentation = nullptr;
  fVis = false;
  fNumx = 0;
  fNumy = 0;
  fFiberCoords.reserve(100000);
}

void ddDRcalo::DRconstructor::construct() {
  // set vis on/off
  fVis = fDescription->visAttributes(fX_det.visStr()).showDaughters();

  implementTowers(fX_barrel, fParamBarrel);
  implementTowers(fX_endcap, fParamEndcap);
}

void ddDRcalo::DRconstructor::implementTowers(xml_comp_t& x_theta, dd4hep::DDSegmentation::DRparamBase_k4geo* param) {
  double currentTheta = x_theta.theta();
  int towerNo = x_theta.start();
  for (xml_coll_t x_dThetaColl(x_theta,_U(deltatheta)); x_dThetaColl; ++x_dThetaColl, ++towerNo ) {
    xml_comp_t x_deltaTheta = x_dThetaColl;

    std::cout << "test Tower No : " << towerNo << std::endl;

    // always use RHS for the reference
    param->SetIsRHS(true);
    param->SetDeltaTheta(x_deltaTheta.deltatheta());

    double currentToC = currentTheta + x_deltaTheta.deltatheta()/2.;
    currentTheta += x_deltaTheta.deltatheta();
    param->SetThetaOfCenter(currentToC);
    param->init();

    dd4hep::Trap assemblyEnvelop( (x_theta.height()+param->GetSipmHeight())/2., 0., 0., param->GetH1(), param->GetBl1(), param->GetTl1(), 0.,
                                  param->GetH2sipm(), param->GetBl2sipm(), param->GetTl2sipm(), 0. );

    dd4hep::Trap tower( x_theta.height()/2., 0., 0., param->GetH1(), param->GetBl1(), param->GetTl1(), 0.,
                        param->GetH2(), param->GetBl2(), param->GetTl2(), 0. );

    dd4hep::Volume towerVol( "tower", tower, fDescription->material(x_theta.materialStr()) );
    towerVol.setVisAttributes(*fDescription, x_theta.visStr());

    implementFibers(x_theta, towerVol, tower, param, towerNo);

    xml_comp_t x_wafer ( fX_sipmDim.child( _Unicode(sipmWafer) ) );

    // Photosensitive wafer
    dd4hep::Trap sipmWaferBox( x_wafer.height()/2., 0., 0., param->GetH2(), param->GetBl2(), param->GetTl2(), 0.,
                               param->GetH2(), param->GetBl2(), param->GetTl2(), 0. );
    dd4hep::Volume sipmWaferVol( "sipmWafer", sipmWaferBox, fDescription->material(x_wafer.materialStr()) );
    if (fVis) sipmWaferVol.setVisAttributes(*fDescription, x_wafer.visStr());
    dd4hep::SkinSurface(*fDescription, *fDetElement, "SiPMSurf_Tower"+std::to_string(towerNo), *fSipmSurf, sipmWaferVol);

    if (x_wafer.isSensitive()) {
      sipmWaferVol.setSensitiveDetector(*fSensDet);
    }

    // Remove sipmLayer, clear fFiberCoords instead of implementSipms()
    fFiberCoords.clear();

    for (int nPhi = 0; nPhi < x_theta.nphi(); nPhi++) {
      placeAssembly(x_theta,x_wafer,param,assemblyEnvelop,towerVol,sipmWaferVol,towerNo,nPhi);

      if ( fX_det.reflect() )
        placeAssembly(x_theta,x_wafer,param,assemblyEnvelop,towerVol,sipmWaferVol,towerNo,nPhi,false);
    }
  }

  param->filled();
  param->SetTotTowerNum( towerNo - x_theta.start() );
}

void ddDRcalo::DRconstructor::placeAssembly(xml_comp_t& x_theta, xml_comp_t& x_wafer, dd4hep::DDSegmentation::DRparamBase_k4geo* param,
                                            dd4hep::Trap& assemblyEnvelop, dd4hep::Volume& towerVol, dd4hep::Volume& sipmWaferVol,
                                            int towerNo, int nPhi, bool isRHS) {
  param->SetIsRHS(isRHS);
  int towerNoLR = param->signedTowerNo(towerNo);
  auto towerId64 = fSegmentation->setVolumeID( towerNoLR, nPhi );
  int towerId32 = fSegmentation->getFirst32bits(towerId64);

  // copy number of assemblyVolume is unpredictable, use dummy volume to make use of copy number of afterwards
  dd4hep::Volume assemblyEnvelopVol( std::string("assembly") + (isRHS ? "" : "_refl") , assemblyEnvelop, fDescription->material("Vacuum") );
  fExperimentalHall->placeVolume( assemblyEnvelopVol, param->GetAssembleTransform3D(nPhi) );

  assemblyEnvelopVol.placeVolume( towerVol, towerId32, dd4hep::Position(0.,0.,-param->GetSipmHeight()/2.) );

  // Remove sipmLayer

  dd4hep::PlacedVolume sipmWaferPhys = assemblyEnvelopVol.placeVolume( sipmWaferVol, towerId32, dd4hep::Position(0.,0.,(x_theta.height()+param->GetSipmHeight()-x_wafer.height())/2.) );
  sipmWaferPhys.addPhysVolID("eta", towerNoLR);
  sipmWaferPhys.addPhysVolID("phi", nPhi);
  sipmWaferPhys.addPhysVolID("module", 0);

  return;
}

void ddDRcalo::DRconstructor::implementFibers(xml_comp_t& x_theta, dd4hep::Volume& towerVol, dd4hep::Trap& trap, dd4hep::DDSegmentation::DRparamBase_k4geo* param, int towerNo) {
  dd4hep::Tube fiberEnv = dd4hep::Tube(0.,fX_cladC.rmax(),x_theta.height()/2.);
  dd4hep::Tube fiber = dd4hep::Tube(0.,fX_cladC.rmax(),x_theta.height()/2.-fX_mirror.height()/2.);
  dd4hep::Tube fiberC = dd4hep::Tube(0.,fX_coreC.rmin(),x_theta.height()/2.-fX_mirror.height()/2.);
  dd4hep::Tube fiberS = dd4hep::Tube(0.,fX_coreS.rmin(),x_theta.height()/2.-fX_mirror.height()/2.);

  auto rootTrap = trap.access();

  float sipmSize = fX_dim.dx();
  float gridSize = fX_dim.distance();
  float towerHeight = x_theta.height();

  float diff = fX_cladC.rmax(); // can be arbitrary small number
  float z1 = towerHeight/2.-2*diff; // can be arbitrary number slightly smaller than towerHeight/2-diff

  fNumx = static_cast<int>( std::floor( ( param->GetTl2()*2. - sipmSize )/gridSize ) ) + 1; // in phi direction
  fNumy = static_cast<int>( std::floor( ( param->GetH2()*2. - sipmSize )/gridSize ) ) + 1; // in eta direction
  int numxBl2 = static_cast<int>( std::floor( ( param->GetBl2()*2. - sipmSize )/gridSize ) ) + 1; // only used for estimating normals

  // full length fibers
  int rmin = 0, rmax = 0, cmin = 0, cmax = 0;
  dd4hep::Box fullBox = calculateFullBox(rootTrap,rmin,rmax,cmin,cmax,rootTrap->GetDz());
  dd4hep::Volume fullBoxVol("fullBox",fullBox,fDescription->material(x_theta.materialStr()));
  fullBoxVol.setVisAttributes(*fDescription, x_theta.visStr());

  dd4hep::Box unitBox = dd4hep::Box(gridSize,gridSize,x_theta.height()/2.);
  dd4hep::Volume unitBoxVol("unitBox",unitBox,fDescription->material(x_theta.materialStr()));

  if (fVis)
    unitBoxVol.setVisAttributes(*fDescription, x_theta.visStr());

  // // Remove cap (mirror or black paint in front of the fiber)
  implementFiber(unitBoxVol, trap, dd4hep::Position(-gridSize/2.,-gridSize/2.,0.), cmin, rmin, fiberEnv, fiber, fiberC, fiberS);
  implementFiber(unitBoxVol, trap, dd4hep::Position(gridSize/2.,-gridSize/2.,0.), cmin+1, rmin, fiberEnv, fiber, fiberC, fiberS);
  implementFiber(unitBoxVol, trap, dd4hep::Position(-gridSize/2.,gridSize/2.,0.), cmin, rmin+1, fiberEnv, fiber, fiberC, fiberS);
  implementFiber(unitBoxVol, trap, dd4hep::Position(gridSize/2.,gridSize/2.,0.), cmin+1, rmin+1, fiberEnv, fiber, fiberC, fiberS);

  bool isEvenRow = false, isEvenCol = false;
  placeUnitBox(fullBoxVol,unitBoxVol,rmin,rmax,cmin,cmax,isEvenRow,isEvenCol);
  towerVol.placeVolume(fullBoxVol);

  // get normals to each side
  double norm1[3] = {0.,0.,0.}, norm2[3] = {0.,0.,0.}, norm3[3] = {0.,0.,0.}, norm4[3] = {0.,0.,0.};
  getNormals(rootTrap,numxBl2,z1,norm1,norm2,norm3,norm4);

  for (int row = 0; row < fNumy; row++) {
    for (int column = 0; column < fNumx; column++) {
      auto localPosition = fSegmentation->localPosition(fNumx,fNumy,column,row);
      dd4hep::Position pos = dd4hep::Position(localPosition);

      if ( row >= rmin && row <= rmax && column >= cmin && column <= cmax ) {
        if ( ( !isEvenRow && row==rmax ) || ( !isEvenCol && column==cmax ) ) {
          double pos_[3] = {pos.x(),pos.y(),-fullBox.z()+TGeoShape::Tolerance()};
          bool check = fullBox.access()->Contains(pos_);

          if (check) {
            implementFiber(fullBoxVol, trap, pos, column, row, fiberEnv, fiber, fiberC, fiberS);
            fFiberCoords.push_back( std::make_pair(column,row) );
          }
        }
      } else {
        // outside tower
        if (!checkContained(rootTrap,pos,z1)) continue;

        double* normX = nullptr;
        double* normY = nullptr;

        // select two closest orthogonal sides
        if (column > fNumx/2) normX = norm2;
        else normX = norm4;

        if (row > fNumy/2) normY = norm3;
        else normY = norm1;

        // compare and choose the shortest fiber length
        float cand1 = calculateFiberLen(rootTrap, pos, normX, z1, diff, towerHeight);
        float cand2 = calculateFiberLen(rootTrap, pos, normY, z1, diff, towerHeight);
        float fiberLen = std::min(cand1,cand2);

        // not enough space to place fiber
        if ( fiberLen < 0. ) continue;

        // trim fiber length in the case calculated length is longer than tower height
        if (fiberLen > towerHeight) fiberLen = towerHeight;
        float centerZ = towerHeight/2. - fiberLen/2.;

        // final check
        if ( checkContained(rootTrap,pos,towerHeight/2.-fiberLen) ) {
        dd4hep::Position centerPos( pos.x(),pos.y(),centerZ );

        dd4hep::Tube shortFiberEnv = dd4hep::Tube(0.,fX_cladC.rmax(),fiberLen/2.);
        dd4hep::Tube shortFiber = dd4hep::Tube(0.,fX_cladC.rmax(),fiberLen/2.-fX_mirror.height()/2.);
        dd4hep::Tube shortFiberC = dd4hep::Tube(0.,fX_coreC.rmin(),fiberLen/2.-fX_mirror.height()/2.);
        dd4hep::Tube shortFiberS = dd4hep::Tube(0.,fX_coreS.rmin(),fiberLen/2.-fX_mirror.height()/2.);

        implementFiber(towerVol, trap, centerPos, column, row, shortFiberEnv, shortFiber, shortFiberC, shortFiberS);
        fFiberCoords.push_back( std::make_pair(column,row) );
        }
      }
    }
  }
}

// Remove cap (mirror or black paint in front of the fiber)
void ddDRcalo::DRconstructor::implementFiber(dd4hep::Volume& towerVol, dd4hep::Trap& trap, dd4hep::Position pos, int col, int row,
                                             dd4hep::Tube& fiberEnv, dd4hep::Tube& fiber, dd4hep::Tube& fiberC, dd4hep::Tube& fiberS) {
  dd4hep::Volume fiberEnvVol("fiberEnv", fiberEnv, fDescription->material(fX_hole.materialStr()));
  towerVol.placeVolume( fiberEnvVol, pos );

  if ( fSegmentation->IsCerenkov(col,row) ) { //c fiber
    dd4hep::Volume cladVol("cladC", fiber, fDescription->material(fX_cladC.materialStr()));
    fiberEnvVol.placeVolume( cladVol, dd4hep::Position(0.,0.,fX_mirror.height()/2.) );
    if (fVis) cladVol.setVisAttributes(*fDescription, fX_cladC.visStr()); // high CPU consumption!

    dd4hep::Volume coreVol("coreC", fiberC, fDescription->material(fX_coreC.materialStr()));
    if (fVis) coreVol.setVisAttributes(*fDescription, fX_coreC.visStr());
    cladVol.placeVolume( coreVol );

    coreVol.setRegion(*fDescription, fX_det.regionStr());
    cladVol.setRegion(*fDescription, fX_det.regionStr());
  } else { // s fiber
    dd4hep::Volume cladVol("cladS", fiber, fDescription->material(fX_coreC.materialStr()));
    fiberEnvVol.placeVolume( cladVol, dd4hep::Position(0.,0.,fX_mirror.height()/2.) );
    if (fVis) cladVol.setVisAttributes(*fDescription, fX_coreC.visStr());

    dd4hep::Volume coreVol("coreS", fiberS, fDescription->material(fX_coreS.materialStr()));
    if (fVis) coreVol.setVisAttributes(*fDescription, fX_coreS.visStr());
    cladVol.placeVolume( coreVol );

    coreVol.setRegion(*fDescription, fX_det.regionStr());
    cladVol.setRegion(*fDescription, fX_det.regionStr());
  }
}

double ddDRcalo::DRconstructor::calculateDistAtZ(TGeoTrap* rootTrap, dd4hep::Position& pos, double* norm, double z) {
  double pos_[3] = {pos.x(),pos.y(),z};

  return rootTrap->DistFromInside(pos_,norm);
}

float ddDRcalo::DRconstructor::calculateFiberLen(TGeoTrap* rootTrap, dd4hep::Position& pos, double* norm, double z1, double diff, double towerHeight) {
  float z2 = z1+diff;
  float y1 = calculateDistAtZ(rootTrap,pos,norm,z1);
  float y2 = calculateDistAtZ(rootTrap,pos,norm,z2);
  float ymin = std::min(y1,y2);

  // return if the distance is smaller than fiber diameter
  if ( ymin < 2.*fX_cladC.rmax() ) return -1.;

  // find the point where the fiber reaches a side of the tower
  float slope = (y2-y1)/diff;
  float y0 = (y1*z2-y2*z1)/diff;
  float z = (fX_cladC.rmax()-y0)/slope;
  float fiberLen = towerHeight/2. - z;

  return fiberLen;
}

bool ddDRcalo::DRconstructor::checkContained(TGeoTrap* rootTrap, dd4hep::Position& pos, double z, bool throwExcept) {
  double pos_[3] = {pos.x(),pos.y(),z};
  bool check = rootTrap->Contains(pos_);

  if ( throwExcept && !check ) throw std::runtime_error("Fiber must be in the tower!");
  return check;
}

void ddDRcalo::DRconstructor::getNormals(TGeoTrap* rootTrap, int numxBl2, double z, double* norm1, double* norm2, double* norm3, double* norm4) {
  dd4hep::Position pos1 = dd4hep::Position( fSegmentation->localPosition(fNumx,fNumy,fNumx/2,0) );
  dd4hep::Position pos2 = dd4hep::Position( fSegmentation->localPosition(fNumx,fNumy,fNumx/2+numxBl2/2-1,fNumy/2) );
  dd4hep::Position pos3 = dd4hep::Position( fSegmentation->localPosition(fNumx,fNumy,fNumx/2,fNumy-1) );
  dd4hep::Position pos4 = dd4hep::Position( fSegmentation->localPosition(fNumx,fNumy,fNumx/2-numxBl2/2+1,fNumy/2) );
  double pos1_[3] = {pos1.x(),pos1.y(),z};
  double pos2_[3] = {pos2.x(),pos2.y(),z};
  double pos3_[3] = {pos3.x(),pos3.y(),z};
  double pos4_[3] = {pos4.x(),pos4.y(),z};
  double dir[3] = {0.,0.,0.};

  rootTrap->ComputeNormal(pos1_,dir,norm1);
  rootTrap->ComputeNormal(pos2_,dir,norm2);
  rootTrap->ComputeNormal(pos3_,dir,norm3);
  rootTrap->ComputeNormal(pos4_,dir,norm4);
  norm1[2] = 0.; // check horizontal distance only
  norm2[2] = 0.;
  norm3[2] = 0.;
  norm4[2] = 0.;
}

dd4hep::Box ddDRcalo::DRconstructor::calculateFullBox(TGeoTrap* rootTrap, int& rmin, int& rmax, int& cmin, int& cmax, double dz) {
  float gridSize = fX_dim.distance();
  double zmin = -rootTrap->GetDz() + TGeoShape::Tolerance();
  float xmin = 0., xmax = 0., ymin = 0., ymax = 0.;

  for (int row = 0; row < fNumy; row++) { // bottom-up
    auto localPosition = dd4hep::Position( fSegmentation->localPosition(fNumx,fNumy,fNumx/2,row) );
    auto pos = localPosition + dd4hep::Position(0.,-gridSize/2.,0.);
    if ( checkContained(rootTrap,pos,zmin) ) {
      ymin = pos.y();
      rmin = row;
      break;
    }
  }

  for (int row = fNumy-1; row !=0 ; row--) { // top-down
    auto localPosition = dd4hep::Position( fSegmentation->localPosition(fNumx,fNumy,fNumx/2,row) );
    auto pos = localPosition + dd4hep::Position(0.,gridSize/2.,0.);
    if ( checkContained(rootTrap,pos,zmin) ) {
      ymax = pos.y();
      rmax = row;
      break;
    }
  }

  for (int col = 0; col < fNumx; col++) { // left-right
    auto localPosition = dd4hep::Position( fSegmentation->localPosition(fNumx,fNumy,col,rmin) );
    auto pos = localPosition + dd4hep::Position(-gridSize/2.,-gridSize/2.,0.);
    if ( checkContained(rootTrap,pos,zmin) ) {
      xmin = pos.x();
      cmin = col;
      break;
    }
  }

  for (int col = fNumx-1; col!=0; col--) { // right-left
    auto localPosition = dd4hep::Position( fSegmentation->localPosition(fNumx,fNumy,col,rmin) );
    auto pos = localPosition + dd4hep::Position(gridSize/2.,-gridSize/2.,0.);
    if ( checkContained(rootTrap,pos,zmin) ) {
      xmax = pos.x();
      cmax = col;
      break;
    }
  }

  return dd4hep::Box( (xmax-xmin)/2., (ymax-ymin)/2., dz );
}

void ddDRcalo::DRconstructor::placeUnitBox(dd4hep::Volume& fullBox, dd4hep::Volume& unitBox, int rmin, int rmax, int cmin, int cmax, bool& isEvenRow, bool& isEvenCol) {
  for (int row = rmin; row < rmax; row+=2) {
    for (int col = cmin; col < cmax; col+=2) {
      auto pos0 = dd4hep::Position( fSegmentation->localPosition(fNumx,fNumy,col,row) );
      auto pos3 = dd4hep::Position( fSegmentation->localPosition(fNumx,fNumy,col+1,row+1) );
      fullBox.placeVolume(unitBox,(pos0+pos3)/2.);
    }
  }

  isEvenRow = (rmax-rmin+1)%2==0;
  isEvenCol = (cmax-cmin+1)%2==0;

  return;
}
