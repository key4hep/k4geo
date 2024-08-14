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
  fX_worldTube( fX_struct.child( _Unicode(worldTube) ) ),
  fX_barrelTube( fX_struct.child( _Unicode(barrelTube) ) ) {
  fExperimentalHall = nullptr;
  fParamBarrel = nullptr;
  fDescription = nullptr;
  fDetElement = nullptr;
  fSensDet = nullptr;
  fSipmSurf = nullptr;
  fSegmentation = nullptr;
  fVis = false;
  fNumx = 0;
  fNumy = 0;
  fFiberCoords.reserve(100000);
  fFiberEnvVec.reserve(4000);
  fFiberCoreCVec.reserve(4000);
  fFiberCoreSVec.reserve(4000);
}

void ddDRcalo::DRconstructor::construct() {
  // set vis on/off
  fVis = fDescription->visAttributes(fX_det.visStr()).showDaughters();

  // Tube to cover all +eta DRC geometry, Rmin = 252 - 1 mm (Rmin in endcap region), Rmax = 4500 + 1 mm, length = 4500 + 1mm
  dd4hep::Tube WorldTube(fX_worldTube.rmin(), fX_worldTube.rmax(), fX_worldTube.height()/2. );
  // Tube to be subtracted from Large tube, this is for barrel region, will be subtracted from z = 0 to z = 2500 -1 mm
  // Rmin = 252 - 1 mm mm, Rmax = 2500 - 1mm, length = 2500 -1 mm
  dd4hep::Tube BarrelTube(fX_barrelTube.rmin(), fX_barrelTube.rmax(), fX_barrelTube.height()/2.);
  // Create World assembly tube by subtracting BarrelTube from WorldTube
  dd4hep::SubtractionSolid AssemblyTube(WorldTube, BarrelTube, dd4hep::Position( 0, 0, -(fX_worldTube.height() - fX_barrelTube.height())/2. ));
  dd4hep::Volume AssemblyTubeVol("AssemblyTube",AssemblyTube, fDescription->material(fX_worldTube.materialStr()) );
  AssemblyTubeVol.setVisAttributes(*fDescription, fX_worldTube.visStr());
  dd4hep::PlacedVolume PlacedAssemblyTubeVol = fExperimentalHall->placeVolume( AssemblyTubeVol, dd4hep::Position( 0, 0, fX_worldTube.height()/2.) );
  PlacedAssemblyTubeVol.addPhysVolID("assembly", 0); // TODO : Have to update box ID encoding

  initiateFibers();

  implementTowers(fX_barrel, fParamBarrel, AssemblyTubeVol);
  implementTowers(fX_endcap, fParamEndcap, AssemblyTubeVol);

  if (fX_det.reflect()) {
    auto refl_pos = dd4hep::Transform3D( dd4hep::RotationZYX(0., 0., M_PI), dd4hep::Position( 0, 0, -(fX_worldTube.height()/2.)) );
    dd4hep::PlacedVolume PlacedAssemblyTubeVol_refl = fExperimentalHall->placeVolume( AssemblyTubeVol, 1, refl_pos);
    PlacedAssemblyTubeVol_refl.addPhysVolID("assembly", 1);
  }
}

void ddDRcalo::DRconstructor::initiateFibers() {
  for (int idx = 1; idx <= 4000; idx++) {
    double length = 0.05*dd4hep::cm * idx; // fiber with length interval of 0.5 mm

    fFiberEnvVec.push_back(   dd4hep::Tube (0., fX_cladC.rmax(), length/2. ) );
    fFiberCoreCVec.push_back( dd4hep::Tube (0., fX_coreC.rmin(), length/2.) );
    fFiberCoreSVec.push_back( dd4hep::Tube (0., fX_coreS.rmin(), length/2.) );
  }
}

void ddDRcalo::DRconstructor::implementTowers(xml_comp_t& x_theta, dd4hep::DDSegmentation::DRparamBase_k4geo* param, dd4hep::Volume& AssemblyBoxVol) {
  double currentTheta = x_theta.theta();
  int towerNo = x_theta.start();
  for (xml_coll_t x_dThetaColl(x_theta,_U(deltatheta)); x_dThetaColl; ++x_dThetaColl, ++towerNo ) {
    xml_comp_t x_deltaTheta = x_dThetaColl;

    // always use RHS for the reference
    param->SetIsRHS(true);
    param->SetDeltaTheta(x_deltaTheta.deltatheta());

    double currentToC = currentTheta + x_deltaTheta.deltatheta()/2.;
    currentTheta += x_deltaTheta.deltatheta();
    param->SetThetaOfCenter(currentToC);
    param->init();

    dd4hep::Trap tower( x_theta.height()/2., 0., 0., param->GetH1(), param->GetBl1(), param->GetTl1(), 0.,
                        param->GetH2(), param->GetBl2(), param->GetTl2(), 0. );

    dd4hep::Volume towerVol( "tower", tower, fDescription->material(x_theta.materialStr()) );
    towerVol.setVisAttributes(*fDescription, x_theta.visStr());
    
    implementFibers(x_theta, towerVol, tower, param);

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
      placeAssembly(param,AssemblyBoxVol,towerVol,sipmWaferVol,towerNo,nPhi);
    }
  }

  param->filled();
  param->SetTotTowerNum( towerNo - x_theta.start() );
}

void ddDRcalo::DRconstructor::placeAssembly(dd4hep::DDSegmentation::DRparamBase_k4geo* param, dd4hep::Volume& AssemblyBoxVol,
                                            dd4hep::Volume& towerVol, dd4hep::Volume& sipmWaferVol, int towerNo, int nPhi, bool isRHS) {
  param->SetIsRHS(isRHS);
  int towerNoLR = param->signedTowerNo(towerNo);
  auto towerId64 = fSegmentation->setVolumeID( towerNoLR, nPhi );
  int towerId32 = fSegmentation->getFirst32bits(towerId64);

  dd4hep::Position towerPos = param->GetTowerPos(nPhi) + dd4hep::Position(0, 0, -(fX_worldTube.height()/2.));
  AssemblyBoxVol.placeVolume( towerVol, towerId32, dd4hep::Transform3D( param->GetRotationZYX(nPhi), towerPos ) );

  // Remove sipmLayer
  dd4hep::Position sipmPos = param->GetSipmLayerPos(nPhi) + dd4hep::Position(0, 0, -(fX_worldTube.height()/2.));
  dd4hep::PlacedVolume sipmWaferPhys = AssemblyBoxVol.placeVolume( sipmWaferVol, towerId32, dd4hep::Transform3D( param->GetRotationZYX(nPhi), sipmPos ) );
  sipmWaferPhys.addPhysVolID("eta", towerNoLR);
  sipmWaferPhys.addPhysVolID("phi", nPhi);
  sipmWaferPhys.addPhysVolID("module", 0);

  return;
}

void ddDRcalo::DRconstructor::implementFibers(xml_comp_t& x_theta, dd4hep::Volume& towerVol, dd4hep::Trap& trap, dd4hep::DDSegmentation::DRparamBase_k4geo* param) {
  auto rootTrap = trap.access();

  float sipmSize = fX_dim.dx();
  float gridSize = fX_dim.distance();
  float towerHeight = x_theta.height();

  float diff = fX_cladC.rmax(); // can be arbitrary small number
  float z1 = towerHeight/2.-2*diff; // can be arbitrary number slightly smaller than towerHeight/2-diff

  fNumx = static_cast<int>( std::floor( ( param->GetTl2()*2. - sipmSize/2. )/gridSize ) ) + 1; // in phi direction
  fNumy = static_cast<int>( std::floor( ( param->GetH2()*2. - sipmSize/2. )/gridSize ) ) + 1; // in eta direction
  int numxBl2 = static_cast<int>( std::floor( ( param->GetBl2()*2. - sipmSize/2. )/gridSize ) ) + 1; // only used for estimating normals

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
  implementFiber(unitBoxVol, dd4hep::Position(-gridSize/2.,-gridSize/2.,0.), cmin, rmin  );
  implementFiber(unitBoxVol, dd4hep::Position(gridSize/2.,-gridSize/2.,0.), cmin+1, rmin );
  implementFiber(unitBoxVol, dd4hep::Position(-gridSize/2.,gridSize/2.,0.), cmin, rmin+1 );
  implementFiber(unitBoxVol, dd4hep::Position(gridSize/2.,gridSize/2.,0.), cmin+1, rmin+1);


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
            implementFiber(fullBoxVol, pos, column, row);
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
          implementFiber(towerVol, centerPos, column, row, fiberLen);
          fFiberCoords.push_back( std::make_pair(column,row) );
        }
      }
    }
  }
}

// Remove cap (mirror or black paint in front of the fiber)
void ddDRcalo::DRconstructor::implementFiber(dd4hep::Volume& towerVol, dd4hep::Position pos, int col, int row, float fiberLen) {
  // Don't implement fiber if the length required is shorter than 0.5 mm
  if (fiberLen < 0.05*dd4hep::cm)
    return;

  int fiberIdx = int( (float) fiberLen / (float) 0.05*dd4hep::cm ) - 1; // index of fiber in fiber vectors
  // Actual length of fiber to be implemented, quantized in 0.5 mm unit
  float approxFiberLen =  0.05*dd4hep::cm * (fiberIdx + 1);
  // Fix Z position of fiber since the length of fiber can differ in [0, 0.5) mm
  dd4hep::Position fixedPos = dd4hep::Position(pos.x(), pos.y(), pos.z() + (fiberLen - approxFiberLen)/2.);

  if ( fSegmentation->IsCerenkov(col,row) ) { //c fiber
    dd4hep::Volume cladVol("cladC", fFiberEnvVec.at(fiberIdx), fDescription->material(fX_cladC.materialStr()));
    towerVol.placeVolume( cladVol, fixedPos );

    if (fVis) cladVol.setVisAttributes(*fDescription, fX_cladC.visStr()); // high CPU consumption!

    dd4hep::Volume coreVol("coreC", fFiberCoreCVec.at(fiberIdx), fDescription->material(fX_coreC.materialStr()));
    if (fVis) coreVol.setVisAttributes(*fDescription, fX_coreC.visStr());
    cladVol.placeVolume( coreVol );

    coreVol.setRegion(*fDescription, fX_det.regionStr());
    cladVol.setRegion(*fDescription, fX_det.regionStr());
  } else { // s fiber
    dd4hep::Volume cladVol("cladS", fFiberEnvVec.at(fiberIdx), fDescription->material(fX_coreC.materialStr()));
    towerVol.placeVolume( cladVol, fixedPos );
    if (fVis) cladVol.setVisAttributes(*fDescription, fX_coreC.visStr());

    dd4hep::Volume coreVol("coreS", fFiberCoreSVec.at(fiberIdx), fDescription->material(fX_coreS.materialStr()));
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
