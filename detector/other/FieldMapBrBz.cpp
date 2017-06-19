#include "FieldMapBrBz.h"

#include <DD4hep/Version.h>
#if DD4HEP_VERSION_GE(0,24)
#include <DD4hep/detail/Handle.inl>
#else
#include <DD4hep/Handle.inl>
#endif

#include <DD4hep/FieldTypes.h>


#include <DD4hep/DetFactoryHelper.h>
#include <XML/Utilities.h>

//Thanks stackoverflow
#include <boost/algorithm/string.hpp>


#include <TFile.h>
#include <TTree.h>


#include <string>
#include <stdexcept>
#include <iostream>
#include <iomanip>

using dd4hep::CartesianField;
using dd4hep::Detector;
using dd4hep::Ref_t;

DD4HEP_INSTANTIATE_HANDLE(FieldMapBrBz);

namespace {
  void checkBranch( int retCode ) {

    if ( retCode != 0 ) {
      std::stringstream error;
      error << "FieldMap[ERROR]: Branch not correctly described ";
      throw std::runtime_error( error.str() );
    }
  }
}

FieldMapBrBz::FieldMapBrBz() {
  type = CartesianField::MAGNETIC;
} //ctor


int FieldMapBrBz::getGlobalIndex(const int rBin, const int zBin)
{

  //Global index in fieldmap array from rho and z axes indexes

  int myrBin = rBin;
  int myzBin = zBin;
  if(rhoOrdering == -1) myrBin = nRho - myrBin - 1; // recalculate rho-axis index in case of high-to-low ordering
  if(zOrdering   == -1) myzBin = nZ   - myzBin - 1; // recalculate z-axis   index in case of high-to-low ordering

  int globalIndex = -1;
  if(coorsOrder == 1)       globalIndex = myrBin + myzBin*nRho;  //RZ coordinates ordering
  else if(coorsOrder == 2)  globalIndex = myzBin + myrBin*nZ;    //ZR coordinates ordering
  
  return  globalIndex;

}

/**
    Use bileanar interpolation to calculate the field at the given position
    This uses large pieces from Mokka FieldX03
 */
void FieldMapBrBz::fieldComponents(const double* pos , double* globalField) {

  //get position coordinates in our system
  const double posRho = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
  const double posZ   = pos[2];
  const double posPhi = atan2(pos[1], pos[0]);

  //Get the rho, z and phi coordinates of the 3D point
  double r   = posRho;
  double z   = posZ;
  double phi = posPhi;

  //Get positive values to do less checks when comparing
  if( z < 0 ) {
    z   *= -1;
    phi += M_PI;
  }

  //APS: Note the mokka field map does not start at 0, so we have to assume that
  //this area is covered, or add some more parameters for the values where the
  //field is supposed to be. Now we just assume it starts at 0/0/0, outside there is no field
  if(r < rhoMin) r = rhoMin;
  if(z < zMin  ) z = zMin;

  //Do nothing if rho and z point are outside fieldmap limits
  if (not (r <= rhoMax && z <= zMax ) ) return;

  //Calculate the bins on the rho and z axis containing the (r,z) point
  int rBin,zBin;
  double r0,z0;

  rBin = int((r - rhoMin)/rhoStep);
  zBin = int((z - zMin)/zStep);
  r0   = rhoMin + rBin*rhoStep;
  z0   = zMin   + zBin*zStep;

  if(r0 > r) {
    r0   -= rhoStep;
    rBin -= 1;
  }
  if(z0 > z) {
    z0   -= zStep;
    zBin -= 1;
  }

  //Get normalized coordinate of (r,z) point in bin
  double rd = (r - r0)/rhoStep;
  double zd = (z - z0)/zStep;

  //Get the field values at the four corners of bin containing the (r,z) point
  int rBin0 = rBin;
  int rBin1 = rBin+1;
  int zBin0 = zBin;
  int zBin1 = zBin+1;
  //Protection in case the sampled coordinate is exactly at maximum value of fieldma
  if(rBin1 > nRho-1) rBin1 = rBin0;
  if(zBin1 > nZ  -1) zBin1 = zBin0;
  const FieldMapBrBz::FieldValues_t& B_r0z0 = fieldMap[rBin0 + zBin0*nRho];
  const FieldMapBrBz::FieldValues_t& B_r1z0 = fieldMap[rBin1 + zBin0*nRho];
  const FieldMapBrBz::FieldValues_t& B_r0z1 = fieldMap[rBin0 + zBin1*nRho];
  const FieldMapBrBz::FieldValues_t& B_r1z1 = fieldMap[rBin1 + zBin1*nRho];

  //field at (r,z) point is linear interpolation of fielmap values at bin corners
  double field[2] = {0.0, 0.0};

  field[0] = (1.0 - rd) * (1.0 - zd) * B_r0z0.Br + 
                    rd  * (1.0 - zd) * B_r1z0.Br + 
             (1.0 - rd) *        zd  * B_r0z1.Br + 
                    rd  *        zd  * B_r1z1.Br;

  field[1] = (1.0 - rd) * (1.0 - zd) * B_r0z0.Bz + 
                    rd  * (1.0 - zd) * B_r1z0.Bz + 
             (1.0 - rd) *        zd  * B_r0z1.Bz + 
                    rd  *        zd  * B_r1z1.Bz;

  globalField[0] += field[0] * sin( phi ) ;
  globalField[1] += field[0] * cos( phi ) ;
  globalField[2] += field[1] ;

  /*
  std::cout << std::endl;
  std::cout << "FieldMapBrBz:: " << std::endl;
  std::cout << "(x,y,z)     = (" << pos[0]/dd4hep::cm  << "," << pos[1]/dd4hep::cm  << "," << pos[2]/dd4hep::cm  << ") cm" << std::endl;
  std::cout << "B(x,y,z)    = (" << globalField[0]/dd4hep::tesla << "," << globalField[1]/dd4hep::tesla << "," << globalField[2]/dd4hep::tesla << ") tesla" << std::endl;
  std::cout << std::endl;
  */

  return;

}


void FieldMapBrBz::fillFieldMapFromTree(const std::string& filename,
                                        double coorUnits, double BfieldUnits) {

  TFile *file = TFile::Open( filename.c_str() );
  if (not file) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: File not found: " << filename;
    throw std::runtime_error( error.str() );
  }

  std::cout << std::endl;
  std::cout << "Ntuple name:   " << ntupleName << std::endl;
  std::cout << "rho  Var name: " << rhoVar    << std::endl;
  std::cout << "z    Var name: " << zVar      << std::endl;
  std::cout << "Brho Var name: " << BrhoVar   << std::endl;
  std::cout << "Bz   Var name: " << BzVar     << std::endl;
  std::cout << std::endl;

  TTree *tree;
  file->GetObject(ntupleName.c_str(), tree);
  if (not tree) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: Tree " << ntupleName << " not found in file: " << filename;
    throw std::runtime_error( error.str() );
  }

  //Set branch adresses
  float r, z, Br, Bz;
  checkBranch( tree->SetBranchAddress(rhoVar.c_str(),  &r) );
  checkBranch( tree->SetBranchAddress(zVar.c_str(),    &z) );
  checkBranch( tree->SetBranchAddress(BrhoVar.c_str(), &Br));
  checkBranch( tree->SetBranchAddress(BzVar.c_str(),   &Bz));

  //Loop over the tree entries. In this loop get,
  // - min, max and step-size values of fieldmap coordinates
  // - coordinates ordering
  zStep       = -1;
  rhoStep     = -1;
  rhoOrdering =  1;
  zOrdering   =  1;
  strCoorsOrder = std::string("");
  const int treeEntries = tree->GetEntries();
  fieldMap.reserve(treeEntries);
  for(int i=0;i<treeEntries;i++) {
    tree->GetEntry(i);

    if(i == 0) {
      rhoMin = r;
      zMin   = z;
    }
    if(i == treeEntries-1) {
      rhoMax = r;
      zMax   = z;
    }

    if(r != rhoMin && rhoStep < 0.0) {
      rhoStep     = TMath::Abs(rhoMin - r);
      strCoorsOrder += std::string("R");
    }
    if(z != zMin && zStep < 0.0) {
      zStep         = TMath::Abs(zMin - z);
      strCoorsOrder += std::string("Z");
    }
  }

  if(strCoorsOrder == std::string("RZ"))      coorsOrder = 1;
  else if(strCoorsOrder == std::string("ZR")) coorsOrder = 2;

  if(rhoStep < 0) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: All rho coordinates in n-tuple have the same value!!!";
    throw std::runtime_error( error.str() );
  }
  if(rhoMax < rhoMin) {
    //rho variable is scanned from high-to-low
    rhoOrdering = -1;
    double aux = rhoMax;
    rhoMax = rhoMin;
    rhoMin = aux;
  }
  if(zStep < 0) { 
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: All z coordinates in n-tuple have the same value!!!";
    throw std::runtime_error( error.str() );
  }
  if(zMax < zMin) {
    //z variable is scanned from high-to-low
    zOrdering = -1;
    double aux = zMax;
    zMax = zMin;
    zMin = aux;
  }

  //Calculate number of bins in fieldmap
  nRho = round(((rhoMax - rhoMin)/rhoStep) + 1);
  nZ   = round(((zMax   - zMin  )/zStep)   + 1);

  //Set coordinates parameters units
  rhoMin  *= coorUnits;
  rhoMax  *= coorUnits;
  rhoStep *= coorUnits;
  zMin    *= coorUnits;
  zMax    *= coorUnits;
  zStep   *= coorUnits;

  const int elements = nRho*nZ;
  if ( elements != treeEntries ) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: Tree does not have the expected number of entries "
          << "nRho*nZ (" << nRho << "*" << nZ << ") == "  << elements
          << "  tree entries == " << treeEntries;
    throw std::runtime_error( error.str() );
  }

  //Fill the array with the Bfield values in the RZ order
  for(int iz=0;iz<nZ;iz++) {
    for(int ir=0;ir<nRho;ir++) {
      tree->GetEntry(getGlobalIndex(ir,iz));
      fieldMap.push_back( FieldMapBrBz::FieldValues_t( double(Br)*bScale*BfieldUnits,
                                                       double(Bz)*bScale*BfieldUnits ) );
    }
  }

  file->Close();
  delete file;

}

static Ref_t create_FieldMap_rzBrBz(Detector& ,
                                    dd4hep::xml::Handle_t handle ) {
  dd4hep::xml::Component xmlParameter(handle);
  bool hasFilename = xmlParameter.hasAttr(_Unicode(filename));

  if (!hasFilename) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: For a FieldMap field at least the filename xml attribute MUST be set.";
    throw std::runtime_error(error.str());
  }

  std::string  filename   = xmlParameter.attr< std::string >(_Unicode(filename));
  std::string  ntupleName = xmlParameter.attr< std::string >(_Unicode(treeName));
  std::string  rhoVar     = xmlParameter.attr< std::string >(_Unicode(rhoVarName));
  std::string  zVar       = xmlParameter.attr< std::string >(_Unicode(zVarName));
  std::string  BrhoVar    = xmlParameter.attr< std::string >(_Unicode(BrhoVarName));
  std::string  BzVar      = xmlParameter.attr< std::string >(_Unicode(BzVarName));

  double rScale      = xmlParameter.attr< double >(_Unicode(rScale));
  double zScale      = xmlParameter.attr< double >(_Unicode(zScale));
  double bScale      = xmlParameter.attr< double >(_Unicode(bScale));

  double coorUnits   = xmlParameter.attr< double >(_Unicode(coorUnits));
  double BfieldUnits = xmlParameter.attr< double >(_Unicode(BfieldUnits));

  CartesianField obj;
  FieldMapBrBz* ptr = new FieldMapBrBz();
  ptr->rScale     = rScale;
  ptr->zScale     = zScale;
  ptr->bScale     = bScale;
  ptr->ntupleName = ntupleName;
  ptr->rhoVar     = rhoVar;
  ptr->zVar       = zVar;
  ptr->BrhoVar    = BrhoVar;
  ptr->BzVar      = BzVar;

  //Read the entries form the file in this place
  ptr->fillFieldMapFromTree(filename, coorUnits, BfieldUnits);

  std::string strRhoOrdering("low-to-high");
  std::string strZOrdering("low-to-high");
  if(ptr->rhoOrdering == -1) strRhoOrdering = "high-to-low";
  if(ptr->zOrdering   == -1) strZOrdering   = "high-to-low";

  std::cout << "rScale      " << std::setw(13) << ptr->rScale                           << std::endl;
  std::cout << "zScale      " << std::setw(13) << ptr->zScale                           << std::endl;
  std::cout << "bScale      " << std::setw(13) << ptr->bScale                           << std::endl;
  std::cout << "rhoMin      " << std::setw(13) << ptr->rhoMin/dd4hep::cm << " cm"       << std::endl;
  std::cout << "rhoMax      " << std::setw(13) << ptr->rhoMax/dd4hep::cm << " cm"       << std::endl;
  std::cout << "rhoStep     " << std::setw(13) << ptr->rhoStep/dd4hep::cm << " cm"      << std::endl;
  std::cout << "nRho        " << std::setw(13) << ptr->nRho                             << std::endl;
  std::cout << "rhoOrdering " << std::setw(13) << strRhoOrdering.c_str()                << std::endl;
  std::cout << "zMin        " << std::setw(13) << ptr->zMin/dd4hep::cm << " cm"         << std::endl;
  std::cout << "zMax        " << std::setw(13) << ptr->zMax/dd4hep::cm << " cm"         << std::endl;
  std::cout << "zStep       " << std::setw(13) << ptr->zStep/dd4hep::cm << " cm"        << std::endl;
  std::cout << "nZ          " << std::setw(13) << ptr->nZ                               << std::endl;
  std::cout << "zOrdering   " << std::setw(13) << strZOrdering.c_str()                  << std::endl;
  std::cout << "CoorsOrder  " << std::setw(13) << ptr->strCoorsOrder.c_str()            << std::endl;
  std::cout << "coorUnits   " << std::setw(13) << coorUnits/dd4hep::cm << " cm"         << std::endl;
  std::cout << "BfieldUnits " << std::setw(13) << BfieldUnits/dd4hep::tesla << " tesla" << std::endl;

  ptr->rhoMin  *= rScale;
  ptr->rhoMax  *= rScale;
  ptr->rhoStep *= rScale;
  ptr->zMin    *= zScale;
  ptr->zMax    *= zScale;
  ptr->zStep   *= zScale;

  obj.assign(ptr, xmlParameter.nameStr(), xmlParameter.typeStr());

  return obj;

}
DECLARE_XMLELEMENT(FieldBrBz,create_FieldMap_rzBrBz)

