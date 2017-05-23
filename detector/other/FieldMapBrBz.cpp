#include "FieldMapBrBz.h"


#include "DD4hep/Handle.inl"
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
  type = DD4hep::Geometry::CartesianField::MAGNETIC;
} //ctor



int FieldMapBrBz::getGlobalIndex(double r, double z)
{

  int rBin = int((r - rhoMin)/rhoStep);
  int zBin = int((z - zMin  )/zStep);
  
  int GlobalIndex = -1;
  if(CoorsOrder == 1)      GlobalIndex = rBin + zBin*nRho;  //RZ coordinates ordering
  else if(CoorsOrder == 2) GlobalIndex = zBin + rBin*nZ;    //ZR coordinates ordering
  
  if(GlobalIndex < 0 || GlobalIndex > nRho*nZ-1) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: GlobalIndex = " << GlobalIndex << " is out of range (0," << nRho*nZ << ")";
    throw std::runtime_error( error.str() );
  }

  return  GlobalIndex;

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

  if (not (r <= rhoMax && z <= zMax ) ) {

    return;
  }

  int rBin,zBin;
  double r0,z0;
  double r1,z1;

  rBin = int((r - rhoMin)/rhoStep);
  zBin = int((z - zMin)/zStep);
  r0   = rhoMin + rBin*rhoStep;
  z0   = zMin   + zBin*zStep;

  if(r0 > r) r0 -= rhoStep;
  r1 = r0 + rhoStep;
  if(z0 > z) z0 -= zStep;
  z1 = z0 + zStep;

  double rd = (r - r0)/(r1 - r0);
  double zd = (z - z0)/(z1 - z0);

  const FieldMapBrBz::FieldValues_t& B_r0z0 = fieldMap[getGlobalIndex(r0,z0)];
  const FieldMapBrBz::FieldValues_t& B_r1z0 = fieldMap[getGlobalIndex(r1,z0)];
  const FieldMapBrBz::FieldValues_t& B_r0z1 = fieldMap[getGlobalIndex(r0,z1)];
  const FieldMapBrBz::FieldValues_t& B_r1z1 = fieldMap[getGlobalIndex(r1,z1)];

  double field[2] = {0.0, 0.0};

  field[0] = (1.0 - rd)*(1.0 - zd)*B_r0z0.Br + 
                     rd*(1.0 - zd)*B_r1z0.Br + 
             (1.0 - rd)*        zd*B_r0z1.Br + 
                     rd*        zd*B_r1z1.Br;

  field[1] = (1.0 - rd)*(1.0 - zd)*B_r0z0.Bz + 
                     rd*(1.0 - zd)*B_r1z0.Bz + 
             (1.0 - rd)*        zd*B_r0z1.Bz + 
                     rd*        zd*B_r1z1.Bz;

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
  std::cout << "Ntuple name:   " << NtupleName << std::endl;
  std::cout << "rho  Var name: " << rhoVar    << std::endl;
  std::cout << "z    Var name: " << zVar      << std::endl;
  std::cout << "Brho Var name: " << BrhoVar   << std::endl;
  std::cout << "Bz   Var name: " << BzVar     << std::endl;
  std::cout << std::endl;

  TTree *tree;
  file->GetObject(NtupleName.c_str(), tree);
  if (not tree) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: Tree " << NtupleName << " not found in file: " << filename;
    throw std::runtime_error( error.str() );
  }

  float r, z, Br, Bz;
  checkBranch( tree->SetBranchAddress(rhoVar.c_str(),  &r) );
  checkBranch( tree->SetBranchAddress(zVar.c_str(),    &z) );
  checkBranch( tree->SetBranchAddress(BrhoVar.c_str(), &Br));
  checkBranch( tree->SetBranchAddress(BzVar.c_str(),   &Bz));

  zStep      = -1;
  rhoStep    = -1;
  StrCoorsOrder = std::string("");
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
      StrCoorsOrder += std::string("R");
    }
    if(z != zMin && zStep < 0.0) {
      zStep         = TMath::Abs(zMin - z);
      StrCoorsOrder += std::string("Z");
    }

    //fieldMap.push_back( FieldMapBrBz::FieldValues_t( double(Br)*bScale*dd4hep::tesla,
	//					     double(Bz)*bScale*dd4hep::tesla ) );
    fieldMap.push_back( FieldMapBrBz::FieldValues_t( double(Br)*bScale*BfieldUnits,
                                                     double(Bz)*bScale*BfieldUnits ) );

  }

  if(StrCoorsOrder == std::string("RZ"))      CoorsOrder = 1;
  else if(StrCoorsOrder == std::string("ZR")) CoorsOrder = 2;

  if(rhoStep < 0) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: All rho coordinates in n-tuple have the same value!!!";
    throw std::runtime_error( error.str() );
  }
  if(rhoMax <= rhoMin) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: rho coordinates in n-tuple are not ordered from lowest to highest!!!";
    throw std::runtime_error( error.str() );
  }
  if(zStep < 0) { 
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: All z coordinates in n-tuple have the same value!!!";
    throw std::runtime_error( error.str() );
  }
  if(zMax <= zMin) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: z coordinates in n-tuple are not ordered from lowest to highest!!!";
    throw std::runtime_error( error.str() );
  }

  nRho = int((rhoMax - rhoMin)/rhoStep) + 1;
  nZ   = int((zMax   - zMin  )/zStep)   + 1;

  const int elements = nRho*nZ;
  if ( elements != treeEntries ) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: Tree does not have the expected number of entries "
          << "nRho*nZ == "  << elements
          << "  tree entries == " << treeEntries;
    throw std::runtime_error( error.str() );
  }

  rhoMin  *= coorUnits;
  rhoMax  *= coorUnits;
  rhoStep *= coorUnits;
  zMin    *= coorUnits;
  zMax    *= coorUnits;
  zStep   *= coorUnits;

  file->Close();
  delete file;

}

static DD4hep::Geometry::Ref_t create_FieldMap_rzBrBz(DD4hep::Geometry::LCDD& ,
						      DD4hep::XML::Handle_t handle ) {
  DD4hep::XML::Component xmlParameter(handle);
  bool hasFilename = xmlParameter.hasAttr(_Unicode(filename));

  if (!hasFilename) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: For a FieldMap field at least the filename xml attribute MUST be set.";
    throw std::runtime_error(error.str());
  }
  std::string  filename   = xmlParameter.attr< std::string >(_Unicode(filename));
  std::string  NtupleName = xmlParameter.attr< std::string >(_Unicode(treeName));
  std::string  rhoVar     = xmlParameter.attr< std::string >(_Unicode(rhoVarName));
  std::string  zVar       = xmlParameter.attr< std::string >(_Unicode(zVarName));
  std::string  BrhoVar    = xmlParameter.attr< std::string >(_Unicode(BrhoVarName));
  std::string  BzVar      = xmlParameter.attr< std::string >(_Unicode(BzVarName));

  double rScale      = xmlParameter.attr< double >(_Unicode(rScale));
  double zScale      = xmlParameter.attr< double >(_Unicode(zScale));
  double bScale      = xmlParameter.attr< double >(_Unicode(bScale));

  double coorUnits   = xmlParameter.attr< double >(_Unicode(coorUnits));
  double BfieldUnits = xmlParameter.attr< double >(_Unicode(BfieldUnits));

  //double rhoMin = xmlParameter.attr< double >(_Unicode(rhoMin));
  //double rhoMax = xmlParameter.attr< double >(_Unicode(rhoMax));
  //double nRho   = xmlParameter.attr< int >   (_Unicode(nRho));
  //double zMin   = xmlParameter.attr< double >(_Unicode(zMin));
  //double zMax   = xmlParameter.attr< double >(_Unicode(zMax));
  //double nZ     = xmlParameter.attr< int >   (_Unicode(nZ));

  DD4hep::Geometry::CartesianField obj;
  FieldMapBrBz* ptr = new FieldMapBrBz();
  ptr->rScale     = rScale;
  ptr->zScale     = zScale;
  ptr->bScale     = bScale;
  ptr->NtupleName = NtupleName;
  ptr->rhoVar     = rhoVar;
  ptr->zVar       = zVar;
  ptr->BrhoVar    = BrhoVar;
  ptr->BzVar      = BzVar;

  //Read the entries form the file in this place
  ptr->fillFieldMapFromTree(filename, coorUnits, BfieldUnits);

  std::cout << "rScale      " << std::setw(13) << ptr->rScale                           << std::endl;
  std::cout << "zScale      " << std::setw(13) << ptr->zScale                           << std::endl;
  std::cout << "bScale      " << std::setw(13) << ptr->bScale                           << std::endl;
  std::cout << "rhoMin      " << std::setw(13) << ptr->rhoMin/dd4hep::cm << " cm"       << std::endl;
  std::cout << "rhoMax      " << std::setw(13) << ptr->rhoMax/dd4hep::cm << " cm"       << std::endl;
  std::cout << "rhoStep     " << std::setw(13) << ptr->rhoStep/dd4hep::cm << " cm"      << std::endl;
  std::cout << "nRho        " << std::setw(13) << ptr->nRho                             << std::endl;
  std::cout << "zMin        " << std::setw(13) << ptr->zMin/dd4hep::cm << " cm"         << std::endl;
  std::cout << "zMax        " << std::setw(13) << ptr->zMax/dd4hep::cm << " cm"         << std::endl;
  std::cout << "zStep       " << std::setw(13) << ptr->zStep/dd4hep::cm << " cm"        << std::endl;
  std::cout << "nZ          " << std::setw(13) << ptr->nZ                               << std::endl;
  std::cout << "CoorsOrder  " << std::setw(13) << ptr->StrCoorsOrder.c_str()            << std::endl;
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

