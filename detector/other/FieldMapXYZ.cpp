#include "FieldMapXYZ.h"


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

DD4HEP_INSTANTIATE_HANDLE(FieldMapXYZ);

namespace {
  void checkBranch( int retCode ) {

    if ( retCode != 0 ) {
      std::stringstream error;
      error << "FieldMap[ERROR]: Branch not correctly described ";
      throw std::runtime_error( error.str() );
    }
  }
}

FieldMapXYZ::FieldMapXYZ() {
  type = DD4hep::Geometry::CartesianField::MAGNETIC;
} //ctor

int FieldMapXYZ::getGlobalIndex(double x, double y, double z)
{

  int xBin = int((x - xMin)/xStep);
  int yBin = int((y - yMin)/yStep);
  int zBin = int((z - zMin)/zStep);
  
  int GlobalIndex = -1;
  if(CoorsOrder == 1)       GlobalIndex = xBin + yBin*(nX) + zBin*(nX*nY); //XYZ coordinates ordering
  else if(CoorsOrder == 2)  GlobalIndex = xBin + zBin*(nX) + yBin*(nX*nZ); //XZY coordinates ordering
  else if(CoorsOrder == 3)  GlobalIndex = yBin + xBin*(nY) + zBin*(nY*nX); //YXZ coordinates ordering
  else if(CoorsOrder == 4)  GlobalIndex = yBin + zBin*(nY) + xBin*(nY*nZ); //YZX coordinates ordering
  else if(CoorsOrder == 5)  GlobalIndex = zBin + xBin*(nZ) + yBin*(nZ*nX); //ZXY coordinates ordering
  else if(CoorsOrder == 6)  GlobalIndex = zBin + yBin*(nZ) + xBin*(nZ*nY); //ZYX coordinates ordering
  
  if(GlobalIndex < 0 || GlobalIndex > nX*nY*nZ-1) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: GlobalIndex = " << GlobalIndex << " is out of range (0," << nX*nY*nZ << ")";
    throw std::runtime_error( error.str() );
  }

  return  GlobalIndex;

}

/**
    Use bileanar interpolation to calculate the field at the given position
    This uses large pieces from Mokka FieldX03
 */
void FieldMapXYZ::fieldComponents(const double* pos , double* globalField) {

  //get position coordinates in our system
  const double x = pos[0];
  const double y = pos[1];
  const double z = pos[2];

  //APS: Note the mokka field map does not start at 0, so we have to assume that
  //this area is covered, or add some more parameters for the values where the
  //field is supposed to be. Now we just assume it starts at 0/0/0, outside there is no field
  if (not ( x >= xMin && x <= xMax &&
            y >= yMin && y <= yMax &&
            z >= zMin && z <= zMax )
     ) {
    return;
  }

  int xBin,yBin,zBin;
  double x0,y0,z0;
  double x1,y1,z1;

  xBin = int((x - xMin)/xStep);
  yBin = int((y - yMin)/yStep);
  zBin = int((z - zMin)/zStep);

  x0   = xMin + xBin*xStep;
  y0   = yMin + yBin*yStep;
  z0   = zMin + zBin*zStep;

  if(x0 >= x) x0 -= xStep;
  x1 = x0 + xStep;
  if(y0 >= y) y0 -= yStep;
  y1 = y0 + yStep;
  if(z0 >= z) z0 -= zStep;
  z1 = z0 + zStep;

  double xd = (x - x0)/(x1 - x0);
  double yd = (y - y0)/(y1 - y0);
  double zd = (z - z0)/(z1 - z0);

  const FieldMapXYZ::FieldValues_t& B_x0y0z0 = fieldMap[getGlobalIndex(x0,y0,z0)];
  const FieldMapXYZ::FieldValues_t& B_x1y0z0 = fieldMap[getGlobalIndex(x1,y0,z0)];
  const FieldMapXYZ::FieldValues_t& B_x0y0z1 = fieldMap[getGlobalIndex(x0,y0,z1)];
  const FieldMapXYZ::FieldValues_t& B_x1y0z1 = fieldMap[getGlobalIndex(x1,y0,z1)];
  const FieldMapXYZ::FieldValues_t& B_x0y1z0 = fieldMap[getGlobalIndex(x0,y1,z0)];
  const FieldMapXYZ::FieldValues_t& B_x1y1z0 = fieldMap[getGlobalIndex(x1,y1,z0)];
  const FieldMapXYZ::FieldValues_t& B_x0y1z1 = fieldMap[getGlobalIndex(x0,y1,z1)];
  const FieldMapXYZ::FieldValues_t& B_x1y1z1 = fieldMap[getGlobalIndex(x1,y1,z1)];

  double B_00,B_01,B_10,B_11,B_0,B_1,B;

  //X-component
  B_00 = (1.0 - xd)*B_x0y0z0.Bx + xd*B_x1y0z0.Bx;
  B_01 = (1.0 - xd)*B_x0y0z1.Bx + xd*B_x1y0z1.Bx;
  B_10 = (1.0 - xd)*B_x0y1z0.Bx + xd*B_x1y1z0.Bx;
  B_11 = (1.0 - xd)*B_x0y1z1.Bx + xd*B_x1y1z1.Bx;
  B_0  = (1.0 - yd)*B_00        + yd*B_10;
  B_1  = (1.0 - yd)*B_01        + yd*B_11;
  B    = (1.0 - zd)*B_0         + zd*B_1;
  globalField[0] += B;

  //Y-component
  B_00 = (1.0 - xd)*B_x0y0z0.By + xd*B_x1y0z0.By;
  B_01 = (1.0 - xd)*B_x0y0z1.By + xd*B_x1y0z1.By;
  B_10 = (1.0 - xd)*B_x0y1z0.By + xd*B_x1y1z0.By;
  B_11 = (1.0 - xd)*B_x0y1z1.By + xd*B_x1y1z1.By;
  B_0  = (1.0 - yd)*B_00        + yd*B_10;
  B_1  = (1.0 - yd)*B_01        + yd*B_11;
  B    = (1.0 - zd)*B_0         + zd*B_1;
  globalField[1] += B;

  //Z-component
  B_00 = (1.0 - xd)*B_x0y0z0.Bz + xd*B_x1y0z0.Bz;
  B_01 = (1.0 - xd)*B_x0y0z1.Bz + xd*B_x1y0z1.Bz;
  B_10 = (1.0 - xd)*B_x0y1z0.Bz + xd*B_x1y1z0.Bz;
  B_11 = (1.0 - xd)*B_x0y1z1.Bz + xd*B_x1y1z1.Bz;
  B_0  = (1.0 - yd)*B_00        + yd*B_10;
  B_1  = (1.0 - yd)*B_01        + yd*B_11;
  B    = (1.0 - zd)*B_0         + zd*B_1;
  globalField[2] += B;
  
  /* 
  std::cout << std::endl;
  std::cout << "FieldMapXYZ:: " << std::endl;
  std::cout << "(x,y,z)     = (" << x/dd4hep::cm  << "," << y/dd4hep::cm  << "," << z/dd4hep::cm  << ") cm" << std::endl;
  std::cout << "(x0,y0,z0)  = (" << x0/dd4hep::cm << "," << y0/dd4hep::cm << "," << z0/dd4hep::cm << ") cm" << std::endl;
  std::cout << "(x1,y1,z1)  = (" << x1/dd4hep::cm << "," << y1/dd4hep::cm << "," << z1/dd4hep::cm << ") cm" << std::endl;
  std::cout << "B(x0,y0,z0) = (" << B_x0y0z0.Bx/dd4hep::tesla << "," << B_x0y0z0.By/dd4hep::tesla << "," << B_x0y0z0.Bz/dd4hep::tesla << ") tesla" << std::endl;
  std::cout << "B(x1,y0,z0) = (" << B_x1y0z0.Bx/dd4hep::tesla << "," << B_x1y0z0.By/dd4hep::tesla << "," << B_x1y0z0.Bz/dd4hep::tesla << ") tesla" << std::endl;
  std::cout << "B(x0,y0,z1) = (" << B_x0y0z1.Bx/dd4hep::tesla << "," << B_x0y0z1.By/dd4hep::tesla << "," << B_x0y0z1.Bz/dd4hep::tesla << ") tesla" << std::endl;
  std::cout << "B(x1,y0,z1) = (" << B_x1y0z1.Bx/dd4hep::tesla << "," << B_x1y0z1.By/dd4hep::tesla << "," << B_x1y0z1.Bz/dd4hep::tesla << ") tesla" << std::endl;
  std::cout << "B(x0,y1,z0) = (" << B_x0y1z0.Bx/dd4hep::tesla << "," << B_x0y1z0.By/dd4hep::tesla << "," << B_x0y1z0.Bz/dd4hep::tesla << ") tesla" << std::endl;
  std::cout << "B(x1,y1,z0) = (" << B_x1y1z0.Bx/dd4hep::tesla << "," << B_x1y1z0.By/dd4hep::tesla << "," << B_x1y1z0.Bz/dd4hep::tesla << ") tesla" << std::endl;
  std::cout << "B(x0,y1,z1) = (" << B_x0y1z1.Bx/dd4hep::tesla << "," << B_x0y1z1.By/dd4hep::tesla << "," << B_x0y1z1.Bz/dd4hep::tesla << ") tesla" << std::endl;
  std::cout << "B(x1,y1,z1) = (" << B_x1y1z1.Bx/dd4hep::tesla << "," << B_x1y1z1.By/dd4hep::tesla << "," << B_x1y1z1.Bz/dd4hep::tesla << ") tesla" << std::endl;
  std::cout << "B(x,y,z)    = (" << globalField[0]/dd4hep::tesla << "," << globalField[1]/dd4hep::tesla << "," << globalField[2]/dd4hep::tesla << ") tesla" << std::endl;
  std::cout << std::endl;
  */

  return;

}

void FieldMapXYZ::fillFieldMapFromTree(const std::string& filename,
                                       const std::string& treeString,
                                       double coorUnits, double BfieldUnits) {


  TFile *file = TFile::Open( filename.c_str() );
  if (not file) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: File not found: " << filename;
    throw std::runtime_error( error.str() );
  }

  std::vector<std::string> strs;
  boost::split(strs, treeString, boost::is_any_of(": "));
  for (std::vector<std::string>::iterator it = strs.begin(); it != strs.end(); ++it) {
    std::cout << "\"" << (*it).c_str() << "\"" << std::endl;
  }

  if( strs.size() != 7 ) {
    std::stringstream error;
    error << "FieldMap[ERROR]: the treeDescription " << treeString
	  << " is not complete. For example 'ntuple:x:y:z:Bx:By:Bz'";
    throw std::runtime_error( error.str() );
  }

  //Getting the tree name and branch variables names
  std::string  NtupleName("");
  std::string  xVar("");
  std::string  yVar("");
  std::string  zVar("");
  std::string  BxVar("");
  std::string  ByVar("");
  std::string  BzVar("");
  for(int i=0;i<int(strs.size());i++) {
    if(strs[i].find(std::string("Bx")) != std::string::npos) {
      BxVar = strs[i];
    }
    else if(strs[i].find(std::string("By")) != std::string::npos) {
      ByVar = strs[i];
    }
    else if(strs[i].find(std::string("Bz")) != std::string::npos) {
      BzVar = strs[i];
    }
    else if(strs[i].find(std::string("x")) != std::string::npos || strs[i].find(std::string("X")) != std::string::npos) {
      xVar = strs[i];
    }
    else if(strs[i].find(std::string("y")) != std::string::npos || strs[i].find(std::string("Y")) != std::string::npos) {
      yVar = strs[i];
    }
    else if(strs[i].find(std::string("z")) != std::string::npos || strs[i].find(std::string("Z")) != std::string::npos) {
      zVar = strs[i];
    }
    else {
      NtupleName = strs[i];
    }
  }

  if(NtupleName == std::string("")) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: ntuple name not set!";
    throw std::runtime_error( error.str() );
  }
  if(xVar == std::string("")) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: x variable name not set!";
    throw std::runtime_error( error.str() );
  }
  if(yVar == std::string("")) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: y variable name not set!";
    throw std::runtime_error( error.str() );
  }
  if(zVar == std::string("")) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: z variable name not set!";
    throw std::runtime_error( error.str() );
  }
  if(BxVar == std::string("")) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: Bx variable name not set!";
    throw std::runtime_error( error.str() );
  }
  if(ByVar == std::string("")) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: By variable name not set!";
    throw std::runtime_error( error.str() );
  }
  if(BzVar == std::string("")) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: Bz variable name not set!";
    throw std::runtime_error( error.str() );
  }

  std::cout << std::endl;
  std::cout << "Ntuple name: " << NtupleName << std::endl;
  std::cout << "x  Var name: " << xVar       << std::endl;
  std::cout << "y  Var name: " << yVar       << std::endl;
  std::cout << "z  Var name: " << zVar       << std::endl;
  std::cout << "Bx Var name: " << BxVar      << std::endl;
  std::cout << "By Var name: " << ByVar      << std::endl;
  std::cout << "Bz Var name: " << BzVar      << std::endl;
  std::cout << std::endl;

  TTree *tree;
  file->GetObject(NtupleName.c_str(), tree);
  if (not tree) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: Tree " << NtupleName << " not found in file: " << filename;
    throw std::runtime_error( error.str() );
  }

  float x, y, z, Bx, By, Bz;
  checkBranch( tree->SetBranchAddress(xVar.c_str(), &x) );
  checkBranch( tree->SetBranchAddress(yVar.c_str(), &y) );
  checkBranch( tree->SetBranchAddress(zVar.c_str(), &z) );
  checkBranch( tree->SetBranchAddress(BxVar.c_str(), &Bx) );
  checkBranch( tree->SetBranchAddress(ByVar.c_str(), &By) );
  checkBranch( tree->SetBranchAddress(BzVar.c_str(), &Bz) );

  xStep      = -1;
  yStep      = -1;
  zStep      = -1;
  StrCoorsOrder = std::string("");
  const int treeEntries = tree->GetEntries();
  fieldMap.reserve(treeEntries);
  for(int i=0;i<treeEntries;i++) {
    tree->GetEntry(i);

    if(i == 0) {
      xMin = x;
      yMin = y;
      zMin = z;
    }
    if(i == treeEntries-1) {
      xMax = x;
      yMax = y;
      zMax = z;
    }
    
    if(x != xMin && xStep < 0.0) {
      xStep       = TMath::Abs(xMin - x);
      StrCoorsOrder += std::string("X");
    }
    if(y != yMin && yStep < 0.0) {
      yStep       = TMath::Abs(yMin - y);
      StrCoorsOrder += std::string("Y");
    }
    if(z != zMin && zStep < 0.0) {
      zStep       = TMath::Abs(zMin - z);
      StrCoorsOrder += std::string("Z");
    }

    //fieldMap.push_back( FieldMapXYZ::FieldValues_t(double(Bx)*bScale*dd4hep::tesla,
    //                                               double(By)*bScale*dd4hep::tesla,
	//					   double(Bz)*bScale*dd4hep::tesla ) );
    fieldMap.push_back( FieldMapXYZ::FieldValues_t(double(Bx)*bScale*BfieldUnits,
                                                   double(By)*bScale*BfieldUnits,
                                                   double(Bz)*bScale*BfieldUnits ) );
  }

  if(StrCoorsOrder == TString("XYZ"))       CoorsOrder = 1;
  else if(StrCoorsOrder == TString("XZY"))  CoorsOrder = 2;
  else if(StrCoorsOrder == TString("YXZ"))  CoorsOrder = 3;
  else if(StrCoorsOrder == TString("YZX"))  CoorsOrder = 4;
  else if(StrCoorsOrder == TString("ZXY"))  CoorsOrder = 5;
  else if(StrCoorsOrder == TString("ZYX"))  CoorsOrder = 6;

  if(xStep < 0) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: All x coordinates in n-tuple have the same value!!!";
    throw std::runtime_error( error.str() );
  }
  if(xMax <= xMin) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: x coordinates in n-tuple are not ordered from lowest to highest!!!";
    throw std::runtime_error( error.str() );
  }
  if(yStep < 0) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: All y coordinates in n-tuple have the same value!!!";
    throw std::runtime_error( error.str() );
  }
  if(yMax <= yMin) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: y coordinates in n-tuple are not ordered from lowest to highest!!!";
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

  nX = int((xMax - xMin)/xStep) + 1;
  nY = int((yMax - yMin)/yStep) + 1;
  nZ = int((zMax - zMin)/zStep) + 1;

  const int elements = nX*nY*nZ;
  if ( elements != treeEntries ) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: Tree does not have the expected number of entries "
          << "nX*nY*nZ == "  << elements
          << "  tree entries == " << treeEntries;
    throw std::runtime_error( error.str() );
  }

  xMin  *= coorUnits;
  xMax  *= coorUnits;
  xStep *= coorUnits;
  yMin  *= coorUnits;
  yMax  *= coorUnits;
  yStep *= coorUnits;
  zMin  *= coorUnits;
  zMax  *= coorUnits;
  zStep *= coorUnits;

  file->Close();
  delete file;

}

static DD4hep::Geometry::Ref_t create_FieldMap_XYZ(DD4hep::Geometry::LCDD& ,
						   DD4hep::XML::Handle_t handle ) {
  DD4hep::XML::Component xmlParameter(handle);
  bool hasFilename = xmlParameter.hasAttr(_Unicode(filename));

  if (!hasFilename) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: For a FieldMap field at least the filename xml attribute MUST be set.";
    throw std::runtime_error(error.str());
  }
  std::string filename   = xmlParameter.attr< std::string >(_Unicode(filename));
  std::string treeString = xmlParameter.attr< std::string >(_Unicode(tree));

  double xScale = xmlParameter.attr< double >(_Unicode(xScale));
  double yScale = xmlParameter.attr< double >(_Unicode(yScale));
  double zScale = xmlParameter.attr< double >(_Unicode(zScale));
  double bScale = xmlParameter.attr< double >(_Unicode(bScale));

  double coorUnits   = xmlParameter.attr< double >(_Unicode(coorUnits));
  double BfieldUnits = xmlParameter.attr< double >(_Unicode(BfieldUnits));

  //double xMin   = xmlParameter.attr< double >(_Unicode(xMin));
  //double xMax   = xmlParameter.attr< double >(_Unicode(xMax));
  //double xStep  = xmlParameter.attr< double >(_Unicode(xStep));

  //double yMin   = xmlParameter.attr< double >(_Unicode(yMin));
  //double yMax   = xmlParameter.attr< double >(_Unicode(yMax));
  //double yStep  = xmlParameter.attr< double >(_Unicode(yStep));

  //double zMin   = xmlParameter.attr< double >(_Unicode(zMin));
  //double zMax   = xmlParameter.attr< double >(_Unicode(zMax));
  //double zStep  = xmlParameter.attr< double >(_Unicode(zStep));
 
  DD4hep::Geometry::CartesianField obj;
  FieldMapXYZ* ptr = new FieldMapXYZ();
  ptr->xScale = xScale;
  ptr->yScale = yScale;
  ptr->zScale = zScale;
  ptr->bScale = bScale;

  //Read the entries form the file in this place
  ptr->fillFieldMapFromTree(filename, treeString,coorUnits,BfieldUnits);

  std::cout << "xScale      " << std::setw(13) << ptr->xScale                           << std::endl;
  std::cout << "zScale      " << std::setw(13) << ptr->zScale                           << std::endl;
  std::cout << "zScale      " << std::setw(13) << ptr->zScale                           << std::endl;
  std::cout << "bScale      " << std::setw(13) << ptr->bScale                           << std::endl;
  std::cout << "xMin        " << std::setw(13) << ptr->xMin/dd4hep::cm << " cm"         << std::endl;
  std::cout << "xMax        " << std::setw(13) << ptr->xMax/dd4hep::cm << " cm"         << std::endl;
  std::cout << "xStep       " << std::setw(13) << ptr->xStep/dd4hep::cm << " cm"        << std::endl;
  std::cout << "nX          " << std::setw(13) << ptr->nX                               << std::endl;
  std::cout << "yMin        " << std::setw(13) << ptr->yMin/dd4hep::cm << " cm"         << std::endl;
  std::cout << "yMax        " << std::setw(13) << ptr->yMax/dd4hep::cm << " cm"         << std::endl;
  std::cout << "yStep       " << std::setw(13) << ptr->yStep/dd4hep::cm << " cm"        << std::endl;
  std::cout << "nY          " << std::setw(13) << ptr->nY                               << std::endl;
  std::cout << "zMin        " << std::setw(13) << ptr->zMin/dd4hep::cm << " cm"         << std::endl;
  std::cout << "zMax        " << std::setw(13) << ptr->zMax/dd4hep::cm << " cm"         << std::endl;
  std::cout << "zStep       " << std::setw(13) << ptr->zStep/dd4hep::cm << " cm"        << std::endl;
  std::cout << "nZ          " << std::setw(13) << ptr->nZ                               << std::endl;
  std::cout << "CoorsOrder  " << std::setw(13) << ptr->StrCoorsOrder.c_str()            << std::endl;
  std::cout << "coorUnits   " << std::setw(13) << coorUnits/dd4hep::cm << " cm"         << std::endl;
  std::cout << "BfieldUnits " << std::setw(13) << BfieldUnits/dd4hep::tesla << " tesla" << std::endl;

  ptr->xMin   *= xScale;
  ptr->xMax   *= xScale;
  ptr->xStep  *= xScale;
  ptr->yMin   *= yScale;
  ptr->yMax   *= yScale;
  ptr->yStep  *= yScale;
  ptr->zMin   *= zScale;
  ptr->zMax   *= zScale;
  ptr->zStep  *= zScale;

  obj.assign(ptr, xmlParameter.nameStr(), xmlParameter.typeStr());

  return obj;

}
DECLARE_XMLELEMENT(FieldXYZ,create_FieldMap_XYZ)

