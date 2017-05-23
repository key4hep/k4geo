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


int FieldMapXYZ::getGlobalIndex(const int xBin, const int yBin, const int zBin)
{

  //Global index in fieldmap array from x, y and z axes indexes

  int MyxBin = xBin;
  int MyyBin = yBin;
  int MyzBin = zBin;
  if(xOrdering == -1) MyxBin = nX - MyxBin - 1; // recalculate x-axis index in case of high-to-low ordering
  if(yOrdering == -1) MyyBin = nY - MyyBin - 1; // recalculete y-axis index in case of high-to-low ordering
  if(zOrdering == -1) MyzBin = nZ - MyzBin - 1; // recalculate z-axis index in case of high-to-low ordering

  int GlobalIndex = -1;
  if(CoorsOrder == 1)       GlobalIndex = MyxBin + MyyBin*(nX) + MyzBin*(nX*nY); //XYZ coordinates ordering
  else if(CoorsOrder == 2)  GlobalIndex = MyxBin + MyzBin*(nX) + MyyBin*(nX*nZ); //XZY coordinates ordering
  else if(CoorsOrder == 3)  GlobalIndex = MyyBin + MyxBin*(nY) + MyzBin*(nY*nX); //YXZ coordinates ordering
  else if(CoorsOrder == 4)  GlobalIndex = MyyBin + MyzBin*(nY) + MyxBin*(nY*nZ); //YZX coordinates ordering
  else if(CoorsOrder == 5)  GlobalIndex = MyzBin + MyxBin*(nZ) + MyyBin*(nZ*nX); //ZXY coordinates ordering
  else if(CoorsOrder == 6)  GlobalIndex = MyzBin + MyyBin*(nZ) + MyxBin*(nZ*nY); //ZYX coordinates ordering
  
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

  //Do nothing if rho and z point are outside fieldmap limits
  if (not ( x >= xMin && x <= xMax &&
            y >= yMin && y <= yMax &&
            z >= zMin && z <= zMax )
     ) {
    return;
  }

  //Calculate the bins on the x, y and z axis containing the (x,y,z) point
  int xBin,yBin,zBin;
  double x0,y0,z0;

  xBin = int((x - xMin)/xStep);
  yBin = int((y - yMin)/yStep);
  zBin = int((z - zMin)/zStep);
  x0   = xMin + xBin*xStep;
  y0   = yMin + yBin*yStep;
  z0   = zMin + zBin*zStep;

  if(x0 > x) {
    x0   -= xStep;
    xBin -= 1;
  }
  if(y0 > y) {
    y0   -= yStep;
    yBin -= 1;
  }
  if(z0 > z) {
    z0   -= zStep;
    zBin -= 1;
  }

  //Get normalized coordinate of (x,y,z) point in bin
  double xd = (x - x0)/xStep;
  double yd = (y - y0)/yStep;
  double zd = (z - z0)/zStep;

  //Get the field values at the eight corners of bin containing the (x,y,z) point
  int xBin0 = xBin;
  int xBin1 = xBin+1;
  int yBin0 = yBin;
  int yBin1 = yBin+1;
  int zBin0 = zBin;
  int zBin1 = zBin+1;
  //Protection in case the sampled coordinate is exactly at maximum value of fieldmap
  if(xBin1 > nX-1) xBin1 = nX-1;
  if(yBin1 > nY-1) yBin1 = nY-1;
  if(zBin1 > nZ-1) zBin1 = nZ-1;
  const FieldMapXYZ::FieldValues_t& B_x0y0z0 = fieldMap[getGlobalIndex(xBin0,yBin0,zBin0)];
  const FieldMapXYZ::FieldValues_t& B_x1y0z0 = fieldMap[getGlobalIndex(xBin1,yBin0,zBin0)];
  const FieldMapXYZ::FieldValues_t& B_x0y0z1 = fieldMap[getGlobalIndex(xBin0,yBin0,zBin1)];
  const FieldMapXYZ::FieldValues_t& B_x1y0z1 = fieldMap[getGlobalIndex(xBin1,yBin0,zBin1)];
  const FieldMapXYZ::FieldValues_t& B_x0y1z0 = fieldMap[getGlobalIndex(xBin0,yBin1,zBin0)];
  const FieldMapXYZ::FieldValues_t& B_x1y1z0 = fieldMap[getGlobalIndex(xBin1,yBin1,zBin0)];
  const FieldMapXYZ::FieldValues_t& B_x0y1z1 = fieldMap[getGlobalIndex(xBin0,yBin1,zBin1)];
  const FieldMapXYZ::FieldValues_t& B_x1y1z1 = fieldMap[getGlobalIndex(xBin1,yBin1,zBin1)];

  //field at (x,y,z) point is linear interpolation of fielmap values at bin corners
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
                                       double coorUnits, double BfieldUnits) {


  TFile *file = TFile::Open( filename.c_str() );
  if (not file) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: File not found: " << filename;
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

  //Set branch adresses
  float x, y, z, Bx, By, Bz;
  checkBranch( tree->SetBranchAddress(xVar.c_str(), &x) );
  checkBranch( tree->SetBranchAddress(yVar.c_str(), &y) );
  checkBranch( tree->SetBranchAddress(zVar.c_str(), &z) );
  checkBranch( tree->SetBranchAddress(BxVar.c_str(), &Bx) );
  checkBranch( tree->SetBranchAddress(ByVar.c_str(), &By) );
  checkBranch( tree->SetBranchAddress(BzVar.c_str(), &Bz) );

  //Loop over the tree entries. In this loop get,
  // - min, max and step-size values of fieldmap coordinates
  // - coordinates ordering
  xStep      = -1;
  yStep      = -1;
  zStep      = -1;
  xOrdering  =  1;
  yOrdering  =  1;
  zOrdering  =  1;
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
  if(xMax < xMin) {
    //x variable is scanned from high-to-low
    xOrdering = -1;
    double aux = xMin;
    xMin = xMax;
    xMax = aux;
  }
  if(yStep < 0) {
    std::stringstream error;
    error << "FieldMapXYZ[ERROR]: All y coordinates in n-tuple have the same value!!!";
    throw std::runtime_error( error.str() );
  }
  if(yMax < yMin) {
    //y variable is scanned from high-to-low
    yOrdering = -1;
    double aux = yMin;
    yMin = yMax;
    yMax = aux;
  }
  if(zStep < 0) {
    std::stringstream error;
    error << "FieldMapBrBz[ERROR]: All z coordinates in n-tuple have the same value!!!";
    throw std::runtime_error( error.str() );
  }
  if(zMax < zMin) {
    //z variable is scanned from high-to-low
    zOrdering = -1;
    double aux = zMin;
    zMin = zMax;
    zMax = aux;
  }

  //Calculate number of bins in fieldmap
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

  //Set coordinates parameters units
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

  std::string  filename   = xmlParameter.attr< std::string >(_Unicode(filename));
  std::string  NtupleName = xmlParameter.attr< std::string >(_Unicode(treeName));
  std::string  xVar       = xmlParameter.attr< std::string >(_Unicode(xVarName));
  std::string  yVar       = xmlParameter.attr< std::string >(_Unicode(yVarName));
  std::string  zVar       = xmlParameter.attr< std::string >(_Unicode(zVarName));
  std::string  BxVar      = xmlParameter.attr< std::string >(_Unicode(BxVarName));
  std::string  ByVar      = xmlParameter.attr< std::string >(_Unicode(ByVarName));
  std::string  BzVar      = xmlParameter.attr< std::string >(_Unicode(BzVarName));

  double xScale = xmlParameter.attr< double >(_Unicode(xScale));
  double yScale = xmlParameter.attr< double >(_Unicode(yScale));
  double zScale = xmlParameter.attr< double >(_Unicode(zScale));
  double bScale = xmlParameter.attr< double >(_Unicode(bScale));

  double coorUnits   = xmlParameter.attr< double >(_Unicode(coorUnits));
  double BfieldUnits = xmlParameter.attr< double >(_Unicode(BfieldUnits));

  DD4hep::Geometry::CartesianField obj;
  FieldMapXYZ* ptr = new FieldMapXYZ();
  ptr->xScale     = xScale;
  ptr->yScale     = yScale;
  ptr->zScale     = zScale;
  ptr->bScale     = bScale;
  ptr->NtupleName = NtupleName;
  ptr->xVar       = xVar;
  ptr->yVar       = yVar;
  ptr->zVar       = zVar;
  ptr->BxVar      = BxVar;
  ptr->ByVar      = ByVar;
  ptr->BzVar      = BzVar;

  //Read the entries form the file in this place
  ptr->fillFieldMapFromTree(filename,coorUnits,BfieldUnits);

  std::string StrXOrdering("low-to-high");
  std::string StrYOrdering("low-to-high");
  std::string StrZOrdering("low-to-high");
  if(ptr->xOrdering == -1) StrYOrdering   = "high-to-low";
  if(ptr->yOrdering == -1) StrYOrdering   = "high-to-low";
  if(ptr->zOrdering == -1) StrZOrdering   = "high-to-low";

  std::cout << "xScale      " << std::setw(13) << ptr->xScale                           << std::endl;
  std::cout << "zScale      " << std::setw(13) << ptr->zScale                           << std::endl;
  std::cout << "zScale      " << std::setw(13) << ptr->zScale                           << std::endl;
  std::cout << "bScale      " << std::setw(13) << ptr->bScale                           << std::endl;
  std::cout << "xMin        " << std::setw(13) << ptr->xMin/dd4hep::cm << " cm"         << std::endl;
  std::cout << "xMax        " << std::setw(13) << ptr->xMax/dd4hep::cm << " cm"         << std::endl;
  std::cout << "xStep       " << std::setw(13) << ptr->xStep/dd4hep::cm << " cm"        << std::endl;
  std::cout << "nX          " << std::setw(13) << ptr->nX                               << std::endl;
  std::cout << "xOrdering   " << std::setw(13) << StrXOrdering.c_str()                  << std::endl;
  std::cout << "yMin        " << std::setw(13) << ptr->yMin/dd4hep::cm << " cm"         << std::endl;
  std::cout << "yMax        " << std::setw(13) << ptr->yMax/dd4hep::cm << " cm"         << std::endl;
  std::cout << "yStep       " << std::setw(13) << ptr->yStep/dd4hep::cm << " cm"        << std::endl;
  std::cout << "nY          " << std::setw(13) << ptr->nY                               << std::endl;
  std::cout << "yOrdering   " << std::setw(13) << StrYOrdering.c_str()                  << std::endl;
  std::cout << "zMin        " << std::setw(13) << ptr->zMin/dd4hep::cm << " cm"         << std::endl;
  std::cout << "zMax        " << std::setw(13) << ptr->zMax/dd4hep::cm << " cm"         << std::endl;
  std::cout << "zStep       " << std::setw(13) << ptr->zStep/dd4hep::cm << " cm"        << std::endl;
  std::cout << "nZ          " << std::setw(13) << ptr->nZ                               << std::endl;
  std::cout << "zOrdering   " << std::setw(13) << StrZOrdering.c_str()                  << std::endl;
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

