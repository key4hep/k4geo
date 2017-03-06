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
    globalField[0] = globalField[1] = globalField[2] = 0.0;
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

  const FieldMapXYZ::FieldValues_t& B_x0y0z0 = fieldMap[GetGlobalIndex(x0,y0,z0)];
  const FieldMapXYZ::FieldValues_t& B_x1y0z0 = fieldMap[GetGlobalIndex(x1,y0,z0)];
  const FieldMapXYZ::FieldValues_t& B_x0y0z1 = fieldMap[GetGlobalIndex(x0,y0,z1)];
  const FieldMapXYZ::FieldValues_t& B_x1y0z1 = fieldMap[GetGlobalIndex(x1,y0,z1)];
  const FieldMapXYZ::FieldValues_t& B_x0y1z0 = fieldMap[GetGlobalIndex(x0,y1,z0)];
  const FieldMapXYZ::FieldValues_t& B_x1y1z0 = fieldMap[GetGlobalIndex(x1,y1,z0)];
  const FieldMapXYZ::FieldValues_t& B_x0y1z1 = fieldMap[GetGlobalIndex(x0,y1,z1)];
  const FieldMapXYZ::FieldValues_t& B_x1y1z1 = fieldMap[GetGlobalIndex(x1,y1,z1)];

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
  
  return; 
 
}

int FieldMapXYZ::GetGlobalIndex(double x, double y, double z)
{

  int xBin = int((x - xMin)/xStep);
  int yBin = int((y - yMin)/yStep);
  int zBin = int((z - zMin)/zStep);

  return  xBin + yBin*nX + zBin*nX*nY;

}

void FieldMapXYZ::fillFieldMapFromTree(const std::string& filename, const std::string& treeString){

  TFile *file = TFile::Open( filename.c_str() );
  if (not file) {
    std::stringstream error;
    error << "FieldMap[ERROR]: File not found: " << filename;
    throw std::runtime_error( error.str() );
  }

  std::vector<std::string> strs;
  boost::split(strs, treeString, boost::is_any_of(": "));
  for (std::vector<std::string>::iterator it = strs.begin(); it != strs.end(); ++it) {
    std::cout << "\"" << (*it).c_str() << "\"" << std::endl;
  }

  if( strs.size() != 5 ) {
    std::stringstream error;
    error << "FieldMap[ERROR]: the treeDescription " << treeString
	  << " is not complete. For example 'fieldmap:rho:z:Brho:Bz'";
    throw std::runtime_error( error.str() );
  }

  TTree *tree;
  file->GetObject(strs[0].c_str(), tree);
  if (not tree) {
    std::stringstream error;
    error << "FieldMap[ERROR]: Tree " << strs[0] << " not found in file: " << filename;
    throw std::runtime_error( error.str() );
  }

  float x, y, z, Bx, By, Bz;
  checkBranch( tree->SetBranchAddress(strs[1].c_str(), &x) );
  checkBranch( tree->SetBranchAddress(strs[2].c_str(), &y) );
  checkBranch( tree->SetBranchAddress(strs[3].c_str(), &z) );
  checkBranch( tree->SetBranchAddress(strs[4].c_str(), &Bx) );
  checkBranch( tree->SetBranchAddress(strs[5].c_str(), &By) );
  checkBranch( tree->SetBranchAddress(strs[6].c_str(), &Bz) );

  xMin -= 0.5*xStep;
  xMax += 0.5*xStep;
  yMin -= 0.5*yStep;
  yMax += 0.5*yStep;
  zMin -= 0.5*zStep;
  zMax += 0.5*zStep;
  
  nX = ((xMax - xMin)/xStep);
  nY = ((yMax - yMin)/yStep);
  nZ = ((zMax - zMin)/zStep);

  const int elements = nX*nY*nZ;
  const int treeEntries = tree->GetEntries();
  if ( elements != treeEntries ) {
    std::stringstream error;
    error << "FieldMap[ERROR]: Tree does not have the same size as described in the XML "
	  << "nX*nY*nZ == "  << elements
	  << "  tree entries == " << treeEntries;
    throw std::runtime_error( error.str() );
  }

  fieldMap.reserve(elements);

  for(int i = 0;i<treeEntries;++i) {
    tree->GetEntry(i);
    fieldMap.push_back( FieldMapXYZ::FieldValues_t(double(Bx)*bScale*dd4hep::tesla,
                                                   double(By)*bScale*dd4hep::tesla,
						   double(Bz)*bScale*dd4hep::tesla ) );

  }

  file->Close();
  delete file;

}

static DD4hep::Geometry::Ref_t create_FieldMap_XYZ(DD4hep::Geometry::LCDD& ,
						   DD4hep::XML::Handle_t handle ) {
  DD4hep::XML::Component xmlParameter(handle);
  bool hasFilename = xmlParameter.hasAttr(_Unicode(filename));

  if (!hasFilename) {
    std::stringstream error;
    error << "FieldMap[ERROR]: For a FieldMap field at least the filename xml attribute MUST be set.";
    throw std::runtime_error(error.str());
  }
  std::string filename   = xmlParameter.attr< std::string >(_Unicode(filename));
  std::string treeString = xmlParameter.attr< std::string >(_Unicode(tree));

  double xScale = xmlParameter.attr< double >(_Unicode(xScale));
  double yScale = xmlParameter.attr< double >(_Unicode(yScale));
  double zScale = xmlParameter.attr< double >(_Unicode(zScale));
  double bScale = xmlParameter.attr< double >(_Unicode(bScale));

  double xMin   = xmlParameter.attr< double >(_Unicode(xMin));
  double xMax   = xmlParameter.attr< double >(_Unicode(xMax));
  double xStep  = xmlParameter.attr< double >(_Unicode(xStep));

  double yMin   = xmlParameter.attr< double >(_Unicode(yMin));
  double yMax   = xmlParameter.attr< double >(_Unicode(yMax));
  double yStep  = xmlParameter.attr< double >(_Unicode(yStep));

  double zMin   = xmlParameter.attr< double >(_Unicode(zMin));
  double zMax   = xmlParameter.attr< double >(_Unicode(zMax));
  double zStep  = xmlParameter.attr< double >(_Unicode(zStep));

  DD4hep::Geometry::CartesianField obj;
  FieldMapXYZ* ptr = new FieldMapXYZ();
  ptr->xScale = xScale;
  ptr->yScale = yScale;
  ptr->zScale = zScale;
  ptr->bScale = bScale;

  ptr->xMin   = xMin*xScale;
  ptr->xMax   = xMax*xScale;
  ptr->xStep  = xStep*xScale;

  ptr->yMin   = yMin*yScale;
  ptr->yMax   = yMax*yScale;
  ptr->yStep  = yStep*yScale;

  ptr->zMin   = zMin*zScale;
  ptr->zMax   = zMax*zScale;
  ptr->zStep  = zStep*zScale;

  std::cout << "xScale  " << std::setw(13) << xScale  << std::endl;
  std::cout << "yScale  " << std::setw(13) << yScale  << std::endl;
  std::cout << "zScale  " << std::setw(13) << zScale  << std::endl;
  std::cout << "bScale  " << std::setw(13) << bScale  << std::endl;

  std::cout << "xMin    " << std::setw(13) << xMin    << std::endl;
  std::cout << "xMax    " << std::setw(13) << xMax    << std::endl;
  std::cout << "xStep   " << std::setw(13) << xStep   << std::endl;

  std::cout << "yMin    " << std::setw(13) << yMin    << std::endl;
  std::cout << "yMax    " << std::setw(13) << yMax    << std::endl;
  std::cout << "yStep   " << std::setw(13) << yStep   << std::endl;

  std::cout << "zMin    " << std::setw(13) << zMin    << std::endl;
  std::cout << "zMax    " << std::setw(13) << zMax    << std::endl;
  std::cout << "zStep   " << std::setw(13) << zStep   << std::endl;

  //Read the entries form the file in this place
  ptr->fillFieldMapFromTree(filename, treeString);

  obj.assign(ptr, xmlParameter.nameStr(), xmlParameter.typeStr());
  return obj;
}
DECLARE_XMLELEMENT(FieldXYZ,create_FieldMap_XYZ)

