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
  const double posRho = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
  const double posZ   = pos[2];
  const double posPhi = atan2(pos[1], pos[0]);

  //Mokka defines these things as x and y and phi
  double x = posRho;
  double y = posZ;
  double phi = posPhi;

  //Get positive values to do less checks when comparing
  if( y < 0 ) {
    y = -y ;
    phi += M_PI ;
  }

  //APS: Note the mokka field map does not start at 0, so we have to assume that
  //this area is covered, or add some more parameters for the values where the
  //field is supposed to be. Now we just assume it starts at 0/0/0, outside there is no field
  if (not ( 0.0 <= x && x <= rhoMax &&
	    0.0 <= y && y <= zMax ) ) {

    return;
  }

  if( x < rhoMin ) x = rhoMin;
  if( y < zMin ) y = zMin ;


  // Position of given point within region, normalized to the range
  // [0,1]
  const double xfraction = ( x - rhoMin ) / ( rhoMax - rhoMin );
  const double yfraction = ( y - zMin )   / ( zMax   - zMin   );

  // Need addresses of these to pass to modf below.
  // modf uses its second argument as an OUTPUT argument.
  double xdindex, ydindex;

  // Position of the point within the cuboid defined by the
  // nearest surrounding tabulated points
  const double xlocal = ( std::modf( xfraction * ( nRho - 1), &xdindex ) );
  const double ylocal = ( std::modf( yfraction * ( nZ   - 1), &ydindex ) );

  // The indices of the nearest tabulated point whose coordinates
  // are all less than those of the given point
  const int xindex = static_cast<int>(xdindex);
  const int yindex = static_cast<int>(ydindex);

  const int index0 =  xindex    * nZ + yindex   ;
  const int index1 =  xindex    * nZ + yindex+1 ;
  const int index2 = (xindex+1) * nZ + yindex   ;
  const int index3 = (xindex+1) * nZ + yindex+1 ;

  const FieldMapXYZ::FieldValues_t& fv0 = fieldMap[index0];
  const FieldMapXYZ::FieldValues_t& fv1 = fieldMap[index1];
  const FieldMapXYZ::FieldValues_t& fv2 = fieldMap[index2];
  const FieldMapXYZ::FieldValues_t& fv3 = fieldMap[index3];

  double field[2] = {0.0, 0.0};

  field[0] =
    fv0.Br * (1-xlocal) * (1-ylocal)  +
    fv1.Br * (1-xlocal) *    ylocal   +
    fv2.Br *    xlocal  * (1-ylocal)  +
    fv3.Br *    xlocal  *    ylocal   ;

  field[1] =
    fv0.Bz * (1-xlocal) * (1-ylocal)  +
    fv1.Bz * (1-xlocal) *    ylocal   +
    fv2.Bz *    xlocal  * (1-ylocal)  +
    fv3.Bz *    xlocal  *    ylocal   ;

  globalField[0] += field[0] * sin( phi ) ;
  globalField[1] += field[0] * cos( phi ) ;
  globalField[2] += field[1] ;


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

  float r, z, Br, Bz;
  checkBranch( tree->SetBranchAddress(strs[1].c_str(), &r) );
  checkBranch( tree->SetBranchAddress(strs[2].c_str(), &z) );
  checkBranch( tree->SetBranchAddress(strs[3].c_str(), &Br) );
  checkBranch( tree->SetBranchAddress(strs[4].c_str(), &Bz) );

  const int elements = nRho*nZ;
  const int treeEntries = tree->GetEntries();
  if ( elements != treeEntries ) {
    std::stringstream error;
    error << "FieldMap[ERROR]: Tree does not have the same size as described in the XML "
	  << "nRho*nZ == "  << elements
	  << "  tree entries == " << treeEntries;
    throw std::runtime_error( error.str() );
  }

  fieldMap.reserve(elements);

  for (int i = 0; i < treeEntries ;++i) {
    tree->GetEntry(i);
    fieldMap.push_back( FieldMapXYZ::FieldValues_t( double(Br)*bScale*dd4hep::tesla,
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
  std::string filename = xmlParameter.attr< std::string >(_Unicode(filename));
  std::string treeString = xmlParameter.attr< std::string >(_Unicode(tree));

  double rScale = xmlParameter.attr< double >(_Unicode(rScale));
  double zScale = xmlParameter.attr< double >(_Unicode(zScale));
  double bScale = xmlParameter.attr< double >(_Unicode(bScale));
  double rhoMin = xmlParameter.attr< double >(_Unicode(rhoMin));
  double rhoMax = xmlParameter.attr< double >(_Unicode(rhoMax));
  double nRho   = xmlParameter.attr< int >   (_Unicode(nRho));
  double zMin   = xmlParameter.attr< double >(_Unicode(zMin));
  double zMax   = xmlParameter.attr< double >(_Unicode(zMax));
  double nZ     = xmlParameter.attr< int >   (_Unicode(nZ));

  DD4hep::Geometry::CartesianField obj;
  FieldMapXYZ* ptr = new FieldMapXYZ();
  ptr->rScale = rScale;
  ptr->zScale = zScale;
  ptr->bScale = bScale;
  ptr->rhoMin = rhoMin*rScale;
  ptr->rhoMax = rhoMax*rScale;
  ptr->nRho   = nRho  ;
  ptr->zMin   = zMin  *zScale;
  ptr->zMax   = zMax  *zScale;
  ptr->nZ     = nZ    ;

  std::cout << "rScale  " << std::setw(13) << rScale  << std::endl;
  std::cout << "zScale  " << std::setw(13) << zScale  << std::endl;
  std::cout << "bScale  " << std::setw(13) << bScale  << std::endl;
  std::cout << "rhoMin  " << std::setw(13) << rhoMin  << std::endl;
  std::cout << "rhoMax  " << std::setw(13) << rhoMax  << std::endl;
  std::cout << "nRho    " << std::setw(13) << nRho    << std::endl;
  std::cout << "zMin    " << std::setw(13) << zMin    << std::endl;
  std::cout << "zMax    " << std::setw(13) << zMax    << std::endl;
  std::cout << "nZ      " << std::setw(13) << nZ      << std::endl;

  //Read the entries form the file in this place
  ptr->fillFieldMapFromTree(filename, treeString);

  obj.assign(ptr, xmlParameter.nameStr(), xmlParameter.typeStr());
  return obj;
}
DECLARE_XMLELEMENT(FieldXYZ,create_FieldMap_XYZ)

