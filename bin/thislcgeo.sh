#!/bin/bash
#################################################################################
#
#  Environment script for DD4hep examples - initializes DD4hep (and ROOT)
#  for package: lcgeo
# 
#  @author F.Gaede, DESY, 2013
#  @author M.Frank, CERN, 2015
#
#################################################################################
# Default of DD4hep is the primary installation directory
if [ ! ${DD4hep_DIR} ]; then
    export DD4hep_DIR=/opt/ilcsoft/v02-02-03-pixel/DD4hep/v01-11-02;
fi;
if [ /opt/ilcsoft/CLHEP/2.3.4.3/lib/CLHEP-2.3.4.3 ]; then
    export CLHEP_DIR=/opt/ilcsoft/CLHEP/2.3.4.3/lib/CLHEP-2.3.4.3;
fi;
source ${DD4hep_DIR}/bin/thisdd4hep.sh;
#
dd4hep_parse_this ${BASH_ARGV[0]} lcgeo;
#
#----PATH---------------------------------------------------------------------
dd4hep_add_path    PATH ${THIS}/bin;
#----PYTHONPATH---------------------------------------------------------------
dd4hep_add_path    PYTHONPATH ${THIS}/;
#----ROOT_INCLUDE_PATH--------------------------------------------------------
dd4hep_add_path    ROOT_INCLUDE_PATH ${THIS}/include;
#----LIBRARY_PATH-------------------------------------------------------------
dd4hep_add_library_path ${THIS}/lib;
# -- need to extend dynamic search path for all external libraries:
if [  ]; then
    for lp in ; do
	dd4hep_add_library_path ${lp};
    done;
fi;
