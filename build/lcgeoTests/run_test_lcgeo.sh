#!/bin/bash
#
# Simple script to run DD4hep tests
# - sources this${PackageName}.sh first and then
#   calls the command (given as first argument)
#   with all following arguments
#

#----- initialize environment for this package - including DD4hep 
export DD4hepExamplesINSTALL=/opt/ilcsoft/v02-02-03-pixel/lcgeo/v00-16-07;
source /opt/ilcsoft/v02-02-03-pixel/lcgeo/v00-16-07/bin/thislcgeo.sh;
#----- parse command line - first argument is the 
#      test to run
command=$1
theargs=""
shift
for i in "$@" ; do
    if [ $i != ${command} ] ; then 
	theargs="${theargs} $i"
    fi
done

echo " #### LD_LIBRARY_PATH = :  ${LD_LIBRARY_PATH}"

echo " ### running test :  '${command} ${theargs}'"
exec ${command} ${theargs}
