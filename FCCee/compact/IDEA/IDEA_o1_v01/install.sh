### Source script to install k4geo 

# using the nightly release
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh                                                                                                                                                                                                                                                                                                              

cd ../../../
mkdir build
cd build
source $DD4hep_DIR/bin/thisdd4hep.sh

cmake .. -DCMAKE_INSTALL_PREFIX=../InstallArea -DBoost_NO_BOOST_CMAKE=ON
make -j4 install
cd ../InstallArea
source bin/thislcgeo.sh

export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib:$PYTHONPATH
