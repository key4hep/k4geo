#!/bin/bash

ILCSOFT=/cvmfs/clicdp.cern.ch/iLCSoft/builds/current/CI_${COMPILER}
source $ILCSOFT/init_ilcsoft.sh

cd /Package
mkdir build
cd build
cmake -C$ILCSOFT/ILCSoft.cmake -GNinja -DUSE_CXX11=ON -DCMAKE_CXX_FLAGS="-fdiagnostics-color=always -Werror" .. && \
    ninja -k 0 && \
    ninja install && \
    ctest --output-on-failure
