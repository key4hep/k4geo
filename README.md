# lcgeo (Linear Collider Geometry)
[![Build Status](https://travis-ci.org/iLCSoft/lcgeo.svg?branch=master)](https://travis-ci.org/iLCSoft/lcgeo)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/12359/badge.svg)](https://scan.coverity.com/projects/ilcsoft-lcgeo)

Implementation of Linear Collider detector models in DD4hep.

lcgeo is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Requirements
DD4hep built with Geant4 and LCIO
## Download and Installation
### Download
  * `git clone https://github.com/iLCSoft/lcgeo.git`

  * `cd lcgeo ; mkdir build ; cd build`
### Initialize dependency
  
  * `source __path_to_DD4hep__/bin/thisdd4hep.sh`
  
  * `cmake -DBoost_NO_BOOST_CMAKE=ON ..`
### Build
  * `make -j4 install`
### Initialize
  * `source ../bin/thislcgeo.sh`

## Running some examples
go to `cd ./example`

### Create single particle input file
Modify the python script `./example/lcio_particle_gun.py` in order to create
an LCIO input file with single particles:
  * `export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib:$PYTHONPATH`
  * `python lcio_particle_gun.py`

There is also example input file with 10 singe muons: ./mcparticles.slcio

### Run the simulation:
   * `ddsim --compactFile ../ILD/compact/ILD_o1_v05/ILD_o1_v05.xml --inputFiles mcparticles.slcio -N 10`

This creates an lcio file simple_lcio.slcio with sim hits and MCParticles.

You can look at it in the usual way

   * `anajob simple_lcio.slcio`

   * `dumpevent simple_lcio.slcio 1`

Change the ddsim command line parameters as needed to read other input files.

## Event displays:

There are several ways for visualizing the detector geometry and the simulated events:

1) CED event display

   `ced2go -d gear_ILD_o1_v05_ORG.xml -t ced2go-template.xml simple_lcio.slcio`


2) teveDisplay

   `ln -s simple_lcio.slcio teve_infile.slcio`
   
   `teveDisplay ../ILD/compact/ILD_o1_v05/ILD_o1_v05.xml`


3) DDEve:

   `root`
   
   `.x $DD4hepINSTALL/examples/DDEve/DDEve.C()`
   
   and then load
  
   `../ILD/compact/ILD_o1_v05/DDEve.xml`

## Running the reconstruction with Marlin (EXPERIMENTAL)
 
 1) create a gear file from the DD4hep detector model, e.g.

    `convertToGear default ../trunk/ILD/compact/ILD_o1_v05.xml gear_ILD_o1_v05_dd4hep.xml`

    [ Friendly reminder: in order to work you need to explicitly build DD4hep with GEAR enabled. Use the flag `-DDD4HEP_USE_GEAR="ON"` for the cmake command. ]

    Note: currently this might still be incomplete and you can copy missing information for example from the original Mokka gear files ( or use `$DDSIM/example/gear_ILD_o1_v05_dd4hep.xml` )


 2) run simulation (as described above), e.g.:

    `python ddsim.py ../ILD/compact/ILD_o1_v05/ILD_o1_v05.xml`



 3) run the standard reconstruction ( tracking only so far ):

    `Marlin ../example/ild_dd4hep_stdreco.xml`  
    `--global.LCIOInputFiles=./simple_lcio.slcio`
    `--MyLCIOOutputProcessor.LCIOOutputFile=./simple_lcio_REC.slcio`

 4) look at the result w/ the CED event display, e.g.
 
    `ced2go -d ../example/gear_ILD_o1_v05_ORG.xml -t ../example/ced2go-template.xml simple_lcio_REC.slcio`

## License and Copyright
Copyright (C), lcgeo Authors

lcgeo is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License long with this program.  If not, see <http://www.gnu.org/licenses/>.
