# lcgeo (Linear Collider Geometry)
![Build status](https://gitlab.cern.ch/CLICdp/lcgeo/badges/master/build.svg)

Implentaion of Linear Collider detector models in DD4hep.

## Requirements
DD4hep built with Geant4 and LCIO
## Download and Installation
### Downlaod
  * `svn co https://svnsrv.desy.de/public/ddsim/lcgeo/trunk lcgeo`
  
  * (or if you have developers rights:   `svn co https://svnsrv.desy.de/basic/ddsim/lcgeo/trunk lcgeo`)
  
  * `cd trunk ; mkdir build ; cd build`
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
   * `python ddsim.py ../ILD/compact/ILD_o1_v05/ILD_o1_v05.xml`

This creates an lcio file simple_lcio.slcio with sim hits and MCParticles.

You can look at it in the ususal way

   * `anajob simple_lcio.slcio`

   * `dumpevent simple_lcio.slcio 1`

Modify ddsim.py as needed to read other input files.

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
  
  

