aciarma - 08/07/24

- MDI_o1_v00
-- Beampipe_o4_v05.xml : "shape based" beam pipe with AlBeMet central chamber and Paraffin cooling. From IDEA_o1_v02.
-- BeamInstrumentation_o1_v01.xml : compensating and screening solenoids
-- HOMAbsorber.xml : old model of the HOM absorbers. Not needed anymore with new low impedance beam pipe.

- MDI_o1_v01
-- Beampipe_CADimport_o1_v02.xml : import CAD models for engineered beam pipe (by F. Fransesini/INFN-LNF)
These .stl files are hosted [here](https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/MDI/MDI_o1_v01/).
The CMake option `INSTALL_BEAMPIPE_STL_FILES=ON` downloads these STL.
```
-- stl_files/Pipe_240430
    ├── AlBeMet162_30042024.stl    : central and elliptoconical chambers, with cooling manifolds 
    ├── Copper_pipe_28092023.stl   : low impedance beam pipe separation region
    ├── Gold_19042024.stl          : 5um coating inside the central chamber 
    ├── Paraffine_19042024.stl     : cooling for central chamber
    ├── Tungsten_mask_02102023.stl : SR masks 2.1m upstream
    └── Water_30042024.stl         : cooling for elliptoconical chambers
-- BeamInstrumentation_o1_v01.xml : compensating and screening solenoids
```

- MDI_o1_CADBased_v01
-- continuation of MDI_o1_v01, fixed issue of a mis-rotation of the beam pipe separation
Added MDI_standalone.xml to help debug MDI geometry