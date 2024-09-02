aciarma - 08/07/24

- MDI_o1_v00
-- Beampipe_o4_v05.xml : "shape based" beam pipe with AlBeMet central chamber and Paraffin cooling. From IDEA_o1_v02.
-- BeamInstrumentation_o1_v01.xml : compensating and screening solenoids
-- HOMAbsorber.xml : old model of the HOM absorbers. Not needed anymore with new low impedance beam pipe.
-- MDI_Dimensions_Materials_o1_v00.xml : Includes the materials and constants necessary for the MDI elements
-- LumiCal_o3_v02_05.xml : lumical xml
-- FFQuads_v01.xml : simple equivalent material geometry for the final focus quadrupoles
-- FFQuads_params_xxx_x.xml : parameters for geometry and fields of the final focus quadrupoles, different lattices and energies
-- fields_antisol_map.xml : imports the magnetic fields for the compensating and screening solenoids from field map (NO 2T DETECTOR FIELD)
-- fields_antisol_ideal.xml : ideal fields for compensating and screening solenoids, no fringe fields (NO 2T DETECTOR FIELD)
only for debug
-- MDI_standalone_o1_v00.xml : master compact for only the MDI elements
-- fields_solenoid.xml : simple 2T detector solenoid field
-- local_elements.xml
-- local_materials.xml


- MDI_o1_v01
-- Beampipe_CADimport_o1_v02.xml : import CAD models for engineered beam pipe (by F. Fransesini/INFN-LNF)
These .stl files are hosted [here](https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/MDI/MDI_o1_v01/).
The CMake option `INSTALL_BEAMPIPE_STL_FILES=ON` downloads these STL.
-- stl_files/Pipe_240430
    ├── AlBeMet162_30042024.stl    : central and elliptoconical chambers, with cooling manifolds 
    ├── Copper_pipe_28092023.stl   : low impedance beam pipe separation region
    ├── Gold_19042024.stl          : 5um coating inside the central chamber 
    ├── Paraffine_19042024.stl     : cooling for central chamber
    ├── Tungsten_mask_02102023.stl : SR masks 2.1m upstream
    └── Water_30042024.stl         : cooling for elliptoconical chambers
-- BeamInstrumentation_o1_v01.xml : compensating and screening solenoids
-- MDI_Dimensions_Materials_o1_v01.xml : Includes the materials and constants necessary for the MDI elements
-- LumiCal_o3_v02_05.xml : lumical xml
-- FFQuads_v01.xml : simple equivalent material geometry for the final focus quadrupoles
-- FFQuads_params_xxx_x.xml : parameters for geometry and fields of the final focus quadrupoles, different lattices and energies
-- fields_antisol_map.xml : imports the magnetic fields for the compensating and screening solenoids from field map (NO 2T DETECTOR FIELD)
-- fields_antisol_ideal.xml : ideal fields for compensating and screening solenoids, no fringe fields (NO 2T DETECTOR FIELD)
only for debug
-- MDI_standalone_o1_v01.xml : master compact for only the MDI elements
-- fields_solenoid.xml : simple 2T detector solenoid field
-- local_elements.xml
-- local_materials.xml



