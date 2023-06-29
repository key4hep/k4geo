CLD_o3_v01 Detector
======================

This option (o3) of CLD includes the Array of RICH Cells (ARC) subdetector. The ARC is placed after the tracker and before the calorimeters. 

# PID assisted by ARC

Array of RICH Cells (ARC) is a novel RICH detector concept. Detailed description of it can be found in [this dedicated README file](../../../detector/PID/README.md).

The identification of the particle crossing the ARC detector is based on the Cherenkov light produced by two components: a radiator gas, and a layer of aerogel.

Similar systems have been used in the past, e.g. at the LHCb experiment. Details of the setup and analysis of the data can be found [here](https://s3.cern.ch/inspire-prod-files-9/92927eb16166b155de56b61339f05521).



## ARC detector compact files

There are 4 compact `.xml` files that are needed for building the detector.

1. `RadiatorCell_FinalOptimisation_o1_v01.xml`. Contains the parameters that determine the position of mirrors and sensors inside the cells. See the *REMINDER* below.

2. `materials_arc_o1_v01.xml`. Contains the description of the different materials needed in the ARC, including their optical properties if needed and the optical surfaces.

3. `elements.xml`. Contains the whole periodic table of the elements. This file is a copy of the one provided by DD4hep.

4. `ARC_o1_v01.xml`. Contains the detector section which builds the detector elements `ARCENDCAP_o1_v01_T` and `ARCBARREL_o1_v01_T`, visualization attributes, the segmentation and readout of the sensors, and specifications to build the detector:

* Material of each component: vessel, radiator gas, light sensor, aerogel layer, cooling plate

* Total size of the light sensor, now implemented as rectangular (X length can be different from Y)

* Thickness of mirror and sensor

* Z position of endcap: The variable `ARC_ENDCAP_ZPOS` defines the middle point of the endcap along the Z axis. Therefore, the endcap spans over `ARC_ENDCAP_ZPOS +- ARC_ENDCAP_LENGTH/2` along the Z axis.

* Readout and segmentation of the light sensors. The readout includes the definition of the Volume ID bit-field. Barrel field corresponds to barrel (0) or endcaps (+/-1), cellnumber fields correspond to a number which is given consecutively as the cells are placed.

*REMINDER*: the cell parameters defined in the file `RadiatorCell_FinalOptimization.xml` were optimized by a dedicated ray-tracing dedicated software for the initial geometry of the ARC detector (radial depth of 20 cm, outer radius of 2.1 m and a length of 4.4 m). If geometry of the ARC changes, these cell parameters should be optimized again.
