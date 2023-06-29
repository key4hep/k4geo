ARC detector o1 v01
=======================

The Array of RICH Cells (ARC) is a novel RICH detector concept. The ARC detector is made up of three `Detector Elements`, a barrel and two endcaps. Both C++ constructors are declared in the same file `ARC_geo_o1_v01.cpp`. 

The barrel is constructed through the repetition of 17 unique cells, 27 times around phi, and mirrored over the XY plane. Each endcap consist of six repetition of a sector. One sector is made out of 21 unique cells, some of them are reflected within the sector. Further details and figures can be found in [Martin Tat's talk](https://indico.cern.ch/event/1231098/contributions/5179993/attachments/2568014/4427756/ARC_Presentation_DD4HEPIntro_15th_December_2022.pdf). The original design included an additional half cell at the edge of the barrel (cell number 18), and some more cells in the endcap. These cells are not included in the current design because they require some extra work to be properly implemented without causing overlaps.


The volume of each cell

* is filled by a radiator gas (e.g. C4F10), as specified in the compact file;

* has a spherical mirror in the outer radial part of the cell, clipped to fit inside the cell volume;

* has a sandwich of three layers in the inner radial part, from inner to outer, a cooling plate, a light sensor, and an aerogel radiator;

* has vessel walls are 1 cm thick.

The current detector has a radial depth of 20 cm, an outer radius of 2.1 m and a length of 4.4 m. These dimensions are hard-coded now. In a later version, global geometry parameters and some other internal parameters will be set up freely in the compact file.

The uniqueness of the cells happens because the position of the sensor and mirror is optimized for each cell, and defined by 5 parameters. They correspond to
* the position of the mirror (2 parameters),
* the position and angle of the sensor (3 parameters)
* an optional parameter (`xxx_DetPositionZ`) can be given to increase the distance of the sensor from the vessel, so that the volumes do not overlap.

The parameters are defined as `DD4hep` constants in the `RadiatorCell_FinalOptimization.xml` file. The value of these parameters is linked to the geometry description. If the geometry of the ARC detector (namely its radius and thickness) is changed, these parameters should be reoptimized by a dedicated ray-tracing software.

The material description is taken from the [Proximity Focusing RICH (pfRICH) detector example in DD4hep](https://github.com/AIDASoft/DD4hep/tree/master/examples/OpticalTracker). The optical surface for the mirror, as well as the materials, are defined in the compact file `materials_arc_o1_v01.xml`. In a later version, the optical surface for the sensor will be used to take into account the possible light detection efficiency. In addition, some optical surfaces for the vessel are needed to properly recreate a possible light background (not included at the moment). The Aerogel material without optical properties is used as template for the bulk part of the vessel walls.

Some presentation about the ARC design can be found here:
* Dec-22 ARC meeting indico: https://indico.cern.ch/event/1231098/
* Mar-28 Discussion about ARC implementation in DD4hep: https://indico.cern.ch/event/1266428/
* Analysis of the data can be adapted from [previous RICH detectors](https://s3.cern.ch/inspire-prod-files-9/92927eb16166b155de56b61339f05521)
