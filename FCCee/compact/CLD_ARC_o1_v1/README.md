ARC detector V0
=======================

Array of RICH Cells (ARC) is a novel RICH detector concept. The ARC detector is made up by two parts, a barrel and two endcaps. The barrel consist in the repetition of 16 unique cells, 27 times around phi, and mirrored over z axis. Each endcap consist in the repetition of 6 times of a sector. One sector is made of 21 unique cells, some of the reflected within the sector. Each unique cell
* is filled by C4F10 radiator gas
* has a spherical mirror in the outer radial part of the cell
* and a sandwich of three layers in the inner radial part, from inner to outer, a cooling plate, a light sensor, and an aerogel radiator.

Each cell has 5 characteristic parameters, which correspond to
* the position of the mirror (2),
* the position and angle of the sensor (3)
* optional parameter (`xxx_DetZOffset`) can be given to increase the distance of the sensor from the vessel, so the volumes do not overlap
There parameters are defined as `DD4hep` constants in the `RadiatorCell_FinalOptimization.xml` file.

The original design included few more cells, 2 more on the edge of the barrel, and some more in the endcap. These cells are not included in the current design because they require some extra work to be properly implemented.

It has a radial depth of 20 cm, outer radius of 2.1 m and a length of 4.4 m. Such dimensions are hardcoded now. In a later version, global dimensions and some other internal parameters will be setup freely in the compact file.

Material description is taken from [Proximity Focusing RICH (pfRICH) detector example in DD4hep](https://github.com/AIDASoft/DD4hep/tree/master/examples/OpticalTracker)

Only the optical surface for the mirror is taken from the material compact file. In a later version, the optical surface for the sensor will be used to take into account the possible light detection efficiency. In addition, some optical surface for the vessel is needed to properly recreate a possible light background.

Some presentation about the ARC design can be found here:
* Dec-22 ARC meeting indico: https://indico.cern.ch/event/1231098/
* Mar-28 Discussion about ARC implementation in DD4hep: https://indico.cern.ch/event/1266428/
