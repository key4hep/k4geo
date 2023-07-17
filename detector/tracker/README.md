# Vertex detectors

## VertexBarrel_detailed and VertexDisks_detailed description
These two detector constructors were derived from ZPlanarTracker.cpp and from VertexEndcap_o1_v06.cpp and adapted to fit the needs of the [IDEA vertex detector engineering design](https://indico.cern.ch/event/1202105/timetable/#242-mechanical-integration-of).

Both the barrel and the disks are made up of staves, that can feature an arbitrary number of cuboid layers to describe the support, readout and sensor structures. The sensors can be described by individual sensitive and insensitive pieces to form large, complex structures such as quad modules, that have an insensitive periphery.
More details can be found in the talks at the [FCC week 2023](https://indico.cern.ch/event/1202105/timetable/#356-idea-vertex-detector-in-ke) with an update the the [MDI meeting of the 10th of July 2023](https://indico.cern.ch/event/1292318/#5-vxd-implementation-in-full-s). Once public, the constructor and the resulting IDEA vertex detector are described in the MDI note for the mid-term review.

# Trackers

## parametrised_SimplifiedDriftChamber
This sub-detector implements a simplied drift chamber (there is no actual wire). It is used in ALLEGRO and IDEA detector concepts.

### o1_v01 
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/main/Detector/DetSensitive/src/SimpleDriftChamber.cpp). 

## DriftChamber
Detailed (i.e. including wires) description of a unique-volume, high granularity, fully stereo, low mass, cylindrical drift chamber. 
Originally designed for the IDEA detector, it will be used also for the ALLEGRO concept with different dimensions. 
A detailed description of the detector it tries to model can be found [here](https://indico.cern.ch/event/932973/contributions/4041314/attachments/2139657/3664808/primavera_FCCworkshop_2020.pdf) or [here](https://indico.cern.ch/event/1283129/contributions/5476695/attachments/2682660/4654170/DeFilippis_DCH_CC.pdf).

### o1_v01
Very first version of the detailed drift chamber with wires, walls, etc. Taken from an early prototype available in FCCDetectors ([link](https://github.com/HEP-FCC/FCCDetectors/blob/main/Detector/DetFCCeeIDEA/src/DriftChamber.cpp)) and started to debug it, clean it and add the sensitive volume definition.
It is clearly not finished but made available to enable further technical developments on the digitization, tracking, particle flow, ...
This version does not rely on any external segmentation, all volumes are defined in the detector builder.
The sensitive cell volume is a rotated (to get the stereo angle and phi position) tube segment that is put inside a hyperboloid for each layer.
The hyperboloid is needed for performance, it helps having a well balanced volume tree.
What to improve:
- sensitive volume definition: fix volume extrusion/overlaps/wholes which can be partially solved by using "intersection solids" but this in turn makes the geometry building very slow and memory demanding (~10GB RAM). Other shapes like twisted tubes (caveat: no TGeo shape equivalent), tessellated solids or extruded volumes should be investigated for the next iteration.
- better automation in the C++ to rely less on user defined parameters from the xml
- better variable naming (started already but could still be improved)
- make sense out of the many layer radiuses defined

Both the barrel and the disks are made up of staves, that can feature an arbitrary number of cuboid layers to describe the support, readout and sensor structures. The sensors can be described by individual sensitive and insensitive pieces to form large, complex structures such as quad modules, that have an insensitive periphery.
More details can be found in the talks at the [FCC week 2023](https://indico.cern.ch/event/1202105/timetable/#356-idea-vertex-detector-in-ke) with an update the the [MDI meeting of the 10th of July 2023](https://indico.cern.ch/event/1292318/#5-vxd-implementation-in-full-s). Once public, the constructor and the resulting IDEA vertex detector are described in the MDI note for the mid-term review.