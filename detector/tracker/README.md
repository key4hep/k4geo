# Vertex detectors

## VertexBarrel_detailed and VertexDisks_detailed description
These two detector constructors were derived from ZPlanarTracker.cpp and from VertexEndcap_o1_v06.cpp and adapted to fit the needs of the [IDEA vertex detector engineering design](https://indico.cern.ch/event/1202105/timetable/#242-mechanical-integration-of).

Both the barrel and the disks are made up of staves, that can feature an arbitrary number of cuboid layers to describe the support, readout and sensor structures. The sensors can be described by individual sensitive and insensitive pieces to form large, complex structures such as quad modules, that have an insensitive periphery.
More details can be found in the talks at the [FCC week 2023](https://indico.cern.ch/event/1202105/timetable/#356-idea-vertex-detector-in-ke) with an update the the [MDI meeting of the 10th of July 2023](https://indico.cern.ch/event/1292318/#5-vxd-implementation-in-full-s). Once public, the constructor and the resulting IDEA vertex detector are described in the MDI note for the mid-term review.

# Trackers

## parametrised_SimplifiedDriftChamber
This sub-detector implements a simplied drift chamber (there is no actual wire). It is used in ALLEGRO and IDEA detector concepts.

### o1_v01 
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/70a989a6fc333610e3b1b979c3596da9c41543d8/Detector/DetFCChhCalDiscs/src/CaloEndcapDiscs_geo.cpp). 

## DriftChamber
Detailed (i.e. including wires) description of a unique-volume, high granularity, fully stereo, low mass, cylindrical drift chamber. 
Originally designed for the IDEA detector, it will be used also for the ALLEGRO concept with different dimensions. 
A detailed description of the detector it tries to model can be found [here](https://indico.cern.ch/event/932973/contributions/4041314/attachments/2139657/3664808/primavera_FCCworkshop_2020.pdf) or [here](https://indico.cern.ch/event/1283129/contributions/5476695/attachments/2682660/4654170/DeFilippis_DCH_CC.pdf).

### o1_v01
First version with all wires, walls, etc. 
This version do not rely on any external segmentation, all volumes are defined in the detector builder. 
The sensitive cell volume is a rotated (to get the stereo angle) tube segment that are put inside a hyperboloid for each layer. 
The hyperboloid is theoretically not needed, but it help a lot to have a well balanced volume tree.
This first version can clearly be improved and still suffers from some overlap but it already enables further technical developments (digitization, tracking, particle flow, ...).
