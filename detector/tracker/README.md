# VertexBarrel_detailed and VertexDisks_detailed description
These two detector constructors were derived from ZPlanarTracker.cpp and from VertexEndcap_o1_v06.cpp and adapted to fit the needs of the [IDEA vertex detector engineering design](https://indico.cern.ch/event/1202105/timetable/#242-mechanical-integration-of).

Both the barrel and the disks are made up of staves, that can feature an arbitrary number of cuboid layers to describe the support, readout and sensor structures. The sensors can be described by individual sensitive and insensitive pieces to form large, complex structures such as quad modules, that have an insensitive periphery.
More details can be found in the talks at the [FCC week 2023](https://indico.cern.ch/event/1202105/timetable/#356-idea-vertex-detector-in-ke) with an update the the [MDI meeting of the 10th of July 2023](https://indico.cern.ch/event/1292318/#5-vxd-implementation-in-full-s). Once public, the constructor and the resulting IDEA vertex detector are described in the MDI note for the mid-term review.

# Trackers

## parametrised_SimplifiedDriftChamber
This sub-detector implements a simplied drift chamber. It is used in ALLEGRO detector concept.

### o1_v01 
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/70a989a6fc333610e3b1b979c3596da9c41543d8/Detector/DetFCChhCalDiscs/src/CaloEndcapDiscs_geo.cpp). 

