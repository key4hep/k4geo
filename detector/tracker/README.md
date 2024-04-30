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


### o1_v02
Reimplementation of the drift chamber concept, based on a fully detailed spreadsheet. The improvements with respect the previous implementation are the following:
- The standalone simulation consumes 400MB of memory, out of which the geometry represents 150MB.
- The DCH is built out of few global parameters, thanks to the prescription given in the spreadsheet
- C++ guidelines are followed
- ASCII drawing is included to make sense of positioning of wires inside the cells
- The cell shape corresponds to a twisted tube, which is missing in ROOT, but DD4hep provides a workaround
- Because of the previous point, this subdetector can not be used by any ROOT-based application, and the shape parameters can be accessed only by DD4hep/Geant4
- Visualization of the full DCH is possible with Geant4+Qt, but not with ROOT-based applications
- The optional tag `<debugGeometry/>` build only 3 sectors of each layer, it must be used only when checking for overlaps.
- Endcap services. A dummy plate with 5% X0 is used to account for such services.
- Vessel wall is a sandwich of Carbon fiber and PE foam. The thickness of the fill material is given as a fraction of the total thickness of the wall. It is adjusted to provide 1.2%X0 radially and 5%X0 longitudinally.
- Material of field and sense wire is averaged for the sake of speedup.
- Guard wires are not implemented yet.
