# Vertex detectors

## VertexBarrel_detailed and VertexEndcap_detailed description

### o1_v01
These two detector constructors were derived from ZPlanarTracker.cpp and from VertexEndcap_o1_v06.cpp and adapted to fit the needs of the [IDEA vertex detector engineering design](https://indico.cern.ch/event/1202105/timetable/#242-mechanical-integration-of).

Both the barrel and the disks are made up of staves, that can feature an arbitrary number of cuboid layers to describe the support, readout and sensor structures. The sensors can be described by individual sensitive and insensitive pieces to form large, complex structures such as quad modules, that have an insensitive periphery.
More details can be found in the talks at the [FCC week 2023](https://indico.cern.ch/event/1202105/timetable/#356-idea-vertex-detector-in-ke) with an update the the [MDI meeting of the 10th of July 2023](https://indico.cern.ch/event/1292318/#5-vxd-implementation-in-full-s). The constructor and the resulting IDEA vertex detector are described in the MDI note for the mid-term review (not public).

### o1_v02
These versions for the vertex barrel and disks fix previous issues with overlaps and enable to define in the xml file mother volumes to speed up the simulation by having less volumes per node in the volume hierarchy (especially relevant for the large-area silicon wrapper that also uses these constructors).
o1_v02 of the barrel also enables the use of curved support and sensor volumes needed to describe the ultra-light vertex detector proposed in the 2024 [FCC Physics Workshop](https://indico.cern.ch/event/1307378/timetable/?view=standard#84-vertex-detector-and-silicon) and [FCC Week](https://indico.cern.ch/event/1298458/timetable/#15-optimization-of-si-tracking).

### o1_v03
Barrel: One can now have a partially insensitive sensor using 'sensor_insensitive_thickness_below' and 'sensor_insensitive_thickness_above', which enables the detailed description of sensors that are not fully sensitive in thickness (e.g. sensors in TPSCo 65 nm process for curved vertex detectors). Curved sensors can now also be approximated by trapezoidal volumes with the 'nsegments' option. This allows to perform track reconstruction with conformal tracking, which currently doesn't support curved sensors (=curved sensitive surfaces) yet. 
For both detector builders, volumes with 0 thickness are ignored now and the naming scheme is more consistent. Some of these changes have been presented in this [FCC Full Sim Meeting](https://indico.cern.ch/event/1649968/#56-update-on-vtx-and-silicon-w).

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

## Straw Tracker

The straw tracker concept uses O(100k) thin wall straw tubes to track particles.
It likely must be paired with a silicon wrapper to achieve the reqiured spatial resolution, similar to the drift chamber.
While the straw will always have more material budget and more dead volume than a drift chamber, the performance degredation is slight and a possible advantage is the modular construction.
The straw tracker has been discussed at FCC meetings [here](https://indico.fnal.gov/event/67484/contributions/313363/attachments/187186/258081/straw_USFCC.pdf) and [here](https://indico.cern.ch/event/1516157/contributions/6444547/attachments/3063289/5418866/straw_TrackerWorkshop_BNL_05082025.pptx%20(2).pdf).
There was also a mini workshop on tracking, focusing on fabrication of straw and drift chambers, including lessons learned from other experiments, which can be found [here](https://indico.cern.ch/event/1408681/).  The detector builder is described in detail [here](https://indico.cern.ch/event/1551837/#71-straw-tube-tracker-implemen).

### o1_v01
Initial design for a tracker based on thin-wall straw tubes.  This will likely be part of the ALLEGRO concept, but can similarly be swapped with the IDEA drift chamber for apples-to-apples comparison studies.
 - Flexible implementation allows for user to customize n tube layers; thickness of wire, mylar, coating, and gas volumes; size of gaps between mutlilayers; etc.
 - Low logical volume count of 61 for a realistic straw concept.
 - 200MB memory usage for straw detector.
 - Stereo angle is implemented, is +/- 2 degrees in ALLEGRO_o2_v01
 - No consideration yet to endcap.

## TPC

### TPC10
Scalable TPC model, intended to be used with plugins/TPCSDAction to combine G4 steps into SimTrackerHits.
TPCSDaction has some ingrained assumptions that G4 steps are passed to it in a logical way, i.e. in the order that they are produced.