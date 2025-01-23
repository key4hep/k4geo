# Calorimeters

## ECalBarrel_NobleLiquid_InclinedTrapezoids

This sub-detectors consists of inclined trapezoidal absorbers, noble liquid sensitive media and straight readout electrodes; everything being placed inside a cryostat. A more detailed technical description of this sub-detector can be found in:
 * [Calorimeters for the FCC-hh](https://arxiv.org/abs/1912.09962)
 * [Noble liquid calorimetry for a future FCC-ee experiment](https://www.sciencedirect.com/science/article/pii/S0168900222004600)

The documentation about its usage is [here](../../doc/detector/calorimeter/ECalBarrel_NobleLiquid_InclinedTrapezoids.md).

### o1_v01
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/main/Detector/DetFCChhECalInclined/src/ECalBarrelInclined_geo.cpp).

### o1_v02
New version, with module-theta based segmentation. In each layer adjacent cells along theta or module directions can be grouped together, with possibly different merging per layer. This is specified with the `mergedCells_Theta` and `mergedModules` vectors in the `segmentation` tag of the xml file. The baseline grouping in theta is by four in all layers except L1 (the strip layer). The baseline grouping in module direction is by two in all layers. The LAr gap has also been slightly adjusted to bring back the number of modules to 1536 (it was 1545 before). The segmentation class needs to know the number of modules, which is passed via the `nModules` parameter of the `segmentation` tag. To ensure that number of modules and layers (length of the mergedXXX vectors) are consistent with number of modules and layers of the detector, the xml defines `ECalBarrelNumLayers` and `ECalBarrelNumPlanes`, and the c++ file doing the detector construction checks that the number of planes and layers calculated dynamically from other parameters matches that in the xml (if not, the code will crash).

Overlaps in the LAr bath volume fixed.

## CaloDisks
This sub-detector makes calorimeter endcaps (original and reflected). It is used in ALLEGRO detector concept.

### o1_v01 
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/70a989a6fc333610e3b1b979c3596da9c41543d8/Detector/DetFCChhCalDiscs/src/CaloEndcapDiscs_geo.cpp).

## ECalEndcap_Turbine

Sub-detector for the ecal endcaps, with the absorbers and readout boards arranged in a "turbine-like" geometry.

### o1_v01
Initial implementation.  A custom segmentation that creates readout cells and constant radius and z is also included (FCCSWEndcapTurbine_k4geo).

### o1_v02
Changes wrt o1_v01: Added flexibility to configure the wheels individual (to allow for example the possibiliity of having different blade angles in each).  This is still a work in progress, so o1_v01 should be used for integrated tests.

## HCalTileBarrel
This sub-detector makes calorimeter barrel. It is used in ALLEGRO detector concept.

### o1_v01
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/70a989a6fc333610e3b1b979c3596da9c41543d8/Detector/DetFCChhHCalTile/src/HCalBarrel_geo.cpp). 

### o1_v02
Changes wrt o1_v01: Added extension (LayeredCalorimeterData) to store radial layer radii and dimensions. Added several checks for the geometry building and made small changes in the code to improve readibility.

## HCalThreePartsEndcap
This sub-detector makes calorimeter endcaps. Each endcap is made up by three cylindrical pieces with different thickness and inner radius, but same outer radius. It is used in ALLEGRO detector concept.

### o1_v01
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/70a989a6fc333610e3b1b979c3596da9c41543d8/Detector/DetFCCeeHCalTile/src/HCalThreePartsEndcap_geo.cpp#L4). 

### o1_v02
Changes wrt o1_v01: Added extension (LayeredCalorimeterData) to store radial layer radii and dimensions. To make this work, the whole code had to be reshuffled, but the way how the geometry and individual volumes are built remains the same as in o1_v01. 

## dual-readout

### o1_v01
This sub-detector makes full 4-pi monolithic fiber dual-readout calorimeter.
Inside the single tower (trapezoidal copper absorber), two types of optical fibers (Cherenkov and scintillation) are implemented. The readout (SiPM) is attached at the rear side of the tower. The tower is repeated in both eta and phi direction to cover both barrel and endcap region.

## dual-readout-tubes

### o1_v01
This folder containes the subdetectors (endcap + barrel) to make a full 4-pi fiber dual-readout calorimeter exploiting the INFN capillary-tubes technology. Each trapezoidal tower is constructed with brass capillary-tubes housing optical fibers (Cherenkov and scintillating). Endcap and barrel calorimeters are implemented ad separate subdetectors.
