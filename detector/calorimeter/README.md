# Calorimeters

## ECalBarrel_NobleLiquid_InclinedTrapezoids

This sub-detectors consists of inclined trapezoidal absorbers, noble liquid sensitive media and straight readout electrodes; everything being placed inside a cryostat. A more detailed technical description of this sub-detector can be found in:
 * [Calorimeters for the FCC-hh](https://arxiv.org/abs/1912.09962)
 * [Noble liquid calorimetry for a future FCC-ee experiment](https://www.sciencedirect.com/science/article/pii/S0168900222004600)

The documentation about its usage is [here](../../doc/detector/calorimeter/ECalBarrel_NobleLiquid_InclinedTrapezoids.md).

### o1_v01
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/main/Detector/DetFCChhECalInclined/src/ECalBarrelInclined_geo.cpp).

### o1_v02
New version adapted to the theta segmentation with the possibility to have different number of cell merged per layer. The main difference is that now one has to set ECalBarrelNumLayers in the xml while before it was just dynamically computed based on other xml parameters (also, to avoid silent mistakes, the number from the xml and the one computed dynamically must match).

## CaloDisks
This sub-detector makes calorimeter endcaps (original and reflected). It is used in ALLEGRO detector concept.

### o1_v01 
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/70a989a6fc333610e3b1b979c3596da9c41543d8/Detector/DetFCChhCalDiscs/src/CaloEndcapDiscs_geo.cpp).

## HCalTileBarrel
This sub-detector makes calorimeter barrel. It is used in ALLEGRO detector concept.

### o1_v01
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/70a989a6fc333610e3b1b979c3596da9c41543d8/Detector/DetFCChhHCalTile/src/HCalBarrel_geo.cpp). 

## HCalThreePartsEndcap
This sub-detector makes calorimeter endcaps. Each endcap is made up by three cylindrical pieces with different thickness and inner radius, but same outer radius. It is used in ALLEGRO detector concept.

### o1_v01
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/70a989a6fc333610e3b1b979c3596da9c41543d8/Detector/DetFCCeeHCalTile/src/HCalThreePartsEndcap_geo.cpp#L4). 



