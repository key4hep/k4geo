# Calorimeters

## ECalBarrel_NobleLiquid_InclinedTrapezoids

This sub-detectors consists of inclined trapezoidal absorbers, noble liquid sensitive media and straight readout electrodes; everything being placed inside a cryostat. A more detailed technical description of this sub-detector can be found in:
 * [Calorimeters for the FCC-hh](https://arxiv.org/abs/1912.09962)
 * [Noble liquid calorimetry for a future FCC-ee experiment](https://www.sciencedirect.com/science/article/pii/S0168900222004600)

The documentation about its usage is [here](../../doc/detector/calorimeter/ECalBarrel_NobleLiquid_InclinedTrapezoids.md).

### o1_v01
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/main/Detector/DetFCChhECalInclined/src/ECalBarrelInclined_geo.cpp).

## CaloDisks
This sub-detector makes calorimeter endcaps (original and reflected). It is used in ALLEGRO detector concept.

### o1_v01 
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/70a989a6fc333610e3b1b979c3596da9c41543d8/Detector/DetFCChhCalDiscs/src/CaloEndcapDiscs_geo.cpp).

Each endcap spans over z = 320 cm and z = 385 cm, with an inner radius of 35 cm and outer radius of 290 cm. Each one consist in 134 layers, each layer is a sandwich of the following material (and thickness):
* LAr (0.2 cm)
* PCB (0.12 cm)
* LAr (0.2 cm)
* lArCaloSteel (0.014 cm)
* lArCaloGlue (0.009 cm)
* Lead (0.104 cm)
* lArCaloGlue (0.009 cm)
* lArCaloSteel (0.014 cm)

## HCalTileBarrel
This sub-detector makes calorimeter barrel. It is used in ALLEGRO detector concept.

### o1_v01
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/70a989a6fc333610e3b1b979c3596da9c41543d8/Detector/DetFCChhHCalTile/src/HCalBarrel_geo.cpp). It is made up by 13 layers, each one is a sandwich of Iron (5 cm) and Air (5 cm). It spans over a radial distance between 280 and 450 cm, and it has a half length of 280 cm.

## HCalThreePartsEndcap
This sub-detector makes calorimeter endcaps. Each endcap is made up by three cylindrical pieces with different thickness and inner radius, but same outer radius. It is used in ALLEGRO detector concept.

### o1_v01
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/70a989a6fc333610e3b1b979c3596da9c41543d8/Detector/DetFCCeeHCalTile/src/HCalThreePartsEndcap_geo.cpp#L4). 



