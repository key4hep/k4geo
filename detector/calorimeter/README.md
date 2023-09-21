# Calorimeters

## ECalBarrel_NobleLiquid_InclinedTrapezoids

This sub-detectors consists of inclined trapezoidal absorbers, noble liquid sensitive media and straight readout electrodes; everything being placed inside a cryostat. A more detailed technical description of this sub-detector can be found in:
 * [Calorimeters for the FCC-hh](https://arxiv.org/abs/1912.09962)
 * [Noble liquid calorimetry for a future FCC-ee experiment](https://www.sciencedirect.com/science/article/pii/S0168900222004600)

The documentation about its usage is [here](../../doc/detector/calorimeter/ECalBarrel_NobleLiquid_InclinedTrapezoids.md).

### o1_v01
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/main/Detector/DetFCChhECalInclined/src/ECalBarrelInclined_geo.cpp).

## CaloDisks
This sub-detector makes calorimeter endcaps (original and reflected). Each endcap spans over z = 320 cm and z = 385 cm, and consist in 134 layers. Each layer is a sandwich of the following material and thickness (in cm):
* LAr, 0.2
* PCB, 0.12
* LAr, 0.2
* lArCaloSteel, 0.014
* lArCaloGlue, 0.009
* Lead, 0.104
* lArCaloGlue, 0.009
* lArCaloSteel, 0.014

### o1_v01 
Original version taken from [FCCDetectors](https://github.com/HEP-FCC/FCCDetectors/blob/70a989a6fc333610e3b1b979c3596da9c41543d8/Detector/DetFCChhCalDiscs/src/CaloEndcapDiscs_geo.cpp)

