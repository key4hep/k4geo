# µRWELL Based Muon System

The IDEA detector concept includes a muon detection system and pre-shower designed using µRWELL technology. Each station consists of a large mosaic of 50 × 50 cm² µRWELL detectors. For more technical details about the new detector technology, see the [µRWELL](https://iopscience.iop.org/article/10.1088/1748-0221/10/02/P02008).

## muonSystemMuRWELL_o1_v01.cpp
The first version of the detailed muon system driver, it can be used to describe both barrel and endcap, if you want to eleminate one of them, just set the number of layers= 0.
The code has been designed to be very flexible, where the user can choose the number of sides in the R-Phi plane, `numSides` (hexagon, octagon, etc), and the detector builder will automatically calculate the number and places of the copied chambers. Some of the code advantages: 
 * If the side length do not fit with an integer number of 50 × 50 cm² , the builder can make a chamber with unusual dimensions, which can fit the excess area at the end of the side.
 * The availability to make multiple layers with different inner radius and barrel length.
 * The code is very general, it can be used to describe any detector system made from repeated tiles (e.g. pre-shower) and has the capability to fill the gaps with unusual dimensions tiles.
