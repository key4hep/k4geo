# IDEA Muon System

The IDEA detector concept includes a muon detection system and pre-shower designed using µRWELL technology. Each station consists of a large mosaic of 50 × 50 cm² µRWELL detectors. For more technical details about the new detector technology, see the [µRWELL](https://iopscience.iop.org/article/10.1088/1748-0221/10/02/P02008).

## muonSystem_v00
The first draft for the detailed version builder for the muon system, it can be used to describe both barrel and endcap parts, if you want eleminate one of the just indicate the number of layers= 0.
The code has been designed to very flexible with the modifications, where the user can choose the number of sides of the shape `numSides` (hexagon, octagon, ....etc), and automatically the builder will calculate the number and places of the copied chambers. Some of the code advantages: 
 * If the side length do not fit with an integer number of 50 × 50 cm² , the builder can make a chamber with unusual dimensions, which can fit the excess area at the end of the side.
 * The availability to make multiple layers with different inner radius and barrel length.
 * The code is very general, which can be used to describe any detector system made from repeated tiles, with having the capability to fill the gaps with unusual dimensions tiles.
 