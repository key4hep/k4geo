IDEA
====

IDEA_o1_v01
------------

IDEA version picked from the latest version in FCCDetectors repo

IDEA_o1_v02
------------

Based on o1_v01 but with a detailed description of the vertex detector, drift chamber, a place holder solenoid and the endplate absorber. Missing: SiWrapper, calorimeter.

IDEA_o1_v03
------------

Based on o1_v02 but replacing the drift chamber (o1, v01) for the lightweight implementation based on twisted tubes (o1,
v02). NB: production threshold and step limit physics have to be tuned for the drift chamber.

July 2024: Added a detailed version of the muon system.

The monolithic fiber dual-readout calorimeter (o1, v01) is added to the directory, but commented in the main file IDEA_o1_v03.xml for the sake of speed. Please remove the comments to include this calorimeter in the full IDEA detector. 

August 2024: Added an updated vertex detector (also with the ultra-light inner vertex option) and a more light-weight 
implementation of the silicon wrapper.

September 2024: Added detailed version of the pre-shower, based on muon system builder.

IDEA_o2_v01
------------

Second option of IDEA detector. The inner part up to the drift-chamber is identical to IDEA_o1, the dual-readout calorimeter uses the INFN capillary-tubes technology and replaces the monolithic calorimeter description. Between the drift-chamber and the dual-readout calorimeter a dual-readout crystal electromagnetic calorimeter will be placed, consequentially the preshower is removed. The muon system is identical to IDEA_o1.

October 2024: first implementation using the dual-readout capillary-tubes endcap geometry.

December 2024: Added the dual-readout capillary-tubes barrel calorimeter.
