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
