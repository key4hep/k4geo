ALLEGRO
========================
ALLEGRO_o1_v01: it is a liquid Noble gas based detector. This version picked from the latest version in FCCDetectors repo.

ALLEGRO_o1_v02: evolves from o1_v01, replacing the barrel ECAL and adding a detailed version drift chamber.
This version has a constant cell size in theta for the ECAL barrel (instead of eta as in o1_v01) and now it is possible to have a different number of cells merged for each longitudinal layer.
Known caveat: the drift chamber has a larger z extent than in the IDEA detector but the wire spacing was not re-optimized. It is ok software-wise but the currently implemented design is not fully compliant with R&D considerations, will need a new drift chamber layout from the detector concept team.

ALLEGRO_o1_v03: with respect to v02 it features an ECal barrel with 11 layers and cell corners projective along phi.
The vertex detector and drift chamber are now taken directly from IDEA_o1_v03, this effectively updates both the vertex detector (which was taken from an old CLD version) and the drift chamber (which was corresponding to IDEA_o1_v02/DriftChamber_o1_v01.xml). The z-extent of the drift chamber is now unchanged w.r.t. the IDEA detector (2 m) since it requires optimization anyway.
Magnetic fields (solenoid + MDI) have been added.
Added "turbine-style" endcap ecal, and invoke this in the top-level xml (replacing the coneCyro geometry).  Updated to allow a more flexible geometry and calibration (details in detector/calorimeter/README.md).
Added HCalBarrel_TileCal_v02.xml which uses HCalTileBarrel_o1_v02_geo.cpp and removed unused readout BarHCal_Readout_phi. 
Added HCalEndcaps_ThreeParts_TileCal_v02.xml which uses HCalThreePartsEndcap_o1_v02_geo.cpp. Additionally, wrt v02 the readout was migrated to the theta-phi segmentation; unused readout *Readout_phi was removed; radial dimensions of layers were modified, so the outer radius of all three cylinders is the same.
For the muon tagger, switched from eta-phi to theta-phi segmentation.

ALLEGRO_o1_v04: same as v03, but material in inner part of absorber in first layer of ECAL is configurable and set by default to G10
