ALLEGRO
========================
ALLEGRO_o1_v01: it is a liquid Noble gas based detector. This version picked from the latest version in FCCDetectors repo.

ALLEGRO_o1_v02: evolves from o1_v01, replacing the barrel ECAL and adding a detailed version drift chamber.
This version has a constant cell size in theta for the ECAL barrel (instead of eta as in o1_v01) and now it is possible to have a different number of cells merged for each longitudinal layer.
Known caveat: the drift chamber has a larger z extent than in the IDEA detector but the wire spacing was not re-optimized. It is ok software-wise but the currently implemented design is not fully compliant with R&D considerations, will need a new drift chamber layout from the detector concept team.

ALLEGRO_o1_v03:
- The vertex detector and drift chamber are now taken directly from IDEA_o1_v03, this effectively updates both the vertex detector (which was taken from an old CLD version) and the drift chamber (which was corresponding to IDEA_o1_v02/DriftChamber_o1_v01.xml).
- The z-extent of the drift chamber is now unchanged w.r.t. the IDEA detector (2 m) since it requires optimization anyway.
- Magnetic fields (solenoid + MDI) have been added.
- ECAL barrel: with respect to v02 it features 11 layers and cell corners projective along phi. Moreover, material in inner part of absorber in first layer of ECAL is configurable (set by default to LAr for backward compatibility, but ideally should be G10)
- Added "turbine-style" endcap ecal, and invoke this in the top-level xml (replacing the coneCyro geometry).
- Updated to allow a more flexible geometry and calibration (details in detector/calorimeter/README.md).
- Added HCalBarrel_TileCal_v03.xml which uses HCalTileBarrel_o1_v02_geo.cpp and removed unused readout BarHCal_Readout_phi. 
- Added HCalEndcaps_ThreeParts_TileCal_v03.xml which uses HCalThreePartsEndcap_o1_v02_geo.cpp. Additionally, wrt v02 the readout was migrated to the theta-phi segmentation; unused readout *Readout_phi was removed; radial dimensions of layers were modified, so the outer radius of all three cylinders is the same.
- Row-phi segmentation is added for HCalBarrel_TileCal_v03.xml and HCalEndcaps_ThreeParts_TileCal_v03.xml.
- Birks constant value is set for Polystyrene scintillator used by HCal. This fixed the abnormal response to hadrons that was observed when migrated from k4SimGeant4 to DDSim.
- For the muon tagger, switched from eta-phi to theta-phi segmentation. Outer R set to 5m as in initial conceptual design (to be replaced in the future by more realistic detector).
- February 2026: The vertex detector and silicon wrapper are now taken from IDEA_o1_v04.

ALLEGRO_o2_v01:
- The drift chamber is replaced with a straw tube tracker.
- The straw tube tracker uses thin wall mylar (12um) straw tubes.
- While the straw tracker will inherently have a higher material budget (worse spatial resolution) and more dead volume (worse dN/dx resolution, provided no compensation is made with pressure) than the drift chamber, its modular design may provide a significant advantage.
- This version of ALLEGRO allows for direct comparison between the drift chamber and the straw tracker, as compared to o1_v03 the only changes are (i) the swap of drift chamber and straw and (ii) elongating the tracking volume, including the silicon wrapper.
- The elongated wrapper is +660mm longer in both the positive and negative z directions with a 60mm buffer between the ends of the straw and the silicon wrapper endcap, to be filled later with straw RO/HV and support.
- February 2026: The vertex detector and silicon wrapper are now taken from IDEA_o1_v04.
