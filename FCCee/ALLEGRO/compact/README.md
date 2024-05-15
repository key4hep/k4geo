ALLEGRO
========================
ALLEGRO_o1_v01: it is a liquid Noble gas based detector. This version picked from the latest version in FCCDetectors repo.

ALLEGRO_o1_v02: evolves from o1_v01, replacing the barrel ECAL and adding a detailed version drift chamber. 
This version has a constant cell size in theta for the ECAL barrel (instead of eta as in o1_v01) and now it is possible to have a different number of cells merged for each longitudinal layer.
Known caveat: the drift chamber has a larger z extent than in the IDEA detector but the wire spacing was not re-optimized. It is ok software-wise but the currently implemented design is not fully compliant with R&D considerations, will need a new drift chamber layout from the detector concept team.

ALLEGRO_o1_v03: with respect to v02 it features an ECal barrel with 11 layers and cell corners projective along phi
