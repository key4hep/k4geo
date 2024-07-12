#!/bin/bash


Marlin --global.LCIOInputFiles="muonGun_tracker10TeV_p001-100_theta89_10k_sim.slcio"  --Output_REC.LCIOOutputFile="muonGun_tracker10TeV_p001-100_theta89_10k_reco.slcio" actsseedckf_steer.xml
