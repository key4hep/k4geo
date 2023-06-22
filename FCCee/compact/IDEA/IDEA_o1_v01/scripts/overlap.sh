#!/bin/bash

### Script to check for overlaps, as described in 
### https://hep-fcc.github.io/fcc-tutorials/full-detector-simulations/Geometry/Geometry.html#overlap-checking

# IDEA
nohup ddsim --compactFile ./FCCee_IDEA_o1_v01.xml --runType run --macroFile scripts/utils/overlap.mac &> overlap.log &
