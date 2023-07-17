#!/bin/bash

### Script to check for overlaps, as described in 
### https://hep-fcc.github.io/fcc-tutorials/full-detector-simulations/Geometry/Geometry.html#overlap-checking

xml=$1
    if [[ -n "$xml" ]]; then # test to see if not empty
    nohup ddsim --compactFile ${xml} --runType run --macroFile ../../../utils/overlap.mac &> overlap.log &
else
    echo "argument error, please provide an xml file as input argument!"
fi