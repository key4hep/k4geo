#!/bin/bash

### Script to check for overlaps, as described in 
### https://hep-fcc.github.io/fcc-tutorials/full-detector-simulations/Geometry/Geometry.html#overlap-checking

script_folder="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )" # Get the location of this script, to execute correctly other commands
xml=$1 # Give xml detector file as input argument!

echo "Running on xml file $xml"

if [[ -n "$xml" ]]; then # test to see if not empty
    nohup ddsim --compactFile ${xml} --runType run --macroFile ${script_folder}/../utils/overlap.mac &> overlap.log &
else
    echo "argument error, please provide an xml file as input argument!"
fi
