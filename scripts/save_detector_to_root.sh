### Small script showing how geometry can be saved to root files, 
### that can be easily viewed using e.g https://root.cern/js/latest/
### as described in https://hep-fcc.github.io/fcc-tutorials/full-detector-simulations/Visualization/Visualization.html#detector-geometry


script_folder="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )" # Get the location of this script, to execute correctly other commands
xml=$1 # Give xml detector file as input argument!
outputName="${2:-detector}" # Name of output file, saved as $outputName_dd4hep2root.root

echo "Running on xml file $xml"

if [[ -n "$xml" ]]; then # test to see if not empty
    python ${script_folder}/../utils/dd4hep2root.py -c ${xml} -o ${outputName}_dd4hep2root.root
else
    echo "argument error, please provide an xml file as input argument!"
fi