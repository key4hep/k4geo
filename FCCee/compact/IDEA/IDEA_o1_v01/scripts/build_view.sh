### Small script showing how geometry can be saved to root files, 
### that can be easily viewed using e.g https://root.cern/js/latest/
### as described in https://hep-fcc.github.io/fcc-tutorials/full-detector-simulations/Visualization/Visualization.html#detector-geometry

chmod u+x ../../../utils/dd4hep2root

xml=$1

if [[ -n "$xml" ]]; then # test to see if not empty
    ./../../../scripts/utils/dd4hep2root -c ${xml} -o fccee_idea.root
else
    echo "argument error, please provide an xml file as input argument!"
fi