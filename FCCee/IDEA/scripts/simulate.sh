### Use ddsim to make a simple particle gun test
### like in the FCC tutorial https://hep-fcc.github.io/fcc-tutorials/full-detector-simulations/Geometry/Geometry.html#modify-an-existing-xml-file

xml=$1

if [[ -n "$xml" ]]; then # test to see if not empty
    ddsim --compactFile ${xml} \                                                                                                                                                                                                                                                                                                                                 
      --enableGun \
      --gun.distribution uniform \
      --gun.energy "10*GeV" \
      --gun.particle mu- \
      --numberOfEvents 500 \
      --outputFile Step1_edm4hep.root
else
    echo "argument error, please provide an xml file as input argument!"
fi