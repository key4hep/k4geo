outputFolder=IDEA_outerTracker_endcap
outputFolder=IDEA_conicalSupport
outputFolder=IDEA_vertex

k4run scripts/utils/material_scan.py 
python scripts/utils/material_plots.py --fname out_material_scan.root --etaMax 3 --x0max 0.06 # eta of 3 is equal to theta of 0.099, so 99 mrad
mkdir $outputFolder
mv out_material_scan.root *.pdf $outputFolder
