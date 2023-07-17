xml=$1

if [[ -n "$xml" ]]; then # test to see if not empty
    outputFolder=IDEA_vertex

    # General settings
    etaBinning=0.01
    thetaBinning=0.25
    cosThetaBinning=0.005
    nPhiBins=400
    cosThetaMin=0.0
    thetaMax=90
    extra=""

    # # for disks only
    # cosThetaMin=0.6
    # thetaMax=40
    # cosThetaBinning=0.001
    # thetaBinning=0.1
    # extra="--x0max 1.4" # for individual disks only

    # # for inner barrel only
    # extra="--x0max 2.4"

    # # for outer barrel only
    # extra="--x0max 0.9"

    # # for whole outer barrel
    # extra="--x0max 1.8"


    ### Running the scans

    # 2D eta
    k4run scripts/utils/material_scan_2D.py --GeoSvc.detector ${xml} --GeoDump.filename ${outputFolder}_scan_2D_eta.root --GeoDump.angleDef eta --GeoDump.angleBinning ${etaBinning} --GeoDump.angleMin -3.0 --GeoDump.angleMax 3.0 --GeoDump.nPhi ${nPhiBins} 

    # 2D theta
    k4run scripts/utils/material_scan_2D.py --GeoSvc.detector ${xml} --GeoDump.filename ${outputFolder}_scan_2D_theta.root --GeoDump.angleDef theta --GeoDump.angleBinning ${thetaBinning} --GeoDump.angleMin 0.0 --GeoDump.angleMax $thetaMax --GeoDump.nPhi ${nPhiBins} 

    # 2D cosTheta
    k4run scripts/utils/material_scan_2D.py --GeoSvc.detector ${xml} --GeoDump.filename ${outputFolder}_scan_2D_cosTheta.root --GeoDump.angleDef cosTheta --GeoDump.angleBinning ${cosThetaBinning} --GeoDump.angleMin $cosThetaMin --GeoDump.angleMax 1.0 --GeoDump.nPhi ${nPhiBins} 

    # 1D cosTheta
    k4run scripts/utils/material_scan.py --GeoSvc.detector ${xml} --GeoDump.filename ${outputFolder}_scan_1D_cosTheta.root --GeoDump.angleDef cosTheta --GeoDump.angleBinning ${cosThetaBinning} --GeoDump.angleMin $cosThetaMin --GeoDump.angleMax 1.0 --GeoDump.nPhiTrials 100

    # 1D eta
    k4run scripts/utils/material_scan.py --GeoSvc.detector ${xml} --GeoDump.filename ${outputFolder}_scan_1D_eta.root --GeoDump.angleDef eta --GeoDump.angleBinning ${etaBinning} --GeoDump.angleMin -3.0 --GeoDump.angleMax 3.0 --GeoDump.nPhiTrials 100

    # 1D theta
    k4run scripts/utils/material_scan.py --GeoSvc.detector ${xml} --GeoDump.filename ${outputFolder}_scan_1D_theta.root --GeoDump.angleDef theta --GeoDump.angleBinning ${thetaBinning} --GeoDump.angleMin 0.0 --GeoDump.angleMax $thetaMax --GeoDump.nPhiTrials 100



    ### Plotting (note: Currently only the stable Key4hep stack can run this, in the nightlies there's a python problem - Armin Ilg, 14.07.2023)

    # 2D eta
    python scripts/utils/material_plots_2D.py --fname ${outputFolder}_scan_2D_eta.root --angleMin -3.0 --angleMax 3.0 --angleDef eta --angleBinning ${etaBinning} --nPhiBins ${nPhiBins}  
    mkdir ${outputFolder}_2D_eta
    mv *.pdf *.png x0.root lambda.root depth.root ${outputFolder}_2D_eta

    # 2D theta
    python scripts/utils/material_plots_2D.py --fname ${outputFolder}_scan_2D_theta.root --angleMin 0.0 --angleMax $thetaMax --angleDef theta --angleBinning ${thetaBinning} --nPhiBins ${nPhiBins} 
    mkdir ${outputFolder}_2D_theta
    mv *.pdf *.png x0.root lambda.root depth.root ${outputFolder}_2D_theta

    # 1D eta
    python scripts/utils/material_plots.py --fname ${outputFolder}_scan_1D_eta.root --angleMax 3.0 --angleMin -3.0 --angleBinning ${etaBinning} --angleDef eta
    mkdir ${outputFolder}_1D_eta
    mv *.pdf *.png x0.root lambda.root depth.root  ${outputFolder}_1D_eta

    # 1D theta
    python scripts/utils/material_plots.py --fname ${outputFolder}_scan_1D_theta.root --angleMin 0.0 --angleMax $thetaMax --angleDef theta --angleBinning ${thetaBinning}
    mkdir ${outputFolder}_1D_theta
    mv *.pdf *.png x0.root lambda.root depth.root ${outputFolder}_1D_theta

    # 2D cosTheta
    python scripts/utils/material_plots_2D.py --fname ${outputFolder}_scan_2D_cosTheta.root --angleMin $cosThetaMin --angleMax 1.0 --angleDef cosTheta --angleBinning ${cosThetaBinning} --nPhiBins ${nPhiBins} ${extra}
    mkdir ${outputFolder}_2D_cosTheta
    mv *.pdf *.png x0.root lambda.root depth.root ${outputFolder}_2D_cosTheta

    # 1D cosTheta
    python scripts/utils/material_plots.py --fname ${outputFolder}_scan_1D_cosTheta.root --angleMin $cosThetaMin --angleMax 1.0 --angleDef cosTheta --angleBinning ${cosThetaBinning} ${extra}
    mkdir ${outputFolder}_1D_cosTheta
    mv *.pdf *.png x0.root lambda.root depth.root ${outputFolder}_1D_cosTheta

else
    echo "argument error, please provide an xml file as input argument!"
fi
