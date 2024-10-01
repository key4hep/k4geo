"""
This script must be called with k4run: 'k4run material_scan_2D.py --{argument} {value}'.
The output files are saved in 'data/{outputDir}/{outputFileBase}.root'.
If no outputDir is specified, it will be 'data/{outputFileBase}.root'.
"""

from os import environ, fspath
from pathlib import Path

from Configurables import ApplicationMgr
from Gaudi.Configuration import *
from k4FWCore.parseArgs import parser

ApplicationMgr().EvtSel = "None"
ApplicationMgr().EvtMax = 1
ApplicationMgr().OutputLevel = INFO

# DD4hep geometry service
from Configurables import GeoSvc

parser.add_argument(
    "--compactFile",
    help="Compact detector file to use",
    type=str,
    default=fspath(
        Path(environ["k4geo_DIR"])
        / "ILD"
        / "compact"
        / "ILD_sl5_v02"
        / "ILD_l5_v02.xml"
    ),
)
parser.add_argument(
    "--outputFileBase",
    help="Base name of all the produced output files",
    default="out_material_scan",
)
parser.add_argument(
    "--angleDef",
    help="angle definition to use: eta, theta, cosTheta or thetaRad, default: eta",
    choices=["eta", "theta", "cosTheta", "thetaRad"],
    default="eta",
)
parser.add_argument(
    "--outputDir",
    "-o",
    type=str,
    default="",
    help="Directory to store the output file in",
)

reco_args = parser.parse_known_args()[0]
compact_file = reco_args.compactFile
angle_def = reco_args.angleDef
output_dir = "data" / Path(reco_args.outputDir)
output_dir.mkdir(
    parents=True, exist_ok=True
)  # Create the directory if it doesn't exist

## parse the given xml file
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [compact_file]
geoservice.OutputLevel = INFO
ApplicationMgr().ExtSvc += [geoservice]

# Using material scan from k4SimGeant4: https://github.com/HEP-FCC/k4SimGeant4/tree/main/Detector/DetComponents/src
from Configurables import MaterialScan_2D_genericAngle

# Material scan is done from the interaction point to the end of world volume.
# In order to use other end boundary, please provide the name of a thin, e.g. cylindrical volume.
# For instance adding envelopeName="BoundaryPostCalorimetry" will perform the scan only till the end of calorimetry.
# BoundaryPostCalorimetry is defined in Detector/DetFCChhECalInclined/compact/envelopePreCalo.xml
materialservice = MaterialScan_2D_genericAngle("GeoDump")
materialservice.filename = fspath(
    output_dir / Path(reco_args.outputFileBase).with_suffix(".root")
)

materialservice.angleDef = angle_def  # eta, theta, cosTheta or thetaRad
materialservice.angleBinning = 0.05
materialservice.angleMax = 3.0
materialservice.angleMin = -3.0
materialservice.nPhi = 100  # number of bins in phi

ApplicationMgr().ExtSvc += [materialservice]
