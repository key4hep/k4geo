import os
from Gaudi.Configuration import *

from Configurables import ApplicationMgr
ApplicationMgr().EvtSel = 'None' 
ApplicationMgr().EvtMax = 1
ApplicationMgr().OutputLevel = INFO

# DD4hep geometry service
from Configurables import GeoSvc
## parse the given xml file
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
                        'FCCee_IDEA_o1_v01.xml'
                       ]
geoservice.OutputLevel = INFO 
ApplicationMgr().ExtSvc += [geoservice]

# Using material scan from k4SimGeant4: https://github.com/HEP-FCC/k4SimGeant4/tree/main/Detector/DetComponents/src
from Configurables import MaterialScan_2D
# Material scan is done from the interaction point to the end of world volume.
# In order to use other end boundary, please provide the name of a thin, e.g. cylindrical volume.
# For instance adding envelopeName="BoundaryPostCalorimetry" will perform the scan only till the end of calorimetry.
# BoundaryPostCalorimetry is defined in Detector/DetFCChhECalInclined/compact/envelopePreCalo.xml
materialservice = MaterialScan_2D("GeoDump")
materialservice.filename = "out_material_scan.root"

materialservice.angleDef = 'eta' # eta, theta, cosTheta or thetaRad
materialservice.angleBinning = 0.05
materialservice.angleMax = 3.0
materialservice.angleMin = -3.0
materialservice.nPhi = 100 # number of bins in phi

ApplicationMgr().ExtSvc += [materialservice]


