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
from Configurables import MaterialScan
# Material scan is done from the interaction point to the end of world volume.
# In order to use other end boundary, please provide the name of a thin, e.g. cylindrical volume.
# For instance adding envelopeName="BoundaryPostCalorimetry" will perform the scan only till the end of calorimetry.
# BoundaryPostCalorimetry is defined in Detector/DetFCChhECalInclined/compact/envelopePreCalo.xml
materialservice = MaterialScan("GeoDump")
materialservice.filename = "out_material_scan.root"
materialservice.angleBinning = 0.05
materialservice.angleMin = -3.0
materialservice.angleMax = 3.0
materialservice.nPhiTrials = 100
materialservice.angleDef = "eta" # eta or cosTheta or theta or thetaRad
ApplicationMgr().ExtSvc += [materialservice]


