

import os
from Gaudi.Configuration import *

from GaudiKernel.SystemOfUnits import MeV, GeV

# Electron momentum in GeV
momentum = 10
# Theta min and max in degrees
thetaMin = 0.
thetaMax = 50.

# Data service
#from Configurables import FCCDataSvc
#podioevent  = FCCDataSvc("EventDataSvc")

from Configurables import k4DataSvc, PodioInput
evtsvc = k4DataSvc('EventDataSvc')
evtsvc.input = "/data1/varnes/FCC/key4hep3/k4geo/ALLEGRO_sim_edm4hep_mu_10GeV.root"
podioevent  = k4DataSvc("EventDataSvc")

inp = PodioInput('InputReader')
inp.collections = [
  'EventHeader',
  'MCParticles',
  'ECalEndcapTurbine',
]

# DD4hep geometry service
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc",
                    OutputLevel = INFO)

path_to_detector = os.environ.get("FCCDETECTORS", "")
detectors_to_use=[ 
    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v02/ALLEGRO_o1_v02.xml',
    #'Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectEmptyMaster.xml',
    #'Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectMaster.xml',
    #'Detector/DetFCCeeECalInclined/compact/original_FCCee_ECalBarrel_calibration.xml',
    #'Detector/DetFCCeeECalInclined/compact/FCCee_ECalBarrel_calibration.xml',
    ]
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]

ecalEndcapHitsName = "ECalEndcapTurbine"
from Configurables import SamplingFractionInLayers
hist = SamplingFractionInLayers("hists",
                                 energyAxis = momentum,
                                 readoutName = "ECalEndcapTurbine",
                                 layerFieldName = "layer",
                                 activeFieldName = "type",
                                 activeFieldValue = 0,
                                 numLayers = 30,
                                 OutputLevel = INFO)
hist.deposits.Path = ecalEndcapHitsName

THistSvc().Output = ["rec DATAFILE='histSF_fccee_turbine.root' TYP='ROOT' OPT='RECREATE'"]
THistSvc().PrintAll=True
THistSvc().AutoSave=True
THistSvc().AutoFlush=False
THistSvc().OutputLevel=INFO

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
#geantsim.AuditExecute = True
hist.AuditExecute = True

from Configurables import PodioOutput
### PODIO algorithm
out = PodioOutput("out",OutputLevel=INFO)
out.outputCommands = ["drop *"]
out.filename = "fccee_samplingFraction_inclinedEcal.root"

from Configurables import EventCounter
event_counter = EventCounter('event_counter')
event_counter.Frequency = 10

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [inp, event_counter, hist, out],
                EvtSel = 'NONE',
                EvtMax = -1,
                # order is important, as GeoSvc is needed by G4SimSvc
                ExtSvc = [geoservice, podioevent, audsvc],
                OutputLevel = INFO,
                StopOnSignal = True
)
