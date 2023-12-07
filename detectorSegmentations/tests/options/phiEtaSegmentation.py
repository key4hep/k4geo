from Gaudi.Configuration import *
from Configurables import ApplicationMgr, HepMCDumper


from Configurables import MomentumRangeParticleGun
from GaudiKernel import PhysicalConstants as constants
guntool = MomentumRangeParticleGun()
guntool.ThetaMin = 0 
guntool.ThetaMax = 2 * constants.pi 
guntool.PdgCodes = [11]
guntool.MomentumMin = 500
guntool.MomentumMax = 1000
from Configurables import FlatSmearVertex
vertexsmeartool = FlatSmearVertex()
vertexsmeartool.xVertexMin = -25.
vertexsmeartool.xVertexMax = 25.
vertexsmeartool.yVertexMin = -25.
vertexsmeartool.yVertexMax = 25.
vertexsmeartool.zVertexMin = -25.
vertexsmeartool.zVertexMax = 25.
from Configurables import GenAlg
gen = GenAlg()
gen.SignalProvider=guntool
gen.VertexSmearingTool=vertexsmeartool
gen.hepmc.Path = "hepmc"

from Configurables import HepMCToEDMConverter
hepmc_converter = HepMCToEDMConverter("Converter")
hepmc_converter.hepmc.Path="hepmc"
hepmc_converter.genparticles.Path="allGenParticles"
hepmc_converter.genvertices.Path="allGenVertices"

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors=['file:Test/TestGeometry/data/Barrel_testCaloSD_phieta.xml'])

from Configurables import SimG4Svc
geantservice = SimG4Svc("SimG4Svc")

from Configurables import SimG4Alg, SimG4SaveCalHits, InspectHitsCollectionsTool
inspecttool = InspectHitsCollectionsTool("inspect", readoutNames=["ECalHits"], OutputLevel = DEBUG)
savecaltool = SimG4SaveCalHits("saveECalHits", readoutNames = ["ECalHits"], OutputLevel = DEBUG)
savecaltool.positionedCaloHits.Path = "positionedCaloHits"
savecaltool.caloHits.Path = "caloHits"
geantsim = SimG4Alg("SimG4Alg", outputs= ["SimG4SaveCalHits/saveECalHits","InspectHitsCollectionsTool/inspect"])

from Configurables import FCCDataSvc, PodioOutput
podiosvc = FCCDataSvc("EventDataSvc")
out = PodioOutput("out", filename="test_phiEtaSegmentation.root")
out.outputCommands = ["keep *"]

ApplicationMgr(EvtSel='NONE',
               EvtMax=10,
               TopAlg=[gen, hepmc_converter, geantsim, out],
               ExtSvc = [podiosvc, geoservice, geantservice],
               )
