from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV, keV, eV
from math import pi

SIM = DD4hepSimulation()
SIM.runType = "batch"

SIM.printLevel = 5
SIM.output.geometry = 7
SIM.output.inputStage = 7
SIM.output.kernel = 5
SIM.output.part = 7
SIM.output.random = 7

SIM.compactFile = ['FCCee/IDEA/compact/IDEA_o2_v01/SCEPCal.xml']
SIM.macroFile = ""

opticalPhysics = True

# SIM.inputFiles = ['examples/wzp6_ee_ZZ_test_ecm240_1k.stdhep']
# SIM.outputFile = 'examples/wzp6_ee_ZZ_test_ecm240_n1_cut0_BEonly.root'
SIM.numberOfEvents = 1
SIM.skipNEvents = 0

SIM.gun.multiplicity = 1
# SIM.gun.position = (0, 0, 0)  #(0, 0, -50.0*mm)  #(0, 0, 50.0*mm)  #(0, 0, 0)
# SIM.gun.direction = (1, 1, 0)
# SIM.gun.isotrop = False
# SIM.gun.distribution = 'uniform'
# SIM.gun.energy = 10*GeV
# SIM.gun.particle = "gamma"
SIM.gun.momentumMin = 10.0*MeV-10*keV  #10.00000*GeV
SIM.gun.momentumMax = 10.0*MeV+10*keV  #10.00001*GeV
# SIM.gun.phiMin = 10*pi/180.0
# SIM.gun.phiMax = 10*pi/180.0
# SIM.gun.thetaMin = (90-10)*pi/180.0  

def setupEDM4hepOutputDR(dd4hepSimulation):
     from DDG4 import EventAction, Kernel
     dd = dd4hepSimulation
     evt_edm4hep = EventAction(Kernel(), 'Geant4Output2EDM4hep_DRC/' + dd.outputFile, True)
     evt_edm4hep.Control = True
     output = dd.outputFile
     if not dd.outputFile.endswith(dd.outputConfig.myExtension):
          output = dd.outputFile + dd.outputConfig.myExtension
     evt_edm4hep.Output = output
     evt_edm4hep.enableUI()
     Kernel().eventAction().add(evt_edm4hep)
     eventPars = dd.meta.parseEventParameters()
     evt_edm4hep.RunHeader = dd.meta.addParametersToRunHeader(dd)
     evt_edm4hep.EventParametersString, evt_edm4hep.EventParametersInt, evt_edm4hep.EventParametersFloat = eventPars
     evt_edm4hep.RunNumberOffset = dd.meta.runNumberOffset if dd.meta.runNumberOffset > 0 else 0
     evt_edm4hep.EventNumberOffset = dd.meta.eventNumberOffset if dd.meta.eventNumberOffset > 0 else 0
     return None
SIM.outputConfig.userOutputPlugin = setupEDM4hepOutputDR
SIM.outputConfig.myExtension = '.root'

SIM.crossingAngleBoost = 0.0
SIM.vertexOffset = [0.0, 0.0, 0.0, 0.0]
SIM.vertexSigma = [0.0, 0.0, 0.0, 0.0]

SIM.filter.filters = {
     'geantino':{'name':'GeantinoRejectFilter/GeantinoRejector',
                 'parameter':{}},
     'edep1keV':{'name':'EnergyDepositMinimumCut/1keV',
                 'parameter':{'Cut':1.0*keV}},
     'edep0'   :{'name':'EnergyDepositMinimumCut/Cut0',
                 'parameter':{'Cut':0.0}}
}

SIM.action.calo = "SCEPCalSDAction_DRHit"
SIM.action.calorimeterSDTypes = ['SegmentedCrystalCalorimeter']

SIM.filter.calo = "edep0"
print(f'Using filter {SIM.filter.calo} !')

SIM.action.trackerSDTypes = ['tracker']
SIM.action.tracker = (
     'Geant4TrackerWeightedAction',
     {
          'HitPositionCombination':2,
          'CollectSingleDeposits':False
     }
)
SIM.filter.tracker = "edep0"

SIM.part.enableDetailedHitsAndParticleInfo = False
SIM.part.keepAllParticles = True
SIM.part.minDistToParentVertex = 2.2e-10
SIM.part.minimalKineticEnergy = 0.0
SIM.part.printEndTracking = False
SIM.part.printStartTracking = False
SIM.part.saveProcesses = ['Decay']
SIM.part.userParticleHandler = ''

SIM.physics.decays = False
SIM.physics.list = "FTFP_BERT"
SIM.physics.pdgfile = None
SIM.physics.rangecut = None
SIM.physics.rejectPDGs = {1, 2, 3, 4, 5, 6,
                          3201, 3203, 4101, 4103,
                          21, 23, 24, 25, 
                          5401, 2203, 5403,
                          3101, 3103, 4403,
                          2101, 5301, 2103, 5303,
                          4301, 1103, 4303, 5201, 
                          5203, 3303, 4201, 4203, 
                          5101, 5103, 5503}
SIM.physics.zeroTimePDGs = {17, 11, 13, 15}

def setupCerenkovScint(kernel):
     from DDG4 import PhysicsList
     seq = kernel.physicsList()

     scint = PhysicsList(kernel, 'Geant4ScintillationPhysics/ScintillationPhys')
     scint.VerboseLevel = 0
     scint.TrackSecondariesFirst = True
     scint.enableUI()
     seq.adopt(scint)

     cerenkov = PhysicsList(kernel, 'Geant4CerenkovPhysics/CerenkovPhys')
     cerenkov.VerboseLevel = 0
     cerenkov.MaxNumPhotonsPerStep = 10
     cerenkov.MaxBetaChangePerStep = 10.0
     cerenkov.TrackSecondariesFirst = True
     cerenkov.enableUI()
     seq.adopt(cerenkov)

     ph = PhysicsList(kernel, 'Geant4OpticalPhotonPhysics/OpticalGammaPhys')
     ph.addParticleConstructor('G4OpticalPhoton')
     ph.VerboseLevel = 0
     ph.enableUI()
     seq.adopt(ph)

     return None
if opticalPhysics:
     SIM.physics.setupUserPhysics(setupCerenkovScint)
     print ("Optical physics is ON !")
else:
     print ("Optical physics is OFF !")

SIM.random.enableEventSeed = False
SIM.random.file = None
SIM.random.luxury = 1
SIM.random.replace_gRandom = True
SIM.random.seed = None
SIM.random.type = None