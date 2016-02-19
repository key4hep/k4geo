from DDSim.DD4hepSimulation import DD4hepSimulation
from SystemOfUnits import mm, GeV, MeV

SIM = DD4hepSimulation()

SIM.runType = "batch"
SIM.numberOfEvents = 10

SIM.skipNEvents = 0
SIM.outputFile = "gun_pion_10GeV_SIM_00.slcio"

SIM.compactFile = "../compact/MainTestBeamSetup.xml"
SIM.dumpSteeringFile = "dumpSteering00.xml"

SIM.field.eps_min = 1*mm

SIM.part.minimalKineticEnergy = 1*MeV

SIM.physicsList = "QGSP_BERT"
   
SIM.enableDetailedShowerMode=True

SIM.enableGun = True

SIM.gun.energy = 10*GeV
SIM.gun.particle = "pi+"
#SIM.gun.multiplicity
SIM.gun.position = "0,0,-1000"
#SIM.gun.isotrop
SIM.gun.direction = "0,0,1"
