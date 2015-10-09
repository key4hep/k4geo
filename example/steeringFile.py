from DDSim.DD4hepSimulation import DD4hepSimulation
from SystemOfUnits import mm, GeV, MeV

SIM = DD4hepSimulation()
SIM.compactFile = "/data/sailer/software/LCGeo_Auth/ILD/compact/ILD_o1_v05/ILD_o1_v05.xml"
SIM.runType = "batch"
SIM.numberOfEvents = 2

SIM.field.eps_min = 1*mm

SIM.gun.multiplicity = 1
SIM.gun.energy = 40*GeV
SIM.gun.direction = (1.0, 1.0, 1.0)

SIM.enableGun = True

SIM.part.minimalKineticEnergy = 1*MeV

SIM.action.mapActions['tpc'] = "TPCSDAction"
