from DDSim.DD4hepSimulation import DD4hepSimulation
from SystemOfUnits import mm, GeV, MeV, keV

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

## if there is a user provided SDAction which needs additional parameters these can be passed as a dictionary
SIM.action.mapActions['ecal'] = ( "CaloPreShowerSDAction", {"FirstLayerNumber": 1} )

## add filter to sensitive detectors:
SIM.filter.mapDetFilter = {"VXD": "edep1kev"} ## default filter
SIM.filter.mapDetFilter = {"FT": "geantino"} ## default filter
SIM.filter.mapDetFilter = {"FTD": "edep3kev"} ## custom filter
SIM.filter.filters['edep3kev'] = dict(name="EnergyDepositMinimumCut/3keV", parameter={"Cut": 3.0*keV} )
