from DDSim.DD4hepSimulation import DD4hepSimulation
from SystemOfUnits import mm, GeV, MeV

SIM = DD4hepSimulation()

SIM.compactFile = "../../DD4HEP/compact/TBModel2015.xml"
SIM.runType = "batch"
SIM.macroFile = "vis.mac"
#SIM.inputFiles = "mcparticles.slcio"
SIM.outputFile = "DD4hep_mu-_8GeV_QGSP_BERT_10k.slcio"

SIM.numberOfEvents = 10000
SIM.skipNEvents = 0
SIM.physicsList = "QGSP_BERT"
SIM.dumpSteeringFile = "TBModel2015_dump.xml"
SIM.enableDetailedShowerMode=True

SIM.random.seed = "0123456789"
SIM.field.eps_min = 1*mm
SIM.part.minimalKineticEnergy = 1*MeV

SIM.action.calo = "Geant4ScintillatorCalorimeterAction"

## set the particle.tbl file to add extra particles to DDsim (B-Baryons)
## use the power of python to get the file from DD4hep wherever it is
import os
if os.path.exists( os.path.join( os.environ.get("DD4hepINSTALL"), "examples/DDG4/examples/particle.tbl") ):
  SIM.physics.pdgfile = os.path.join( os.environ.get("DD4hepINSTALL"), "examples/DDG4/examples/particle.tbl")


SIM.enableGun = True
SIM.gun.particle = "mu-"
SIM.gun.energy = 8*GeV
SIM.gun.position = "0, 0, -1000"
SIM.gun.direction = "0,0,1"

#SIM.gun.isotrop
SIM.gun.multiplicity = 1
