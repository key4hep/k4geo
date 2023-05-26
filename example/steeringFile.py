from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV, keV

SIM = DD4hepSimulation()
SIM.compactFile = "/data/sailer/software/LCGeo_Auth/ILD/compact/ILD_o1_v05/ILD_o1_v05.xml"
SIM.runType = "batch"
SIM.numberOfEvents = 2

SIM.field.eps_min = .001*mm

SIM.gun.multiplicity = 1
SIM.gun.energy = 40*GeV
SIM.gun.direction = (1.0, 1.0, 1.0)

SIM.enableGun = True

SIM.part.minimalKineticEnergy = 1*MeV

SIM.action.mapActions['tpc'] = "TPCSDAction"

## if there is a user provided SDAction which needs additional parameters these can be passed as a dictionary
SIM.action.mapActions['ecal'] = ( "CaloPreShowerSDAction", {"FirstLayerNumber": 1} )

## add filter to sensitive detectors:
# Either assign dict
SIM.filter.mapDetFilter = {"VXD": "edep1kev"} ## default filter
# or use bracket operator
SIM.filter.mapDetFilter['VXD'] = "edep1kev" ## default filter
SIM.filter.mapDetFilter['FTD'] = ["edep3kev","geantino"] ## custom filter

#create your own filter with a dict containing name and parameter
SIM.filter.filters['edep3kev'] = dict(name="EnergyDepositMinimumCut/3keV", parameter={"Cut": 3.0*keV} )

#chose the physicslist and range cut
SIM.physics.list = "FTFP_BERT"
SIM.physics.rangecut = 1*mm

## set the particle.tbl file to add extra particles to DDsim (B-Baryons)
## use the power of python to get the file from DD4hep wherever it is
import os
if os.path.exists( os.path.join( os.environ.get("DD4hepINSTALL"), "examples/DDG4/examples/particle.tbl") ):
  SIM.physics.pdgfile = os.path.join( os.environ.get("DD4hepINSTALL"), "examples/DDG4/examples/particle.tbl")

## Add parameters to the run header
## This will store the Parameter "MyNewParameter" in the runHeader of the slcio file
## all members are added to the runheader. Isn't Python beautiful?
SIM.MyNewParameter = "Value"
