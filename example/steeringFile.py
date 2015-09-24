from RunProg import DD4hepSimulation
from SystemOfUnits import GeV
##We need to do this because DD4hepSimulation is calling execfile in a function and thus these things like GeV are only added to the local namespace of the function readSteeringFile
global GeV


def myGunOptions( gun ):
  """set the starting properties of the DDG4 particle gun"""
  print "CREATING MY OWN GUN!!!"
  gun.energy      = 3*GeV
  gun.particle    = "mu-"
  gun.multiplicity = 3
  gun.position     = (0.0,0.0,0.0)
  gun.isotrop      = False
  gun.direction    = (0,0,1)
  return gun

SIM = DD4hepSimulation()
SIM.compactFile = "/data/sailer/software/LCGeo_Auth/ILD/compact/ILD_o1_v05/ILD_o1_v05.xml"
SIM.runType = "vis"
SIM.setGunOptions = myGunOptions
