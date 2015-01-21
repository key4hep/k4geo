import sys

"""

   DD4hep simulation example setup using the python configuration

   @author  M.Frank
            F.Gaede - modified for lcgeo
   @version 1.0

"""

#------------------------------------------------
# read the compact xml file from the command line
#
try:
  compactFile = sys.argv[1]
except IndexError:
  print " usage:  python runSim.py compact.xml"
  print
  sys.exit(1)

#-----------------------------------------------

import os, time, DDG4
from DDG4 import OutputLevel as Output
from SystemOfUnits import *


#==================================================================================
###################################################################################
# some configuration variables (could go to a seperate 'steering' file)
#

numberOfEvents = 3

lcioInputFile  = 'mcparticles.slcio'

lcioOutputFile = 'simple_lcio.slcio'

physicsList    = 'FTFP_BERT'  # 'QGSP_BERT'

#---subset of detectors to initialize (all if list is empty)
# this does not work, as the complete geometry form the compact.xml files
# is already instantiated ...
#detectorList = []
#detectorList = ['VTX','SIT','FTD']

#
###################################################################################
#==================================================================================




def getDetectorLists( lcdd ):
  ''' get lists of trackers and calorimeters that are defined in lcdd (the compact xml file)'''
#  if len(detectorList):
#    print " subset list of detectors given - will only instantiate these: " , detectorList
  t,c = [],[]
  for i in lcdd.detectors():
    det = DDG4.DetElement(i.second)
    name = det.name()
    sd =  lcdd.sensitiveDetector( name )
    if sd.isValid():
      type = sd.type()
#      if len(detectorList) and not(name in detectorList):
#        continue
      print 'getDetectorLists - found active detctor ' ,  name , ' type: ' , type
      if type == "tracker":
        t.append( det.name() )
      if type == "calorimeter":
        c.append( det.name() )

  return t,c

#==================================================================================

def run():

  kernel = DDG4.Kernel()

  try:
    install_dir = os.environ['DD4hepINSTALL']
    lcgeo_dir   = os.environ['LCGEO']
  except (KeyError):
    print " please set the environment variables DD4hepINSTALL and  "
    print "        LCGEO to your DD4hep/lcgeo installation pathes ! "
    exit(1)

  example_dir = lcgeo_dir+'/example';

#  kernel.loadGeometry("file:"+lcgeo_dir+"/ILD/compact/ILD_o1_v05.xml")
  kernel.loadGeometry("file:"+ compactFile )

  ##fg: fixme: this defines setup parameters for the B-field
  ##    would like to do this from pyhton ...
  kernel.loadXML("file:"+example_dir+"/physics.xml")

  lcdd = kernel.lcdd() 

  DDG4.importConstants( lcdd ) 

#----------------------------------------------------------------------------------

  simple = DDG4.Simple( kernel, tracker='Geant4TrackerAction',calo='Geant4CalorimeterAction')

  simple.printDetectors()

  # Configure UI
  #simple.setupCshUI()
  simple.setupUI()

  kernel.UI=""
  kernel.NumEvents=numberOfEvents 

#-----------------------------------------------------------------------------------

  # Configure Run actions
  run1 = DDG4.RunAction(kernel,'Geant4TestRunAction/RunInit')
  kernel.registerGlobalAction(run1)
  kernel.runAction().add(run1)

  # Configure I/O 
  evt_lcio = simple.setupLCIOOutput('LcioOutput', lcioOutputFile )
  
  gen = DDG4.GeneratorAction(kernel,"LCIOInputAction/LCIO1")
  gen.Input="LCIOFileReader|"+lcioInputFile

  simple.buildInputStage( [gen] , output_level=DDG4.OutputLevel.INFO )

#================================================================================================

  # And handle the simulation particles.
  part = DDG4.GeneratorAction(kernel,"Geant4ParticleHandler/ParticleHandler")
  kernel.generatorAction().adopt(part)
  #part.SaveProcesses = ['conv','Decay']
  part.SaveProcesses = ['Decay']
  part.MinimalKineticEnergy = 100*MeV
  part.OutputLevel = Output.INFO #generator_output_level
  part.enableUI()
  user = DDG4.Action(kernel,"Geant4TCUserParticleHandler/UserParticleHandler")
  user.TrackingVolume_Zmax = DDG4.tracker_region_zmax
  user.TrackingVolume_Rmax = DDG4.tracker_region_rmax
#  user.enableUI()
  part.adopt(user)

#=================================================================================
  

  # Setup global filters fur use in sensintive detectors

  f1 = DDG4.Filter(kernel,'GeantinoRejectFilter/GeantinoRejector')
  kernel.registerGlobalFilter(f1)

  f4 = DDG4.Filter(kernel,'EnergyDepositMinimumCut')
  f4.Cut = 1.*keV
  kernel.registerGlobalFilter(f4)


#=================================================================================
# get lists of trackers and calorimeters in lcdd

  trk,cal = getDetectorLists( lcdd )

# ---- add the trackers:
# fixme: this assumes the same filters for all trackers ...

  for t in trk:
    print 'simple.setupTracker(  ' , t , ')'
    seq,act = simple.setupTracker( t )
    seq.add(f1)

# ---- add the calorimeters:

  for c in cal:
    print 'simple.setupCalorimeter(  ' , c , ')'
    seq,act = simple.setupCalorimeter( c )


#=================================================================================
  # Now build the physics list:
  phys = simple.setupPhysics( physicsList )

  #fg: do we need these really ?
  #fg:  ph = DDG4.PhysicsList(kernel,'Geant4PhysicsList/Myphysics')
  #fg:  ph.addParticleConstructor('G4BosonConstructor')
  #fg:  ph.addParticleConstructor('G4LeptonConstructor')
  #fg:  ph.addParticleProcess('e[+-]','G4eMultipleScattering',-1,1,1)
  #fg:  ph.addPhysicsConstructor('G4OpticalPhysics')
  #fg:  ph.enableUI()
  #fg:  phys.add(ph)
  #fg:  phys.dump()


  kernel.configure()
  kernel.initialize()

  #DDG4.setPrintLevel(Output.DEBUG)
  kernel.run()
  kernel.terminate()

if __name__ == "__main__":
  run()
