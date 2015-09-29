#!/usr/bin/python
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
  print " usage:  python ddsim.py compact.xml"
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

numberOfEvents = 10

# ---------------------------------------
# input file: either .slcio or .stdhep
# --------------------------------------

lcioInputFile  = 'mcparticles.slcio'
#lcioInputFile  = 'bbudsc_3evt.stdhep'
#lcioInputFile  = 'bbudsc_3evt_dd4hep.slcio'
#lcioInputFile  = 'mcparticles_single_muon_5GeV_10deg.slcio'
#lcioInputFile  = 'mcparticles_single_muon_5GeV_7deg.slcio'
#lcioInputFile  = 'mcparticles_single_muon_5GeV_85deg.slcio'

#lcioOutputFile  = 'bbudsc_3evt_dd4hep_sim.slcio'
#lcioOutputFile = 'simpleCLIC_single_muon_5GeV_85deg.slcio'
#lcioOutputFile = 'simpleILD_single_muon_5GeV_85deg.slcio'
lcioOutputFile = 'simple_lcio.slcio'
#lcioOutputFile = lcioInputFile[:len(lcioInputFile)-len('.slcio')]+'_SIM.slcio'

physicsList    = 'FTFP_BERT'  # 'QGSP_BERT'

#---subset of detectors to initialize (all if list is empty)
# this does not work, as the complete geometry form the compact.xml files
# is already instantiated ...
#detectorList = []
#detectorList = ['VTX','SIT','FTD']


#-------- dictionary for configuration variables used below ----------------------
cfgDict = {}

#---- B field stepping -------
cfgDict['field.stepper']            = "HelixSimpleRunge"
cfgDict['field.equation']           = "Mag_UsualEqRhs"
cfgDict['field.eps_min']            = 5e-05*mm
cfgDict['field.eps_max']            = 0.001*mm
cfgDict['field.min_chord_step']     = 0.01*mm
cfgDict['field.delta_chord']        = 0.25*mm
cfgDict['field.delta_intersection'] = 1e-05*mm
cfgDict['field.delta_one_step']     = 1e-04*mm


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
#    lcgeo_dir   = os.environ['LCGEO']
  except (KeyError):
    print " please set the environment variable  DD4hepINSTALL  "
    print "        to your DD4hep installation path ! "
    exit(1)

#  example_dir = lcgeo_dir+'/example';

  kernel.loadGeometry("file:"+ compactFile )

  lcdd = kernel.lcdd() 

  DDG4.importConstants( lcdd ) 

#----------------------------------------------------------------------------------

#  simple = DDG4.Geant4( kernel, tracker='Geant4TrackerCombineAction',calo='Geant4ScintillatorCalorimeterAction')
## Apply BirksLaw effect for Scintillator Calorimeter by using 'Geant4ScintillatorCalorimeterAction'.
  simple = DDG4.Geant4( kernel, tracker='Geant4TrackerAction',calo='Geant4ScintillatorCalorimeterAction')

  simple.printDetectors()

  # Configure UI
  #simple.setupCshUI()
  simple.setupUI()

  kernel.UI=""
  kernel.NumEvents=numberOfEvents 

#-----------------------------------------------------------------------------------
#  setup the magnetic field:

  field = simple.addConfig('Geant4FieldTrackingSetupAction/MagFieldTrackingSetup')
  field.stepper            = cfgDict['field.stepper']
  field.equation           = cfgDict['field.equation']          
  field.eps_min            = cfgDict['field.eps_min']           
  field.eps_max            = cfgDict['field.eps_max']           
  field.min_chord_step     = cfgDict['field.min_chord_step']    
  field.delta_chord        = cfgDict['field.delta_chord']       
  field.delta_intersection = cfgDict['field.delta_intersection']
  field.delta_one_step     = cfgDict['field.delta_one_step']    

#----------------------------------------------------------------------------------

  # Configure Run actions
  run1 = DDG4.RunAction(kernel,'Geant4TestRunAction/RunInit')
  kernel.registerGlobalAction(run1)
  kernel.runAction().add(run1)

  # Configure I/O 
  evt_lcio = simple.setupLCIOOutput('LcioOutput', lcioOutputFile )
  
  gen = DDG4.GeneratorAction(kernel,"LCIOInputAction/LCIO1")

  if( lcioInputFile[ (len(lcioInputFile)-6 ) : ] == ".slcio" ):
    gen.Input="LCIOFileReader|"+lcioInputFile
  else:
    gen.Input="LCIOStdHepReader|"+lcioInputFile
    
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
    if( 'tpc' in t.lower() ):
      seq,act = simple.setupTracker( t ,type='TPCSDAction')
    else:
      seq,act = simple.setupTracker( t )
    seq.add(f1)
    act.HitCreationMode = 2

# ---- add the calorimeters:

  for c in cal:
    print 'simple.setupCalorimeter(  ' , c , ')'
    if( 'ecal' in c.lower() ):
      seq,act = simple.setupCalorimeter( c ,type='CaloPreShowerSDAction')
      act.FirstLayerNumber = 1
    else:
      seq,act = simple.setupCalorimeter( c )
    act.HitCreationMode = 2


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
