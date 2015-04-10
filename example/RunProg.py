#!/usr/bin/env python

__RCSID__ = "$Id$"
"""

DD4hep simulation with some argument parsing
Based on M. Frank and F. Gaede runSim.py
   @author  A.Sailer
   @version 0.1

"""

import DDG4, DD4hep
from DDG4 import OutputLevel as Output
from SystemOfUnits import *
import argparse

#------------------------------------------------
# read the compact xml file from the command line
#


class DD4hepSimulation(object):
  """Class to hold all the parameters and functions to run simulation"""

  def __init__(self):
    self.compactFile = ""
    self.inputFile = ""
    self.outputFile = "dummyOutput.slcio"
    self.runType = "batch"
    self.printLevel = Output.INFO

    self.numberOfEvents = 0
    self.skipNEvents = 0
    self.physicsList = "FTFP_BERT"
    self.crossingAngleBoost = 0.0
    self.macroFile = ''
    self.gun = False
    self.magneticFieldDict = {}

  @staticmethod
  def getOutputLevel( level ):
    """return output.LEVEL"""
    level = int (level)
    levels = { 1: Output.VERBOSE,
               2: Output.DEBUG,
               3: Output.INFO,
               4: Output.WARNING,
               5: Output.ERROR,
               6: Output.FATAL,
               7: Output.ALWAYS }
    return levels[level]

  def parseOptions(self):
    """parse the command line options"""
    parser = argparse.ArgumentParser("RunProg: Run DDG4")
    parser.add_argument("compactFile", action="store")
    parser.add_argument("runType", action="store", choices=("batch","vis","run","shell"))

    parser.add_argument("--inputFile", "-I", action="store")
    parser.add_argument("--outputFile","-O", action="store", default=self.outputFile)
    parser.add_argument("-v", "--printLevel", action="store", default=3, dest="printLevel")
    parser.add_argument("--nEvents", "-N", action="store", dest="numberOfEvents", default=self.numberOfEvents)
    parser.add_argument("--skipNEvents", action="store", dest="skipNEvents", default=self.skipNEvents)
    parser.add_argument("--physicsList", action="store", dest="physicsList", default=self.physicsList)
    parser.add_argument("--crossingAngleBoost", action="store", dest="crossingAngleBoost", default = self.crossingAngleBoost)
    parser.add_argument("--macroFile", "-M", action="store", dest="macroFile", default = self.macroFile)

    parser.add_argument("--gun", "-G", action="store_true", dest="gun", default=self.gun)

    parsed = parser.parse_args()
    print "parsed things", parsed
    self.compactFile = parsed.compactFile
    self.inputFile = parsed.inputFile
    self.outputFile = parsed.outputFile
    self.runType = parsed.runType
    self.printLevel = self.getOutputLevel(parsed.printLevel)

    self.numberOfEvents = parsed.numberOfEvents
    self.skipNEvents = parsed.skipNEvents
    self.physicsList = parsed.physicsList
    self.crossingAngleBoost = parsed.crossingAngleBoost
    self.macroFile = parsed.macroFile
    self.gun = parsed.gun

    if self.runType == "batch":
      if not self.numberOfEvents:
        print "ERROR: Batch mode requested, but did not set number of events"
        exit(1)
      if not self.inputFile and not self.gun:
        print "ERROR: Batch mode requested, but did not set inputfile or gun"
        exit(1)

  def setupMagneticField(self, simple):
    self.magneticFieldSetup()
    field = simple.addConfig('Geant4FieldTrackingSetupAction/MagFieldTrackingSetup')
    field.stepper            = self.magneticFieldDict['field.stepper']
    field.equation           = self.magneticFieldDict['field.equation']
    field.eps_min            = self.magneticFieldDict['field.eps_min']
    field.eps_max            = self.magneticFieldDict['field.eps_max']
    field.min_chord_step     = self.magneticFieldDict['field.min_chord_step']
    field.delta_chord        = self.magneticFieldDict['field.delta_chord']
    field.delta_intersection = self.magneticFieldDict['field.delta_intersection']
    field.delta_one_step     = self.magneticFieldDict['field.delta_one_step']


  def magneticFieldSetup(self):
    """options for the magneticFieldStepper"""

    #---- B field stepping -------
    self.magneticFieldDict['field.stepper']            = "HelixSimpleRunge"
    self.magneticFieldDict['field.equation']           = "Mag_UsualEqRhs"
    self.magneticFieldDict['field.eps_min']            = 5e-05*mm
    self.magneticFieldDict['field.eps_max']            = 0.001*mm
    self.magneticFieldDict['field.min_chord_step']     = 0.01*mm
    self.magneticFieldDict['field.delta_chord']        = 0.25*mm
    self.magneticFieldDict['field.delta_intersection'] = 1e-05*mm
    self.magneticFieldDict['field.delta_one_step']     = 1e-04*mm

  @staticmethod
  def getDetectorLists( lcdd ):
    ''' get lists of trackers and calorimeters that are defined in lcdd (the compact xml file)'''
  #  if len(detectorList):
  #    print " subset list of detectors given - will only instantiate these: " , detectorList
    trackers,calos = [],[]
    for i in lcdd.detectors():
      det = DDG4.DetElement(i.second)
      name = det.name()
      sd =  lcdd.sensitiveDetector( name )
      if sd.isValid():
        detType = sd.type()
  #      if len(detectorList) and not(name in detectorList):
  #        continue
        print 'getDetectorLists - found active detctor ' ,  name , ' type: ' , type
        if detType == "tracker":
          trackers.append( det.name() )
        if detType == "calorimeter":
          calos.append( det.name() )

    return trackers,calos

#==================================================================================

  def run(self):
    """setup the geometry and dd4hep and geant4 and do what was asked to be done"""

    kernel = DDG4.Kernel()
    DD4hep.setPrintLevel(self.printLevel)
    #kernel.setOutputLevel('Compact',1)

    kernel.loadGeometry("file:"+ self.compactFile )
    lcdd = kernel.lcdd()

    DDG4.importConstants( lcdd )

  #----------------------------------------------------------------------------------

    simple = DDG4.Geant4( kernel, tracker='Geant4TrackerAction',calo='Geant4CalorimeterAction')

    simple.printDetectors()

    if self.runType =="vis":
      simple.setupUI(typ="csh", vis=True, macro=self.macroFile)
    elif self.runType =="run":
      simple.setupUI(typ="csh", vis=False, macro=self.macroFile, ui=False)
    elif self.runType =="shell":
      simple.setupUI(typ="csh", vis=False, macro=None, ui=True)
    elif self.runType == "batch":
      simple.setupUI(typ="csh", vis=False, macro=None, ui=False)
    else:
      print "ERROR: unknown runType"
      exit(1)

    #kernel.UI="csh"
    kernel.NumEvents=self.numberOfEvents

    #-----------------------------------------------------------------------------------
    # setup the magnetic field:
    self.setupMagneticField(simple)

    #----------------------------------------------------------------------------------

    # Configure Run actions
    run1 = DDG4.RunAction(kernel,'Geant4TestRunAction/RunInit')
    kernel.registerGlobalAction(run1)
    kernel.runAction().add(run1)

    # Configure I/O

    evt_lcio = simple.setupLCIOOutput('LcioOutput', self.outputFile)

    actionList = []

    if self.gun and not self.inputFile:
      gun = DDG4.GeneratorAction(kernel,"Geant4ParticleGun/"+"Gun")
      gun.energy      = 10*GeV
      gun.particle    = "e-"
      gun.multiplicity = 1
      gun.position     = (0.0,0.0,0.0)
      gun.isotrop      = True
      gun.direction    = (0,0,1)
      #actionList.append(gun)
      kernel.generatorAction().add(gun)

    if self.inputFile:
      gen = DDG4.GeneratorAction(kernel,"LCIOInputAction/LCIO1")
      gen.Sync = self.skipNEvents
      gen.Input="LCIOFileReader|"+self.inputFile
      actionList.append(gen)

    if self.crossingAngleBoost:
      lbo = DDG4.GeneratorAction(kernel, "Geant4InteractionVertexBoost")
      lbo.Angle = 0.01
      actionList.append(lbo)

    if actionList:
      simple.buildInputStage( actionList , output_level=DDG4.OutputLevel.DEBUG )

    #================================================================================================

    # And handle the simulation particles.
    part = DDG4.GeneratorAction(kernel,"Geant4ParticleHandler/ParticleHandler")
    kernel.generatorAction().adopt(part)
    #part.SaveProcesses = ['conv','Decay']
    part.SaveProcesses = ['Decay']
    part.MinimalKineticEnergy = 1*MeV
    #part.OutputLevel = Output.INFO #generator_output_level
    part.enableUI()
    user = DDG4.Action(kernel,"Geant4TCUserParticleHandler/UserParticleHandler")
    try:
      user.TrackingVolume_Zmax = DDG4.tracker_region_zmax
      user.TrackingVolume_Rmax = DDG4.tracker_region_rmax
    except AttributeError as e:
      print "No Attribute: ", str(e)

    #  user.enableUI()
    part.adopt(user)

    #=================================================================================


    # Setup global filters for use in sensitive detectors

    f1 = DDG4.Filter(kernel,'GeantinoRejectFilter/GeantinoRejector')
    kernel.registerGlobalFilter(f1)

    f4 = DDG4.Filter(kernel,'EnergyDepositMinimumCut')
    f4.Cut = 1.*keV
    kernel.registerGlobalFilter(f4)


    #=================================================================================
    # get lists of trackers and calorimeters in lcdd

    trk,cal = self.getDetectorLists( lcdd )

  # ---- add the trackers:
  # fixme: this assumes the same filters for all trackers ...

    for tracker in trk:
      print 'simple.setupTracker(  ' , tracker , ')'
      seq,act = simple.setupTracker( tracker )
      seq.add(f1)
      act.HitCreationMode = 2

  # ---- add the calorimeters:

    for calo in cal:
      print 'simple.setupCalorimeter(  ' , calo , ')'
      seq,act = simple.setupCalorimeter( calo )


  #=================================================================================
    # Now build the physics list:
    phys = simple.setupPhysics( self.physicsList )

    #fg: do we need these really ?
    #fg:  ph = DDG4.PhysicsList(kernel,'Geant4PhysicsList/Myphysics')
    #fg:  ph.addParticleConstructor('G4BosonConstructor')
    #fg:  ph.addParticleConstructor('G4LeptonConstructor')
    #fg:  ph.addParticleProcess('e[+-]','G4eMultipleScattering',-1,1,1)
    #fg:  ph.addPhysicsConstructor('G4OpticalPhysics')
    #fg:  ph.enableUI()
    #fg:  phys.add(ph)
    #fg:  phys.dump()

    DD4hep.setPrintLevel(self.printLevel)

    kernel.configure()
    kernel.initialize()

    #DDG4.setPrintLevel(Output.DEBUG)
    kernel.run()
    kernel.terminate()

if __name__ == "__main__":
  RUNNER = DD4hepSimulation()
  RUNNER.parseOptions()
  RUNNER.run()
