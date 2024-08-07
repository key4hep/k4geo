# ==========================================================================
#  AIDA Detector description implementation
# --------------------------------------------------------------------------
# Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
# All rights reserved.
#
# For the licensing terms see $DD4hepINSTALL/LICENSE.
# For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
#
# ==========================================================================

from __future__ import absolute_import, unicode_literals
import os
import time
import DDG4
from DDG4 import OutputLevel as Output
from g4units import GeV, MeV, m

def run():
  args = DDG4.CommandLine()
  kernel = DDG4.Kernel()
  logger = DDG4.Logger('DRCfastSD')

  kernel.loadGeometry(str("./IDEA_o1_v03/IDEA_o1_v03.xml"))

  detectorDescription = kernel.detectorDescription()

  DDG4.importConstants(detectorDescription, debug=False)

  geant4 = DDG4.Geant4(kernel)
  geant4.printDetectors()

  # Configure UI
  ui = geant4.setupCshUI()
  ui.Commands = ['/run/beamOn ' + str("10"), '/ddg4/UI/terminate']
  #
  # Configure field
  geant4.setupTrackingField()
  #
  # Configure G4 geometry setup
  seq, act = geant4.addDetectorConstruction('Geant4DetectorGeometryConstruction/ConstructGeo')
  #
  # Assign sensitive detectors according to the declarations 'tracker' or 'calorimeter', etc
  # seq, act = geant4.addDetectorConstruction('Geant4DetectorSensitivesConstruction/ConstructSD')
  # Assign sensitive detectors in Geant4 by matching a regular expression in the detector sub-tree
  seq, act = geant4.addDetectorConstruction('Geant4RegexSensitivesConstruction/ConstructSDRegEx')
  act.Detector = 'DRcalo'
  act.OutputLevel = Output.ALWAYS
  act.Match = ['(.*)(core|clad)(.*)'] # regex is expensive, always try to avoid complex match
                                      # overhead time corresponds to the O(n) of the list size
  #
  # Configure I/O
  evt_root = DDG4.EventAction(kernel, 'Geant4Output2EDM4hep_DRC/fastSDtest.root', True)
  evt_root.Control = True
  output = "fastSDtest.root"
  evt_root.Output = output
  evt_root.enableUI()
  kernel.eventAction().add(evt_root)
  #
  # Setup particle gun
  gun = geant4.setupGun('Gun', particle='e-', energy=20 * GeV, multiplicity=1)
  gun.direction = (0.02, 1, 0.009)
  gun.OutputLevel = Output.INFO
  gun.enableUI()
  #
  # And handle the simulation particles.
  part = DDG4.GeneratorAction(kernel, 'Geant4ParticleHandler/ParticleHandler')
  part.SaveProcesses = ['Decay']
  part.MinimalKineticEnergy = 50 * MeV
  kernel.generatorAction().adopt(part)
  #
  # Map sensitive detectors
  sd = geant4.description.sensitiveDetector(str('DRcalo'))
  logger.info(f'+++ DRcalo: SD type: {str(sd.type())}')
  seq, act = geant4.setupCalorimeter('DRcalo','DRCaloSDAction')
  act.skipScint = True # we skip the optical photon propagation for scint channels
  #
  # Now build the physics list:
  phys = geant4.setupPhysics(str('FTFP_BERT'))
  phys.dump()

  seq = kernel.physicsList()
  cerenkov = DDG4.PhysicsList(kernel, 'Geant4CerenkovPhysics/CerenkovPhys')
  cerenkov.TrackSecondariesFirst = True
  cerenkov.VerboseLevel = 1
  cerenkov.enableUI()
  seq.adopt(cerenkov)

  opt = DDG4.PhysicsList(kernel, 'Geant4OpticalPhotonPhysics/OpticalGammaPhys')
  opt.addParticleConstructor('G4OpticalPhoton')
  opt.VerboseLevel = 1
  opt.BoundaryInvokeSD = True
  opt.enableUI()
  seq.adopt(opt)

  geant4.execute()

  #kernel.configure()
  #kernel.initialize()
  #kernel.run()
  #kernel.terminate()


if __name__ == "__main__":
  run()
