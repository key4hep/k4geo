"""
Test of ARC detector
The lines 27-95 setup the simulation using DD4hep
The lines 99-127 run the simulation and produce some plots
An exception is thrown in case the simulation failed or the output file is too small.
        Last modification 2023-Jun-23
        Author: Alvaro Tolosa-Delgado
"""

from __future__ import absolute_import, unicode_literals
import logging
import sys
import os

from DDSim.DD4hepSimulation import DD4hepSimulation

import ROOT

if __name__ == "__main__":
    logging.basicConfig(
        format="%(name)-16s %(levelname)s %(message)s",
        level=logging.INFO,
        stream=sys.stdout,
    )
    logger = logging.getLogger("DDSim")

    SIM = DD4hepSimulation()

    # Default is enable visualization of tracks
    SIM.runType = "qt"
    SIM.macroFile = "vis.mac"

    # Ensure that Cerenkov and optical physics are always loaded
    def setupCerenkov(kernel):
        from DDG4 import PhysicsList

        seq = kernel.physicsList()
        cerenkov = PhysicsList(kernel, "Geant4CerenkovPhysics/CerenkovPhys")
        cerenkov.MaxNumPhotonsPerStep = 10
        cerenkov.MaxBetaChangePerStep = 10.0
        cerenkov.TrackSecondariesFirst = False
        cerenkov.VerboseLevel = 0
        cerenkov.enableUI()
        seq.adopt(cerenkov)
        ph = PhysicsList(kernel, "Geant4OpticalPhotonPhysics/OpticalGammaPhys")
        ph.addParticleConstructor("G4OpticalPhoton")
        ph.VerboseLevel = 0
        # BoundaryInvokeSD is disabled:
        # if the SD action is triggered at a boundary, the default cellID calculation fails
        # because DD4hep uses the step midpoint, which lies outside the active volume,
        # making a valid cellID assignment impossible.
        # Temporary workaround: ARC light sensors are given a fake refractive index so
        # photons enter the active volume and the default DD4hep SD action for
        # opticalTracker can be used.
        # A custom SD action would be required to safely enable this option.
        ph.BoundaryInvokeSD = False
        ph.enableUI()
        seq.adopt(ph)
        return None

    SIM.physics.setupUserPhysics(setupCerenkov)

    # Disable filtering for ARC detectors (no hits are filtered out)
    SIM.filter.mapDetFilter["ARCBARREL"] = None
    SIM.filter.mapDetFilter["ARCENDCAP"] = None

    # Use DD4hep optical tracker
    SIM.action.mapActions["ARCBARREL"] = "Geant4OpticalTrackerAction"
    SIM.action.mapActions["ARCENDCAP"] = "Geant4OpticalTrackerAction"

    # Disable user tracker particle handler because no tracking region is defined in arc_full_v0.xml
    SIM.part.userParticleHandler = ""

    # Particle gun settings: pions with fixed energy and theta, varying phi
    SIM.numberOfEvents = 1000
    SIM.enableGun = True
    SIM.gun.energy = "50*GeV"
    SIM.gun.particle = "pi+"
    # SIM.gun.thetaMin = "30*deg"
    # SIM.gun.thetaMax = "150*deg"
    # SIM.gun.phiMin = "0*deg"
    # SIM.gun.phiMax = "240.1*deg"
    SIM.gun.distribution = "uniform"
    SIM.gun.multiplicity = 1
    SIM.gun.position = "0 0 0*cm"

    # Default compact file
    SIM.compactFile = "./compact/arc_full_v0.xml"

    # Output file (assuming CWD)
    SIM.outputFile = "arcsim.root"
    if hasattr(SIM, "outputConfig") and hasattr(SIM.outputConfig, "forceDD4HEP"):
        SIM.outputConfig.forceDD4HEP = True  # use DD4hep root format, not EDM4HEP

    # Override with user options
    SIM.parseOptions()

    # Run the simulation
    try:
        SIM.run()
        if os.path.getsize(SIM.outputFile) < 1000000:
            raise RuntimeError("Output file not found or size less than 1MB")
        # Create some images
        outImagePrefix = "arcsim_"
        ROOT.gROOT.SetBatch(1)
        rootfile = ROOT.TFile(SIM.outputFile)
        number_of_events = []
        EVENT = rootfile.Get("EVENT")
        n1 = EVENT.Draw(
            "ArcCollection.position.Z():ArcCollection.position.phi()",
            "((ArcCollection.cellID>>5)&0x7)==0",
        )
        number_of_events.append(n1)
        ROOT.gPad.SaveAs(outImagePrefix + "barrel_" + SIM.gun.particle + ".png")
        n2 = EVENT.Draw(
            "ArcCollection.position.Y():ArcCollection.position.X()",
            "((ArcCollection.cellID>>5)&0x7)==1",
        )
        number_of_events.append(n2)
        ROOT.gPad.SaveAs(outImagePrefix + "endcapZpos_" + SIM.gun.particle + ".png")
        n3 = EVENT.Draw(
            "ArcCollection.position.Y():ArcCollection.position.X()",
            "((ArcCollection.cellID>>5)&0x7)==2",
        )
        number_of_events.append(n3)
        ROOT.gPad.SaveAs(outImagePrefix + "endcapZneg_" + SIM.gun.particle + ".png")
        EVENT.Draw(
            "ArcCollection.position.Y():ArcCollection.position.X()",
            "((ArcCollection.cellID>>5)&0x7)==2&&ArcCollection.position.Y()>0&&ArcCollection.position.X()>0",
        )
        ROOT.gPad.SaveAs(outImagePrefix + "endcapZneg_" + SIM.gun.particle + "zoom.png")
        rootfile.Close()

        if any(0 == val for val in number_of_events):
            logger.fatal("TEST: failed. At least 1 ARC subsystem did not record any hit")
            raise RuntimeError("TEST: failed. At least 1 ARC subsystem did not record any hit")
        else:
            logger.info("TEST: passed")

    except NameError as e:
        logger.fatal("TEST: failed")
        if "global name" in str(e):
            globalToSet = str(e).split("'")[1]
            logger.fatal(
                "Unknown global variable, please add\nglobal %s\nto your steeringFile"
                % globalToSet
            )
