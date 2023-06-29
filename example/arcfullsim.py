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
    SIM.macroFile ='vis.mac'

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
        ph.BoundaryInvokeSD = True
        ph.enableUI()
        seq.adopt(ph)
        return None

    SIM.physics.setupUserPhysics(setupCerenkov)

    # Allow energy depositions to 0 energy in trackers (which include optical detectors)
    SIM.filter.tracker = "edep0"

    # Some detectors are only sensitive to optical photons
    SIM.filter.filters["opticalphotons"] = dict(
        name="ParticleSelectFilter/OpticalPhotonSelector",
        parameter={"particle": "opticalphoton"},
        )
    SIM.filter.mapDetFilter["ARCBARREL"] = "opticalphotons"
    SIM.filter.mapDetFilter["ARCENDCAP"] = "opticalphotons"

    # Use the optical tracker for the PFRICH
    SIM.action.mapActions["ARCBARREL"] = "Geant4OpticalTrackerAction"
    SIM.action.mapActions["ARCENDCAP"] = "Geant4OpticalTrackerAction"

    # Disable user tracker particle handler, so hits can be associated to photons
    SIM.part.userParticleHandler = ""

    # Particle gun settings: pions with fixed energy and theta, varying phi
    SIM.numberOfEvents = 1000
    SIM.enableGun = True
    SIM.gun.energy = "50*GeV"
    SIM.gun.particle = "pi+"
    #SIM.gun.thetaMin = "30*deg"
    #SIM.gun.thetaMax = "150*deg"
    #SIM.gun.phiMin = "0*deg"
    #SIM.gun.phiMax = "240.1*deg"
    SIM.gun.distribution = "uniform"
    SIM.gun.multiplicity = 1
    SIM.gun.position = "0 0 0*cm"


    # Default compact file
    SIM.compactFile = "./compact/arc_full_v0.xml"

    # Output file (assuming CWD)
    SIM.outputFile = "arcsim.root"
    #SIM.outputFile = "arcsim_edm4hep.root"

    # Override with user options
    SIM.parseOptions()

    # Run the simulation
    try:
        SIM.run()
        if os.path.getsize( SIM.outputFile ) < 1000000 :
            raise RuntimeError("Output file not found or size less than 1MB")
        # Create some images
        outImagePrefix = "arcsim_"
        ROOT.gROOT.SetBatch(1)
        rootfile = ROOT.TFile(SIM.outputFile)
        if not "edm4hep" in SIM.outputFile :
            EVENT = rootfile.Get("EVENT")
            EVENT.Draw("ARC_HITS.position.Z():ARC_HITS.position.phi()","((ARC_HITS.cellID>>5)&0x7)==0")
            ROOT.gPad.SaveAs( outImagePrefix + "barrel_" + SIM.gun.particle + ".png")
            EVENT.Draw("ARC_HITS.position.Y():ARC_HITS.position.X()","((ARC_HITS.cellID>>5)&0x7)==1")
            ROOT.gPad.SaveAs( outImagePrefix + "endcapZpos_" + SIM.gun.particle + ".png")
            EVENT.Draw("ARC_HITS.position.Y():ARC_HITS.position.X()","((ARC_HITS.cellID>>5)&0x7)==2")
            ROOT.gPad.SaveAs( outImagePrefix + "endcapZneg_" + SIM.gun.particle + ".png")
            EVENT.Draw("ARC_HITS.position.Y():ARC_HITS.position.X()","((ARC_HITS.cellID>>5)&0x7)==2&&ARC_HITS.position.Y()>0&&ARC_HITS.position.X()>0")
            ROOT.gPad.SaveAs( outImagePrefix + "endcapZneg_" + SIM.gun.particle + "zoom.png")
        else :
            EVENT = rootfile.Get("events")
            EVENT.Draw("ARC_HITS.position.z:atan(ARC_HITS.position.y/ARC_HITS.position.x)","((ARC_HITS.cellID>>5)&0x7)==0&& ARC_HITS.position.x>0")
            ROOT.gPad.SaveAs( outImagePrefix + "barrel_" + SIM.gun.particle + ".png")
            EVENT.Draw("ARC_HITS.position.y:ARC_HITS.position.x","((ARC_HITS.cellID>>5)&0x7)==1")
            ROOT.gPad.SaveAs( outImagePrefix + "endcapZpos_" + SIM.gun.particle + ".png")
            EVENT.Draw("ARC_HITS.position.y:ARC_HITS.position.x","((ARC_HITS.cellID>>5)&0x7)==2")
            ROOT.gPad.SaveAs( outImagePrefix + "endcapZneg_" + SIM.gun.particle + ".png")
            EVENT.Draw("ARC_HITS.position.y:ARC_HITS.position.x","((ARC_HITS.cellID>>5)&0x7)==2&& ARC_HITS.position.x>0&& ARC_HITS.position.y>0")
            ROOT.gPad.SaveAs( outImagePrefix + "endcapZneg_" + SIM.gun.particle + "zoom.png")

        rootfile.Close()

        logger.info("TEST: passed")

    except NameError as e:
        logger.fatal("TEST: failed")
        if "global name" in str(e):
            globalToSet = str(e).split("'")[1]
            logger.fatal("Unknown global variable, please add\nglobal %s\nto your steeringFile" % globalToSet)
