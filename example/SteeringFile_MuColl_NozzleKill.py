"""ddsim steering for Muon Collider detectors using the NozzleKillAction plugin.

NozzleKillAction (plugins/NozzleKillAction.cpp) stops every track inside two cones
with their tips on the beam axis at z = +-KillZ, opening away from the interaction
point with half-angle OpeningAngle, and ending at the end of the nozzle.
Values are Geant4 internal units; use g4units to express them.

At the beginning of the run the plugin reads the nozzle extent and outer radius
from the volumes of the Geant4 region NozzleRegion.  The cones stop where the
nozzle ends, and if OpeningAngle would bring a cone closer than NozzleMargin to
the outer radius, the angle is reduced to fit (with a warning).  Set
NozzleMargin < 0 to skip that check.
"""

from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import cm, deg

SIM = DD4hepSimulation()

SIM.action.step = (
    "NozzleKillAction",
    {
        "KillZ": 200 * cm,  # cone tips at +-KillZ
        "OpeningAngle": 4 * deg,  # cone half-angle from the beam axis
        "NozzleMargin": 5 * cm,  # keep the cones this far inside the nozzle outer radius
        "NozzleRegion": "NozzleRegion",  # region holding the nozzle volumes
    },
)
