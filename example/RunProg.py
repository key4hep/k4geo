#!/usr/bin/env python
"""

DD4hep simulation with some argument parsing
Based on M. Frank and F. Gaede runSim.py
   @author  A.Sailer
   @version 0.1

"""

__RCSID__ = "$Id$"

from DD4hepSimulation import DD4hepSimulation

#------------------------------------------------

if __name__ == "__main__":
  RUNNER = DD4hepSimulation()
  RUNNER.parseOptions()
  RUNNER.run()
