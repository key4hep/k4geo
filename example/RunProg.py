#!/usr/bin/env python

__RCSID__ = "$Id$"
"""

DD4hep simulation with some argument parsing
Based on M. Frank and F. Gaede runSim.py
   @author  A.Sailer
   @version 0.1

"""

from DD4hepSimulation import DD4hepSimulation

#------------------------------------------------
# read the compact xml file from the command line
#

if __name__ == "__main__":
  RUNNER = DD4hepSimulation()
  RUNNER.parseOptions()
  RUNNER.run()
