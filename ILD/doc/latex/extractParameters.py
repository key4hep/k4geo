import sys

"""
   Extracts global parameters from a DD4hep model
   into a python file with a dictionary.

   The parameters have to be spefified in a simple text file (one per line)

   @author  F.Gaede, CERN/DESY
   @version 1.0

"""
#-----------------------------------------------

#------------------------------------------------
# read the compact xml file from the command line
#
try:
  compactFile = sys.argv[1]
  paramFile   = sys.argv[2]
  dictFile     = sys.argv[3]
except IndexError:
  print " usage:  python extractParameters.py compact.xml param_names.txt pyDict.py"
  print
  sys.exit(1)

#-----------------------------------------------

import os, time, DDG4
from DDG4 import OutputLevel as Output
from SystemOfUnits import *

#-----------------------------------------------

def run():

  kernel = DDG4.Kernel()

  try:
    install_dir = os.environ['DD4hepINSTALL']

  except (KeyError):
    print " please set the environment variable  DD4hepINSTALL  "
    print "        to your DD4hep installation path ! "
    exit(1)

  kernel.loadGeometry("file:"+ compactFile )
  lcdd = kernel.detectorDescription()
  DDG4.importConstants( lcdd ) 

  #--------

  inf = open( paramFile , 'r' )
  outf = open( dictFile , 'w' )

  names = readNames( inf ) 
  
  writeDictionary( names, outf )
  
  inf.close()
  outf.close()


#-----------------------------------------------

def readNames( inf ):
  
  """ 
  read ascii file with parameter names
  """
  names = []
  
  for line in inf:
    cols = line.split()
    for n in cols:
      names.append( n )
      
  return names
#-----------------------------------------------


def writeDictionary( names, of ):
  
  of.write(  '""" \n' )
  of.write(  ' python dictionary with parameters extracted from: ' + compactFile + '\n' )
  of.write(  '""" '+ '\n')
  of.write(  'values={}'+ '\n')
  
  for n in names:
    of.write(  'values["'+n+'"] = ' + str(  getattr( DDG4, n ) ) + '\n')
    
#-----------------------------------------------


def printEnvelopeParameters( det ):
    
    print " ========== ", det , " ================== "
    for p in dict[ det ]:
        print "  ", p ,  getattr( DDG4, p )


#-----------------------------------------------

if __name__ == "__main__":

    run()


    


