#!/usr//bin/env python
#######################################################
#
# simple script to read guineapig pair background files
# and convert them to LCIO files 
#
# @author F.Gaede, DESY
# @date 18/01/2017
#
# initialize environment:
#  export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib
#
######################################################
import math
import random
from array import array
import sys

# --- LCIO dependencies ---
from pyLCIO import UTIL, EVENT, IMPL, IO, IOIMPL



#================================================
if len( sys.argv ) < 2:
    print " usage: python guineapig_to_lcio.py eepair_guineapig.pair  "
    sys.exit(0)

#=================================================
infile = sys.argv[1]
outfile = infile + ".slcio"


wrt = IOIMPL.LCFactory.getInstance().createLCWriter( )
wrt.open( outfile , EVENT.LCIO.WRITE_NEW ) 
col = IMPL.LCCollectionVec( EVENT.LCIO.MCPARTICLE ) 

evt = IMPL.LCEventImpl() 
evt.setEventNumber( 0 ) 

evt.addCollection( col , "MCParticle" )

# --- 

f = open( infile , 'r' )

for lin in f:

    d = lin.split()

#    print d 

    energy = float( d[0] ) 
    betaX  = float( d[1] ) 
    betaY  = float( d[2] ) 
    betaZ  = float( d[3] ) 
    posX   = float( d[4] ) 
    posY   = float( d[5] ) 
    posZ   = float( d[6] )   

# ----

    pdg = 11
    charge = -1.

    if energy < 0:
        pdg = -11
        energy = - energy 
        charge = 1.

#---

    px = betaX * energy
    py = betaY * energy
    pz = betaZ * energy

    momentum  = array('f',[ px, py, pz ] )  

    nm2mm = 1e-6 

    vtxx = nm2mm * posX
    vtxy = nm2mm * posY
    vtxz = nm2mm * posZ

    vertex = array('d',[vtxx, vtxy, vtxz ] )  



#--------------- create MCParticle -------------------

    mcp = IMPL.MCParticleImpl() 

    mcp.setGeneratorStatus(1) 

    mcp.setMass( 0.0005109989461 )
    
    mcp.setPDG( pdg ) 

    mcp.setMomentum( momentum )
    mcp.setVertex( vertex ) 
    mcp.setCharge( charge ) 


#-------------------------------------------------------


    col.addElement( mcp )


wrt.writeEvent( evt ) 
wrt.close() 

