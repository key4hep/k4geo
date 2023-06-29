#!/bin/bash
#
# Display detector using the Qt interface of Geant4, via ddsim
# Arguments:
#	1: compact file, with absolute path
ddsim --compactFile $1 --runType qt --macroFile ../example/vis.mac --part.userParticleHandler=''
