#!/bin/bash

# Simple test script to demonstrate how to use the 2D U(1) code
# The size of the lattice (L) is hardcoded in the main.cpp file
# to make writing new code simpler. Please please edit and remake
# if you wish to vary L.

#edit this path to be the location of your source code
SRC=/Users/deanhowarth/RPI-github/2D-Schwinger-Model

# configure preamble
#---------------------------------------------------------------
LX=16
LY=32
# Construct the correct executable
cp ${SRC}/wilson/main.cpp main.cpp
cp ${SRC}/wilson/Makefile Makefile

sed -i '.bak' -e s/__LX__/${LX}/g main.cpp Makefile
sed -i '.bak' -e s/__LY__/${LY}/g main.cpp Makefile
sed -i '.bak' -e s/__SRC__/${SRC}/g Makefile

rm *.bak
make

#create directories for data
mkdir -p {gauge,data/{data,creutz,polyakov,rect,top,pion,vacuum}}
#---------------------------------------------------------------

# The value of the coupling in the U(1) 2D theory
BETA=$1

# The total number of HMC iterations to perform.
HMC_ITER=10000
# The number of HMC iterations for thermalisation.
HMC_THERM=250
# The number of HMC iterations to skip bewteen measurements.
HMC_SKIP=1
# Dump the gauge field every HMC_CHKPT iterations after thermalisation.
HMC_CHKPT=5000
# If non-zero, read in the HMC_CHKPT_START gauge field. 
HMC_CHKPT_START=0
# HMC time steps in the integration 
HMC_NSTEP=25
# HMC trajectory time
HMC_TAU=1.0

# Number of APE smearing hits to perform when measuring topology
APE_ITER=1
# The alpha value in the APE smearing
APE_ALPHA=0.5

# The RNG seed
RNG_SEED=1234

# DYNAMIC (1) or QUENCHED (0)
DYN_QUENCH=1

# Lock the Z gauge to unit (1) or allow z dynamics (0)
ZLOCKED=1

# Dynamic fermion parameters
# Fermion mass
MASS=0.1
# Maximum CG iterations
MAX_CG_ITER=1000
# CG tolerance
CG_EPS=1e-16

# Eigensolver parameters
# Tolerance on the residual
TOL=1e-8
# Maximum ARPACK iterations
ARPACK_MAXITER=100000

#polyACC (experimental)
USE_ACC=0
AMAX=11
AMIN=1.0
N_POLY=100

# Measuremets: 1 = measure, 0 = no measure
# Polyakov loops
MEAS_PL=1
# Wilson loops and Creutz ratios
MEAS_WL=1
# Pion Correlation function
MEAS_PC=1
# Vacuum trace
MEAS_VT=0

#construct command 
command="./2D-Wilson-LX$LX-LY$LY $BETA $HMC_ITER $HMC_THERM $HMC_SKIP $HMC_CHKPT 
         $HMC_CHKPT_START $HMC_NSTEP $HMC_TAU $APE_ITER $APE_ALPHA $RNG_SEED 
	 $DYN_QUENCH $MASS $MAX_CG_ITER $CG_EPS $TOL $ARPACK_MAXITER $USE_ACC $AMAX 
    	 $AMIN $N_POLY $MEAS_PL $MEAS_WL $MEAS_PC $MEAS_VT $PULSE"

#echo to stdout for reference
echo $command

$command
