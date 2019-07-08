#!/bin/bash

# Simple test script to demonstrate how to use the 2D U(1) code
# The size of the lattice (L) is hardcoded in the main.cpp file
# to make writing new code simpler. Please please edit and remake
# if you wish to vary L.

#create directories for data
mkdir -p {gauge,data/{data,creutz,polyakov,rect,top,pion,vacuum}}
#---------------------------------------------------------------

# The value of the coupling in the U(1) 2D theory
BETA=4.0

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

# DYNAMIC (true) or QUENCHED (false)
DYN_QUENCH='true'

# Lock the Z gauge to unit (1) or allow z dynamics (0)
ZLOCKED=1

# Dynamic fermion parameters
# Fermion mass
MASS=0.1
# Maximum CG iterations
MAX_CG_ITER=1000
# CG tolerance
CG_EPS=1e-16

# Polyakov loops
MEAS_PL=true
# Wilson loops and Creutz ratios
MEAS_WL=true
# Pion Correlation function
MEAS_PC=true
# Vacuum trace
MEAS_VT=false

#construct command 
command="./2D-Wilson --beta ${BETA} \
		     --iterHMC ${HMC_ITER} \
		     --therm ${HMC_THERM} \
		     --skip ${HMC_SKIP} \
		     --checkpoint ${HMC_CHKPT} \
		     --checkpointStart ${HMC_CHKPT_START} \
		     --nSteps ${HMC_NSTEP} \
		     --tau ${HMC_TAU} \
                     --smearIter ${APE_ITER} \
		     --smearAlpha ${APE_ALPHA} \
	 	     --dynamic ${DYN_QUENCH} \
		     --mass ${MASS} \
		     --tol ${CG_EPS} \
		     --maxIterCG ${MAX_CG_ITER} \
		     --measPL ${MEAS_PL} \
		     --measWL ${MEAS_WL} \
		     --measPC ${MEAS_PC} \
		     --measVT ${MEAS_VT}"

#echo to stdout for reference
echo $command

$command
