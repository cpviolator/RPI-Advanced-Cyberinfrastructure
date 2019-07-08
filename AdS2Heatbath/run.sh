#!/bin/bash

LEVELS=$1
TRIANGULATION=$2
THERM=10000
MEASUREMENTS=10000000
DELTA_PHI=1.0

MSQR=$3
LAMBDA=$4

SEED=1234

COMMAND="./ads_heatbath --levels ${LEVELS} --q ${TRIANGULATION} --nTherm ${THERM} \
			--nMeas ${MEASUREMENTS} --deltaPhi ${DELTA_PHI} 
			--lambda ${LAMBDA} --mSqr ${MSQR} --seed ${SEED}"

echo $COMMAND

$COMMAND
