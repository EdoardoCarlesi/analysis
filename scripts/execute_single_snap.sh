#!/bin/bash
# Script to execute all the analysis settings.

n_procs=4
n_files=16

model='cde000'
zeta='0.046'
masscut='4e+13'

#sys='taurus'
#sys='marenos'
#sys='comodo'
sys='castor'

if [ "$sys" == "taurus" ] ; then
base_analysis=${HOME}/Analysis/
fi

if [ "$sys" == "castor" ] ; then
base_analysis=${HOME}/Analysis/
fi

if [ "$sys" == "marenos" ] ; then
base_analysis=${HOME}/projects/Analysis/
fi

if [ "$sys" == "comodo" ] ; then
base_analysis=${HOME}/Analysis/
fi

cd $base_analysis
make clean; make

$base_analysis'scripts/general_execute.sh' $n_procs $n_files $model $zeta $masscut $sys
