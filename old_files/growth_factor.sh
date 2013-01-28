#!/bin/bash				

echo "Making the growth factor calculation tool"

if [ $1 -eq 1 ] ; then

echo "*********************"
analysisdir=$2
datadir=$3
k=$4
model=$5
number=$6
info=$7
a_url=$8
output_prefix=$6$7'growth_factor.dat'
use_pk=$9
else 

# Where the Analysis program is found 
analysisdir=${HOME}'/Analysis/'
base_dir=${HOME}'/Gadget-devel/gadget-dedm/'
k=0.01
model=cde3
info=''
number=04-

# Prefix name for the output file
output_prefix=$analysisdir'output/'$number$model$info'_'

# Where the gadget snapshots are stored
datadir=$base_dir'cosmo_'$model'/'

# Where the snapshots' redshifts are written
a_url=$base_dir'/parameterfiles/outputs_lcdm_gas.txt'

fi

cd $analysisdir'src/'
make clean; make growth_factor
echo "Growth factor routine compiled."

echo $analysisdir/bin/growth_factor $k $datadir $a_url $analysisdir $output_prefix
$analysisdir/bin/growth_factor $k $datadir $a_url $analysisdir $output_prefix

echo "Growth factor calculated!"
