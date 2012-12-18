#!/bin/bash
# Fit a NFW profile

echo 'Fit halos NFW profiles'

if [ $1 == 1 ] ; then
analysisdir=$4

else 

analysisdir=${HOME}/Analysis

fi

echo $analysisdir

rm -rf bin/fit_nfw
cd src/
make clean; make fit_nfw
cd ..

if [ -e "bin/fit_nfw" ] ; then

echo 'Built bin/fit_nfw'

if [ $1 == 1 ] ; then

file=$2
profile=$3
box_size=$5
lines_skip=$6
bins=$7
massCut=$8
mp=$9
npart=$10

else 

box_size=50
file=${HOME}/data/dedm/50/lcdm/ahf/MergedSnapshots/merged_snapshot__99_z0.AHF_halos
lines_skip=1
massCut=2.5e+14
bins=15
mp=1.e+11
npart=512
profile=${HOME}/data/dedm/50/lcdm/ahf/MergedSnapshots/merged_snapshot__99_z0.AHF_profiles

fi

echo $analysisdir/bin/fit_nfw $file $profile $analysisdir $box_size $lines_skip $bins $massCut $mp $npart
$analysisdir/bin/fit_nfw $file $profile $analysisdir $box_size $lines_skip $bins $massCut $mp $npart

else 

echo 'Could not compile fit_nfw.'

fi

rm -rf 'temp/h*.tmp'

exit 0
