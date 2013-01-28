##!/bin/bash
# halo_statistics.sh 

if [ $1 -eq 1 ] ; then

halo_file=$2
profile_file=$3
analysisdir=$4
box=$5
mf_lines=$6
m_thresh=$7
bins=$8
virial=$9
output=${10}
rmin=${11} 
rmax=${12}

else

analysisdir=${HOME}/Analysis
amf=$analysisdir/bin/amf
# Simulations - file setting
number=99
box=50
grid=512
z=0
model=lcdm

# General settings
bins=50
pk_lines=11
mf_lines=1

# Cosmological settings
mass=1.e+11
h=0.7
sigma8=0.8

# Do a fit? no=0; yes=1
fit=0

# Mass range to calculate the mass function
m_max=1.e+16
m_min=1.e+10

halo_file=${HOME}/data/dedm/$box/$model/MERGED/merged_snapshots__$number'_z'$z'.AHF_halos'
profile_file=${HOME}/data/dedm/$box/$model/MERGED/merged_snapshots__$number'_z'$z'.AHF_profiles'

fi

cd $analysisdir'/src'
make clean
make halo_statistics
halo_statistics=$analysisdir/bin/halo_statistics

$halo_statistics $halo_file $profile_file $analysisdir $box $mf_lines $m_thresh $bins $virial $output $rmin $rmax

exit 0
