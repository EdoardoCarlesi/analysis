##!/bin/bash
# mass_function.sh 
# initialize and start the theoretical mass function fit and calculation routine

if [ $1 -eq 1 ] ; then

grid=$2
mf_lines=$3
pk_lines=$4
mass=$5
box=$6
h=$7
sigma8=$8
bins=$9
halo_file=${10}
pk_file=${11}
fit=${12}
m_max=${13}
m_min=${14}
z=${15}
analysisdir=${17}
mass_function=${18}
output_file=${16}'mass_function.dat'

else

analysisdir=${HOME}/Analysis
mass_function=$analysisdir/bin/amf
# Simulations - file setting
number=60
box=50
grid=512
z=0
model='cde2'

# General settings
bins=50
pk_lines=0
mf_lines=1

# Cosmological settings
mass=0.6975e+8
h=1
sigma8=0.8

# Do a fit? no=0; yes=1
fit=0

# Mass range to calculate the mass function
m_max=1.e+16
m_min=1.e+9

model='ude'
halo_file='/home/edoardo/50/cde2/MERGED/merged_snapshot__60_z0.AHF_halos'
#pk_file='/home/edoardo/Gadget-devel/Pk/pk_'$model'.dat'
#pk_file='/home/edoardo/devel/build/Tables/Om_027/uDE/cdm_pk_z0.dat'
#pk_file='/home/edoardo/devel/build/cdm_z0.dat'
pk_file='/home/edoardo/pk_latta.dat'
output_file=$analysisdir'/output/test_'$model'_mass_function.dat'

fi 

echo $analysisdir

cd $analysisdir/src/
#make clean
make mass_function
mass_function=$analysisdir/bin/mass_function

$mass_function $grid $mf_lines $pk_lines $mass $box $h $sigma8 $bins $halo_file $pk_file $fit $m_max $m_min $z $output_file

exit 0
