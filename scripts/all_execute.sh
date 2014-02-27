#!/bin/bash
# Script to execute all the analysis settings.
# $1 (input argument) must be set to: 	- 1 (fit_nfw)
#					- 2 (make_pk.sh)	
#					- 3 (growth_factor)
#					- 4 (mass_function)
#					- 5 (halo_statistics)
#					- 6
#					- 7
#					- 8
#					- 9
#					- 10 (theoretical_mass_function)
#					- 11 (test)
#
#sys='taurus'
sys='comodo'

# Input parameters
n_procs=$1
n_files=$2

# MPI Settings - if using MPI
use_mpi='1'

if [ $use_mpi -eq 1 ] ; then 
use_multiple_cat=1
else
use_multiple_cat=0
fi

# Model
model='lcdm'

# Analyze catalog at a given redshift
zzzz='.z'

if [ $use_multiple_cat -eq 1 ] ; then
zzzz='z0.503'
fi

# Simulation ettings
box_size=250
particle_number=1024
web_size=256
tot_snaps=61

# Catalogue settings when using one halo catalogue only
catalogue_z=0
catalogue_number=61

# Number of bins for general distributions and for the radial alignment
n_bins=50
n_bins_th=200
r_bins=11

# Scale of P(k) for growth factor calculation
k=0.13

# Lines to skip when reading P(k) and mass function files
pk_skip=1
mf_skip=1

# Minimum and maximum mass for the mass function computation
m_min=1.e+9
m_max=1.e+15

# Mass threshold for printing individual halo profiles
m_print=1.e+14

# Minimum particles per halo or minimum mass per halo, spin and virial criterion
n_min=20
m_th=9.e+11
m_print=1.e+10
virial=1.5
spin=0.15

# Minimum eigenvalue for the velocity shear tensor
l_web=0.1

# use_n_min = 1 means we use particle number instead of mass as threshold criterion 
use_n_min=0

# use_n_haloes = 0
use_n_haloes=0

# use_criterion = 1 use mass/num, 2 use virialization, 3 use spin, 4 use concentration, 5 use all combined
use_criterion=1

#Computer settings
swap=0

# Radius for the halo alignement
r_min=3
r_max=100

# Cosmological parameters
h=0.7
s8=0.8
om=0.27
ol=0.73
dc=1.686

# Integration redshifts
z0=0
z1=1
zMax=3

# Fit data = 1, use th. only = 0
fit=0

# Base directories
base_ahf=${HOME}/AHF/latest_ahf/
base_data=${HOME}/data/dedm/
base_analysis=${HOME}/Analysis/
base_out=$base_analysis/output/
base_temp=$base_analysis/temp/

DATA=$base_data/$particle_number/$box_size/

ahf_dir=$DATA'ahf/'$model'/'
snaps_dir=$DATA'snaps/'$model'/'
pk_dir=$DATA'pk/'$model'/'
halo_dir=$DATA'catalogues/'$model'/'
web_dir=$DATA'vweb/'$model'/'

# General file prefixes 
prefix=$base_out$model'-'$box_size'-'$particle_number'-'$zzzz

# Expansion factor outputs
outputs=$base_data/$particle_number/$box_size/outputs_new.txt

cat_zero='00'

if [ $catalogue_number -gt 9 ] ; then
cat_zero='0'
fi

halo_name=`ls $halo_dir*$zzzz*_halos`
profile_name=`ls $halo_dir*$zzzz*_profiles`

if [ $use_mpi -eq 1 ] ; then 
halo_name=`ls $halo_dir*0000*$zzzz*_halos`
profile_name=`ls $halo_dir*0000*$zzzz*_profiles`
fi

pk_file=`ls $pk_dir*snap*$catalogue_number`

echo pk_file=`ls $pk_dir*$catalogue_number`

if [ -z $pk_file ] ; then
pk_file='pk_not_found.dummy'
echo 'pk file ' $pk_file
fi

web_dm_file=`ls $web_dir*_dm*$web_size*ascii`
web_gas_file=`ls $web_dir*_gas*$web_size*ascii`

if [ -z $web_dm_file ] ; then
web_dm_file='web_not_found.dummy'
echo 'web dm file' $web_dm_file
fi

if [ -z $web_gas_file ] ; then
web_gas_file='web_not_found.dummy'
echo 'web gas file' $web_gas_file
fi

halo_file=$halo_dir/$halo_name
profile_file=$halo_dir/$profile_name

# File where the output lists are printed
pk_list=$base_temp/pk.list
halo_list=$base_temp/halo.list
profile_list=$base_temp/profile.list
subhalo_list=$base_temp/subhalo.list

ls -r $snaps_dir/Pk* > $pk_list 
ls -r $halo_dir/*$zzzz*halos > $halo_list 
ls -r $halo_dir/*$zzzz*profiles > $profile_list
ls -r $halo_dir/*$zzzz*substructure > $subhalo_list

cd $base_analysis/src/
make clean

url_var=$n_files' '$outputs' '$halo_file' '$profile_file' '$pk_file' '$web_dm_file' '$web_gas_file
set_var1=$box_size' '$particle_number' '$web_size' '$n_bins' '$n_bins_th' '$r_bins' '$pk_skip' '$mf_skip' '$catalogue_number
set_var2=$fit' '$catalogue_z' '$m_th' '$m_min' '$m_max' '$r_min' '$r_max' '$l_web' '$n_min' '$use_n_min' '$use_n_haloes' '$use_criterion' '$m_print
cosmo_var=$h' '$s8' '$om' '$ol' '$dc' '$spin' '$virial' '$k' '$zMax
evolution_var=$halo_list' '$profile_list' '$subhalo_list' '$pk_list' '$tot_snaps' '$zzzz

all_variables=$url_var' '$set_var1' '$set_var2' '$cosmo_var' '$prefix' '$evolution_var

execute=$base_analysis

# Check if the Makefile has implemented WITH_MPI or not
$base_analysis/scripts/mpi_check.sh $base_analysis $use_mpi

if [ $use_mpi -eq 1 ] ; then

#execute='mpiexec -n '$n_procs' valgrind -v '$base_analysis
#execute='mpiexec -mca opal_set_max_sys_limits 1 -n '$n_procs' '$base_analysis
#execute='mpiexec -mca opal_set_max_sys_limits 1 -n '$n_procs' valgrind -v '$base_analysis

if [ "$sys" == "taurus" ] ; then
execute='mpiexec.openmpi -n '$n_procs' '$base_analysis
fi

if [ "$sys" == "comodo" ] ; then
execute='mpiexec -n '$n_procs' '$base_analysis
fi

fi 

# Finally execute the program

make halo_statistics
$execute/bin/halo_statistics $all_variables
