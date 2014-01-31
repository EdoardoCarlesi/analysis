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
# MPI Settings - if using MPI
use_mpi=$2
n_procs=$3
use_multiple_cat=$4

if [ $use_mpi -eq 1 ] ; then 
use_multiple_cat=1
fi

# Model and simulation ettings
model1='cde099'
model2='nothing'
box_size=250
particle_number=1024
web_size=256
tot_snaps=61

# Catalogue settings when using one halo catalogue only
catalogue_z=0
catalogue_number=61

# Number of bins for general distributions and for the radial alignment
n_bins=15
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

# Mass threshold for printing halo profiles
m_print=1.e+14

# Minimum particles per halo or minimum mass per halo, spin and virial criterion
n_min=20
m_th=9.e+13
virial=1.5
spin=0.15

# Minimum eigenvalue for the velocity shear tensor
l_web=0.1

# use_n_min = 1 means we use particle number instead of mass as threshold criterion 
use_n_min=0

# use_n_haloes = 0
use_n_haloes=0

# use_criterion = 1 use mass/num, 2 use virialization, 3 use spin, 4 use concentration, 5 use all combined
use_criterion=5

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

dir_ahf=ahf
dir_snaps=snaps
dir_pk=pk
dir_halo=catalogues
dir_web=vweb

# General file prefixes 
prefix=$base_out$model1'-'$box_size'-'$particle_number'-'

# Where the output redshift file for the snapshot is
if [ $particle_number -eq 32 ] ; then
outputs=$base_data/output_z20.txt
fi

if [ $particle_number -eq 64 ] ; then
outputs=$base_data/output_z30.txt
fi

if [ $particle_number -eq 128 ] ; then
outputs=$base_data/output_z40.txt
fi

if [ $particle_number -eq 256 ] ; then
outputs=$base_data/$particle_number/$box_size/outputs_lcdm_gas.txt
fi

if [ $particle_number -eq 512 ] ; then
outputs=$base_data/$particle_number/$box_size/outputs_new.txt
fi

if [ $particle_number -eq 1024 ] ; then
outputs=$base_data/$particle_number/$box_size/outputs_new.txt
fi

DATA1=$base_data/$particle_number/$box_size/$model1/
DATA2=$base_data/$particle_number/$box_size/$model2/

snaps_dir1=$DATA1$dir_snaps/
snaps_dir2=$DATA2$dir_snaps/
halo_dir1=$DATA1$dir_halo/
halo_dir2=$DATA2$dir_halo/
halo_dir_ahf1=$DATA1$dir_ahf/
halo_dir_ahf2=$DATA2$dir_ahf/
pk_dir1=$DATA1$dir_pk/
pk_dir2=$DATA2$dir_pk/
web_dir1=$DATA1$dir_web/

zzzz='.z'
cat_zero='00'

if [ $use_multiple_cat -eq 1 ] ; then
zzzz='0000.z'
#zzzz='0000'
fi

if [ $catalogue_number -gt 9 ] ; then
cat_zero='0'
fi

cd $halo_dir1
halo_name1=`ls *$cat_zero$catalogue_number*$zzzz*_halos`
halo_name2=`ls *$cat_zero$catalogue_number*$zzzz*_halos`

cd $halo_dir1
profile_name1=`ls *$cat_zero$catalogue_number*$zzzz*_profiles`
profile_name2=`ls *$cat_zero$catalogue_number*$zzzz*_profiles`

pk_file1=`ls $pk_dir1*snap*$catalogue_number`
pk_file2=`ls $pk_dir2*snap*$catalogue_number`
echo pk_file1=`ls $pk_dir1*snap*$catalogue_number`
echo 'pk file 1'$pk_file1'*'

if [ -z $pk_file1 ] ; then
pk_file1='pk_not_found.dummy'
echo 'pk file 1' $pk_file1
fi

web_dm_file1=`ls $web_dir1*_dm*$web_size*ascii`
web_gas_file1=`ls $web_dir1*_gas*$web_size*ascii`

if [ -z $web_dm_file1 ] ; then
web_dm_file1='web_not_found.dummy'
echo 'web dm file 1' $web_dm_file1
fi

if [ -z $web_gas_file1 ] ; then
web_gas_file1='web_not_found.dummy'
echo 'web gas file 1' $web_gas_file1
fi

halo_file1=$halo_dir1/$halo_name1
halo_file2=$halo_dir2/$halo_name2
profile_file1=$halo_dir1/$profile_name1
profile_file2=$halo_dir2/$profile_name2

# File where the output Pks list is printed
pk_list=$base_temp/pk.list
halo_list=$base_temp/halo.list
profile_list=$base_temp/profile.list
subhalo_list=$base_temp/subhalo.list

ls -r $snaps_dir1/Pk* > $pk_list 
ls -r $halo_dir1/*$zzzz*halos > $halo_list 
ls -r $halo_dir1/*$zzzz*profiles > $profile_list
ls -r $halo_dir1/*$zzzz*substructure > $subhalo_list

cd $base_analysis/src/
make clean

url_var=$outputs' '$halo_file1' '$profile_file1' '$pk_file1' '$web_dm_file1' '$web_gas_file1
set_var1=$box_size' '$particle_number' '$web_size' '$n_bins' '$n_bins_th' '$r_bins' '$pk_skip' '$mf_skip' '$catalogue_number
set_var2=$fit' '$catalogue_z' '$m_th' '$m_min' '$m_max' '$r_min' '$r_max' '$l_web' '$n_min' '$use_n_min' '$use_n_haloes' '$use_criterion' '$m_print
cosmo_var=$h' '$s8' '$om' '$ol' '$dc' '$spin' '$virial' '$k' '$zMax
evolution_var=$halo_list' '$profile_list' '$subhalo_list' '$pk_list' '$tot_snaps
halo2_var=$pk_file2' '$halo_file2' '$profile_file2' '$snaps_dir2' '$halo_dir2

all_variables=$url_var' '$set_var1' '$set_var2' '$cosmo_var' '$prefix' '$evolution_var' '

execute=$base_analysis

# Check if the Makefile has implemented WITH_MPI or not
$base_analysis/scripts/mpi_check.sh $base_analysis $use_mpi

if [ $use_mpi -eq 1 ] ; then
#execute='mpiexec -n '$n_procs' valgrind -v '$base_analysis
#execute='mpiexec -mca opal_set_max_sys_limits 1 -n '$n_procs' '$base_analysis
#execute='mpiexec -mca opal_set_max_sys_limits 1 -n '$n_procs' valgrind -v '$base_analysis
execute='/home/carlesi/bin/mpiexec -n '$n_procs' '$base_analysis
fi

if [ $1 -eq 3 ] ; then
execute=$base_analysis
fi

if [ $1 -eq 1 ] ; then
make fit_nfw
$execute/bin/fit_nfw $all_variables
fi

if [ $1 -eq 2 ] ; then
echo $execute/scripts/make_pk.sh $swap $base_ahf $snaps_dir1 0 $tot_snaps $particle_number
$execute/scripts/make_pk.sh $swap $base_ahf $snaps_dir1 0 $tot_snaps $particle_number
fi

if [ $1 -eq 3 ] ; then
find $snaps_dir1'/Pk-'$particle_number* -print > $list_file
make growth_factor
$execute/bin/growth_factor $all_variables
fi

if [ $1 -eq 4 ] ; then
make mass_function
$execute/bin/mass_function $all_variables
fi

if [ $1 -eq 5 ] ; then
make halo_statistics
$execute/bin/halo_statistics $all_variables
fi

if [ $1 -eq 6 ] ; then
make number_density
$execute/bin/number_density $all_variables
fi

if [ $1 -eq 7 ] ; then
make halo_evolution
$execute/bin/halo_evolution $all_variables
fi

if [ $1 -eq 8 ] ; then
make halo_comparison
$execute/bin/halo_comparison $all_variables
fi

if [ $1 -eq 9 ] ; then
make subhalo_statistics 
$execute/bin/subhalo_statistics $all_variables
fi

if [ $1 -eq 10 ] ; then
make theoretical_mass_function
$execute/bin/theoretical_mass_function  $all_variables
fi

if [ $1 -eq 11 ] ; then
make test
$execute/bin/test  $all_variables
fi
