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

# Model and simulation ettings
model1='lcdm'
model2='ude'
box_size=75
particle_number=256
tot_snaps=29

# Catalogue settings when using one halo catalogue only
catalogue_z=0
catalogue_number=028

# Number of bins for general distributions and for the radial alignment
n_bins=20
r_bins=7

# Scale of P(k) for growth factor calculation
k=0.13

# Lines to skip when reading P(k) and mass function files
pk_skip=11
mf_skip=1

#Minimum and maximum mass for the mass function computation
m_min=1.e+11
m_max=1.e+15
m_th=1.e+10

#Minimum particles per halo
n_min=200

# use_n_min = 1 means we use particle number instead of mass as threshold criterion 
use_n_min=1

# use_n_haloes = 0 means use all haloes, otherwise as specified
use_n_haloes=0

#Computer settings
swap=0

# Radius for the halo alignement
r_min=1
r_max=50

# Cosmological parameters
h=0.7
s8=0.8
om=0.27
ol=0.73
dc=1.686
virial=2.5
spin=0.15

# Integration redshifts
z0=0
z1=1
zMax=3

# Fit data = 1, use th. only = 0
fit=0

# Base directories
base_ahf=${HOME}/ahf-v1.0/
base_data=${HOME}/data/dedm/
base_analysis=${HOME}/Analysis/
base_out=$base_analysis/output/
base_temp=$base_analysis/temp/

dir_ahf=ahf
dir_snaps=snaps
dir_halo=catalogues01

# General file prefixes 
prefix1=$base_out$model1'-'$box_size'-'$particle_number'-'
prefix2='k_'$k'-'
prefix3='z_'$catalogue_z'-'
prefix4='M_'$m_min'-'
prefix5='N_'$n_min'_'
prefix6='snaps_'$tot_snaps

if [ $fit -eq 1 ] ; then 
prefix1=$prefix1'fit-'
fi

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

zzzz='.z'

if [ $use_mpi -eq 1 ] ; then
zzzz='0000.z'
fi

cd $halo_dir1
halo_name1=`ls *$catalogue_number*$zzzz*_halos`
halo_name2=`ls *$catalogue_number*$zzzz*_halos`

cd $halo_dir1
profile_name1=`ls *$catalogue_number*$zzzz*_profiles`
profile_name2=`ls *$catalogue_number*$zzzz*_profiles`

pk_file_base1=$snaps_dir1/Pk*$particle_number
pk_file_base2=$snaps_dir2/Pk*$particle_number
pk_file1=`ls $pk_file_base1*snap*$catalogue_number`
pk_file2=`ls $pk_file_base2*snap*$catalogue_number`

halo_file1=$halo_dir1/$halo_name1
halo_file2=$halo_dir2/$halo_name2
profile_file1=$halo_dir1/$profile_name1
profile_file2=$halo_dir2/$profile_name2

# File where the output Pks list is printed
pk_list=$base_temp/pk.list
halo_list=$base_temp/halo.list
profile_list=$base_temp/profile.list
subhalo_list=$base_temp/subhalo.list

ls $pk_file_base1/* > $pk_list 
ls -r $halo_dir1/*$zzzz*halos > $halo_list 
ls -r $halo_dir1/*$zzzz*profiles > $profile_list
ls -r $halo_dir1/*$zzzz*substructure > $subhalo_list

cd $base_analysis/src/
make clean

url_variables=$base_analysis' '$outputs' '$pk_file1' '$halo_file1' '$profile_file1' '$pk_file_base1' '$snaps_dir1' '$halo_dir1
set_variables=$box_size' '$particle_number' '$n_bins' '$pk_skip' '$mf_skip' '$fit' '$catalogue_z' '$m_th' '$m_min' '$m_max' '$r_min' '$r_max' '$r_bins' '$n_min' '$use_n_min' '$use_n_haloes
cosmo_variables=$h' '$s8' '$om' '$ol' '$dc' '$spin' '$virial
extra_variables=$k' '$zMax
halo_evolution=$halo_list' '$profile_list' '$subhalo_list' '$pk_list' '$tot_snaps
halo2_variables=$pk_file2' '$halo_file2' '$profile_file2' '$pk_file_base2' '$snaps_dir2' '$halo_dir2

all_variables=$url_variables' '$set_variables' '$cosmo_variables' '$extra_variables' '

execute=$base_analysis

if [ $use_mpi -eq 1 ] ; then
#execute='mpiexec -n '$n_procs' valgrind -v '$base_analysis
execute='mpiexec -n '$n_procs' '$base_analysis
fi

if [ $1 -eq 1 ] ; then
make fit_nfw
$execute/bin/fit_nfw $all_variables$prefix1$prefix3$prefix4
fi

if [ $1 -eq 2 ] ; then
echo $execute/scripts/make_pk.sh $swap $base_ahf $snaps_dir1 0 $tot_snaps $particle_number
$execute/scripts/make_pk.sh $swap $base_ahf $snaps_dir1 0 $tot_snaps $particle_number
fi

if [ $1 -eq 3 ] ; then
find $snaps_dir1'/Pk-'$particle_number* -print > $list_file
make growth_factor
$execute/bin/growth_factor $all_variables$prefix1$prefix2 $pk_list_file
fi

if [ $1 -eq 4 ] ; then
make mass_function
$execute/bin/mass_function $all_variables$prefix1$prefix3
fi

if [ $1 -eq 5 ] ; then
make halo_statistics
$execute/bin/halo_statistics $all_variables$prefix1$prefix3$prefix4 $halo_evolution
fi

if [ $1 -eq 6 ] ; then
make number_density
$execute/bin/number_density $all_variables$prefix1$prefix4
fi

if [ $1 -eq 7 ] ; then
make halo_evolution
$execute/bin/halo_evolution $all_variables$prefix1$prefix6 $pk_list_file $halo_evolution
fi

if [ $1 -eq 8 ] ; then
make halo_comparison
$execute/bin/halo_comparison $all_variables$prefix1$prefix3$prefix4' '$halo2_variables
fi

if [ $1 -eq 9 ] ; then
make subhalo_statistics 
$execute/bin/subhalo_statistics $all_variables$prefix1$prefix3$prefix5
fi

if [ $1 -eq 10 ] ; then
make theoretical_mass_function
$execute/bin/theoretical_mass_function  $all_variables$prefix1
fi

if [ $1 -eq 11 ] ; then
make test
$execute/bin/test  $all_variables$prefix1$prefix4
fi
