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
#					- 11

#letter='l'
#model1='lcdm_features_'$letter
model1='lcdm'
model2='ude'
box_size=75
particle_number=256
n_bins=20
r_bins=7
catalogue_z=0
catalogue_number=28
snap_zero=0
tot_snaps=29
use_snaps=32
k=0.13
fit=0
pk_skip=11
mf_skip=1

#Minimum and maximum mass for the mass function computation
m_min=1.e+11
m_max=1.e+15
m_th=1.e+9

#Particle number threshold
num_th=1.e+2
#Minimum particles per halo
n_min=200

#Computer settings
cpus=6
swap=0

# Radius for the alignement
r_min=1
r_max=50

# Cosmo parameters
h=0.7
s8=0.8
om=0.27
ol=0.73
dc=1.686
virial=2.5
spin=0.15

# Integration redshifts
zMax=3
z0=0
z1=1

# Base directories
#base_ahf=${HOME}/amiga_old/
base_ahf=${HOME}/ahf-v1.0/
base_data=${HOME}/data/dedm/
base_analysis=${HOME}/Analysis/
base_out=$base_analysis/output/
base_temp=$base_analysis/temp/

# General file prefixes 
prefix1=$base_out$model1'-'$box_size'-'$particle_number'-'
prefix2='k_'$k'-'
prefix3='z_'$catalogue_z'-'
prefix4='M_'$m_min'-'
prefix5='N_'$n_min'_'
prefix6='snaps_'$use_snaps

if [ $fit -eq 1 ] ; then 
prefix1=$prefix1'fit-'
fi

# Where the output redshift file for the snapshot is
if [ $particle_number -eq 32 ] ; then
#outputs=$base_data/outputs_lcdm_gas.txt
outputs=$base_data/output_z20.txt
fi
if [ $particle_number -eq 64 ] ; then
#outputs=$base_data/outputs_lcdm_gas.txt
outputs=$base_data/output_z30.txt
fi
if [ $particle_number -eq 128 ] ; then
#outputs=$base_data/outputs_lcdm_gas.txt
outputs=$base_data/output_z40.txt
fi
if [ $particle_number -eq 256 ] ; then
outputs=$base_data/$particle_number/$box_size/outputs_lcdm_gas.txt
fi
if [ $particle_number -eq 512 ] ; then
outputs=$base_data/$particle_number/$box_size/outputs_new.txt
fi

# Use pre-calculated pk = 0 ; calculate pk again = 1
use_pk=0
# Fit data = 1, use th. only = 0
fit=0

DATA1=$base_data/$particle_number/$box_size/$model1/
DATA2=$base_data/$particle_number/$box_size/$model2/

snaps_dir1=$DATA1/snaps/
snaps_dir2=$DATA2/snaps/
halo_dir1=$DATA1/MERGED/
halo_dir2=$DATA2/MERGED/
halo_dir_ahf1=$DATA1/ahf/
halo_dir_ahf2=$DATA2/ahf/

cd $halo_dir1
halo_name1=`ls *0$catalogue_number*_halos`
halo_name2=`ls *0$catalogue_number*_halos`
echo $halo_name1

cd $halo_dir1
profile_name1=`ls *0$catalogue_number*_profiles`
profile_name2=`ls *0$catalogue_number*_profiles`
echo $halo_name1

pk_file_base1=$snaps_dir1/Pk*$particle_number
pk_file_base2=$snaps_dir2/Pk*$particle_number
#pk_file1=${HOME}/Spettri_Latta/P_$letter'_matterpower.dat'
#pk_file1=${HOME}/Spettri_Latta/lcdm.dat
pk_file1=${HOME}/Gadget-devel/N-GenIC/pks/cde099.dat
#pk_file1=`ls $pk_file_base1*snap*$catalogue_number`
pk_file2=`ls $pk_file_base2*snap*$catalogue_number`

halo_file1=$halo_dir1/$halo_name1
halo_file2=$halo_dir2/$halo_name2
profile_file1=$halo_dir1/$profile_name1
profile_file2=$halo_dir2/$profile_name2


# File where the output Pks list is printed
pk_list_file=$base_analysis/temp/pk_files.list
halo_list=$base_temp/halo.list
profile_list=$base_temp/profile.list
subhalo_list=$base_temp/subhalo.list

echo $pk_file_base1
pkfilelist=`ls $pk_file_base1* > $pk_list_file` 
halolist=`ls -r $halo_dir1/*halos > $halo_list` 
profilelist=`ls -r $halo_dir1/*profiles > $profile_list` 
subhalolist=`ls -r $halo_dir1/*substructure > $subhalo_list`

cd $base_analysis/src/
make clean

url_variables=$base_analysis' '$outputs' '$pk_file1' '$halo_file1' '$profile_file1' '$pk_file_base1' '$snaps_dir1' '$halo_dir1
set_variables=$box_size' '$particle_number' '$n_bins' '$pk_skip' '$mf_skip' '$fit' '$catalogue_z' '$m_th' '$m_min' '$m_max' '$r_min' '$r_max' '$r_bins' '$n_min
cosmo_variables=$h' '$s8' '$om' '$ol' '$dc' '$spin' '$virial
extra_variables=$k' '$zMax
halo2_variables=$pk_file2' '$halo_file2' '$profile_file2' '$pk_file_base2' '$snaps_dir2' '$halo_dir2
halo_evolution=$halo_list' '$profile_list' '$subhalo_list' '$use_snaps

all_variables=$url_variables' '$set_variables' '$cosmo_variables' '$extra_variables' '

#echo $all_variables

if [ $1 -eq 1 ] ; then
make fit_nfw
$base_analysis/bin/fit_nfw $all_variables$prefix1$prefix3$prefix4
fi

if [ $1 -eq 2 ] ; then
echo $base_analysis/scripts/make_pk.sh $swap $base_ahf $snaps_dir1 $snap_zero $tot_snaps $particle_number
$base_analysis/scripts/make_pk.sh $swap $base_ahf $snaps_dir1 $snap_zero $tot_snaps $particle_number
fi

if [ $1 -eq 3 ] ; then
find $snaps_dir1'/Pk-'$particle_number* -print > $list_file
make growth_factor
$base_analysis/bin/growth_factor $all_variables$prefix1$prefix2 $pk_list_file
fi

if [ $1 -eq 4 ] ; then
make mass_function
$base_analysis/bin/mass_function $all_variables$prefix1$prefix3
fi

if [ $1 -eq 5 ] ; then
make halo_statistics
$base_analysis/bin/halo_statistics $all_variables$prefix1$prefix3$prefix4
fi

if [ $1 -eq 6 ] ; then
make number_density
$base_analysis/bin/number_density $all_variables$prefix1$prefix4
fi

if [ $1 -eq 7 ] ; then
make halo_evolution
$base_analysis/bin/halo_evolution $all_variables$prefix1$prefix6 $pk_list_file $halo_evolution
fi

if [ $1 -eq 8 ] ; then
make halo_comparison
$base_analysis/bin/halo_comparison $all_variables$prefix1$prefix3$prefix4' '$halo2_variables
fi

if [ $1 -eq 9 ] ; then
make subhalo_statistics 
$base_analysis/bin/subhalo_statistics $all_variables$prefix1$prefix3$prefix5
fi

if [ $1 -eq 10 ] ; then
make theoretical_mass_function
$base_analysis/bin/theoretical_mass_function  $all_variables$prefix1
fi


if [ $1 -eq 11 ] ; then
make test
$base_analysis/src/test  $all_variables$prefix1$prefix4
fi

#rm -rf $base_analysis/temp/*
