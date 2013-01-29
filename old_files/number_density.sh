##!/bin/bash
# number_density.sh initialize and start the
# analytical and numerical halo number density calculation routine

analysisdir=${HOME}/Analysis
a_number_density=$analysisdir/bin/analytical_number_density
cd $analysisdir

rm -rf $a_number_density
cd $analysisdir'/src'
make clean; make analytical_number_density

# General settings
bins=40
num=512
pk_lines=11
mf_lines=1
z_0=0
z_1=4
M=0
invert=0
fit=1

box=50
typ=1
halo_dir=/home/carlesi/data/lcdm/500/ahf/MergedSnapshots/
mass=6.95e+10
h=0.7
oLam=0.73
oMat=0.27
s8=0.8
echo 'z1: ' 
read z_1

echo ''
echo $a_number_density $bins $num $pk_lines $mf_lines $typ $halo_dir $box $mass $z_1 $h $M $z_0 $oLam $oMat $s8
$a_number_density $bins $num $pk_lines $mf_lines $typ $halo_dir $box $mass $z_1 $h $M $z_0 $oLam $oMat $s8 $invert $fit
echo ''

exit 0

#h=0.62
#oLam=0.612
#oMat=0.388
#s8=0.83

# VDE 500/1000 Mpc settings
#box=1000
#typ=3
#mass=1.e+11
#halo_dir=/home/carlesi/data/vde/500/ahf/MergedSnapshots/
#halo_dir=/home/carlesi/data/vde/1000/ahf/MergedSnapshots/
#mass=8.02e+11
#h=0.62
#oLam=0.612
#oMat=0.388
#s8=0.83

#echo 'z0: '
#read z_0




