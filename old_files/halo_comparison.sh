#!/bin/bash

echo 'Halo Comparison'
echo ''

if [ $1 -eq 1 ] ; then 

analysisdir=$2
model1=$3
model2=$4
num=$4
box=$5
fname=$6
file1=$7
file2=$8
profile1=$9
profile2=$10
cc_list=$11
lines=$12
halos=$13
bins=$14

else

# Where the Analysis program is found 
analysisdir=${HOME}'/Analysis/'
model1=lcdm
model2=cde1
num=99
box=50
fname=merged_snapshot_$num_z0.AHF_

file1=/home/carlesi/data/dedm/$box/$model1/ahf/MERGED/$fname'halos'
file2=/home/carlesi/data/dedm/$box/$model2/ahf/MERGED/$fname'halos'
profile1=/home/carlesi/data/dedm/$box/$model1/ahf/MERGED/$fname'profiles'
profile2=/home/carlesi/data/dedm/$box/$model2/ahf/MERGED/$fname'profiles'

cc_list=/home/carlesi/amiga-v0.0/MergedTreeVarious/lcdm_lcdmVDE_1.cc_idx

lines=1 
halos=1000
bins=12

fi 

rm -rf $analysisdir/bin/halo_comparison
cd $analysisdir/src
make clean; make halo_comparison
cd ..

if [ -e "bin/halo_comparison" ] ; then

$analysisdir/bin/halo_comparison $file1 $file2 $cc_list $profile1 $profile2 $lines $halos $bins $box

else 

echo 'Could not compile the halo_comparison binary.' 

fi

exit 0
