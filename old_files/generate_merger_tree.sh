#!/bin/bash
# merger_tree.sh initialize and start the
# halo number density calculation routine

if [ $1 -eq 1 ] ; then

lines=$2
nMin=$3
bins=$4
vir=$5
mtype=$6
mtype2=$7
mass=$8
model=$9
model2=$10
box=$11
base_dir=$12
out=$13
out_bin=$14
halo_dir2=$15
halo_dir1=$16
cc_list=$17

else

analysisdir=${HOME}/Analysis/
merger_tree=$analysisdir/bin/merger_tree
# model type:
# 1=lcdm_1, 3=vde_1, 4=lcdmvde_1, 5=lcdm_05, 6=vde_05, 8=lcdmvde_05

lines=1
nMin=1000
bins=7
vir=10.6
mtype=5
mtype2=6
mass=0.5e+14
model='/lcdm'
model2='/ude'
box=50
base_dir=${HOME}'/data/dedm/'

out='../output/'$model$box'_z_formation'
out_bin='../output/'$model$box'_z_formation_binned'

#halo_dir1='/home/carlesi/data/vde/1000/ahf/MergedSnapshots/'
halo_dir2=$base_dir$box$model2'/MERGED/'
halo_dir1=$base_dir$box$model'/MERGED/'

#cc_list=/home/carlesi/amiga-v0.0/MergedTreeVarious/lcdm_lcdmVDE.ccomp_idx
#cc_list=/home/carlesi/amiga-v0.0/MergedTreeVarious/lcdm_vde_500.ccomp_idx
#cc_list=/home/carlesi/amiga-v0.0/MergedTreeVarious/lcdm_lcdmVDE_1.cc_idx
cc_list=/home/carlesi/amiga-v0.0/MergedTreeVarious/lcdm-vde-1gpc-b.ccomp_idx

fi

rm -rf $merger_tree
cd $analysisdir'/src'
make clean; make merger_tree

if [ -e "../bin/merger_tree" ] ; then
echo 'Built bin/merger_tree'

$merger_tree $lines $nMin $mtype $mtype2 $bins $vir $mass $analysisdir $halo_dir1 $halo_dir2 $cc_list $out $out_bin

else
echo 'Could not compile merger_tree'
fi

exit 0
