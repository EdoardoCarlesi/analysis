#!/bin/bash
# self_merger_tree.sh
# Script to generate satellite catalogues from
# self-merging particles' halos catalogues
# $1 = output file name, $2 = directory containing _particles catalogues
# $3 = satellite catalogue file name

echo ''
echo 'Self Merger Tree, finding subhalo catalogue'
echo ''

echo $3

# Define files and directories (local and global)
merge=/home/carlesi/amiga-v0.0/bin/MergerTree-file
tempdir=${HOME}/Analysis/temp/

particles_file=$tempdir$1
merging_output=$particles_file"-sat_catalogue"
temp_out_file=$tempdir"temp_satellite_file"

# Go to the particles directories catalogues
cd $2

# Generate i/o files list
ls *_particles > $particles_file

#echo  $particles_file >  $temp_out_file
#echo  $particles_file >>  $temp_out_file
echo  $3 >  $temp_out_file
echo  $3 >>  $temp_out_file

sed s/_particles/_satellites/g <$temp_out_file >$merging_output

echo $merge 2 $temp_out_file $merging_output
$merge 2 $temp_out_file $merging_output

rm -rf $temp_out_file $merging_output $particles_file
