#!/bin/bash
# merger_tree.sh
# Script to generate satellite catalogues from
# self-merging particles' halos catalogues
# $1 = output file name, 
# $2 = directory containing _particles catalogues
# $3,4 = satellite catalogue file name

echo ''
echo 'Merger Tree, finding parent halos'
echo ''

echo $3
echo $4

# Define files and directories (local and global)
merge=/home/carlesi/amiga-v0.0/bin/MergerTree-file
tempdir=${HOME}/Analysis/temp/

temp_out_file=$tempdir"temp_merger_file"
out_merge=$tempdir"temp_out_merge"

# Go to the particles directories catalogues
cd $2

echo $3 >  $temp_out_file
echo $4 >>  $temp_out_file

echo $1 > $out_merge
echo $1 >> $out_merge

echo $merge 2 $temp_out_file $out_merge
$merge 2 $temp_out_file $out_merge

#rm -rf $temp_out_file $merging_output 
