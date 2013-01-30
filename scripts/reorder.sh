#!/bin/bash

N=$1
file='halo_cic_grid_N_'$1'_cells.dat'
echo $file
tmp=f_tmp

sort $file > $tmp
head -n 20 $tmp > sort_halo_$N  
