#!/bin/bash

home=${HOME}/Analysis/scripts/
out=${HOME}/Analysis/output/
all_ex=$home/all_execute.sh

N_Mod=3
N_Out=3

N_proc=4
N_file=256

model[0]='cde000'
model[1]='cde099'
model[2]='lcdm'

z[0]='0.000'
z[1]='0.008'
z[2]='0.027'

mass[0]=7.e+13
mass[1]=9.e+13
mass[2]=3.e+14

for (( i=0; i<$N_Out; i++ ))
do

zeta=${z[$i]}
mass=${mass[$i]}

for (( j=2; j<$N_Mod; j++ ))
do

mod=${model[$j]}

echo $all_ex $N_proc $N_file $mod $zeta $mass
$all_ex $N_proc $N_file $mod $zeta $mass

done

done
