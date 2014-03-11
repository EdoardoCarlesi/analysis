#!/bin/bash

home=${HOME}/Analysis/scripts/
out=${HOME}/Analysis/output/
all_ex=$home/all_execute.sh

N_Mod=3
N_Out=1
N_Red=14

N_proc=4
N_file=256

model[0]='cde000'
model[1]='cde099'
model[2]='lcdm'

z[0]='0.000'
z[1]='0.008'
z[2]='0.027'
z[3]='0.046'
z[4]='0.065'
z[5]='0.085'
z[6]='0.105'
z[7]='0.126'
z[8]='0.148'
z[9]='0.170'
z[10]='0.241'
z[11]='0.503'
z[12]='0.697'
z[13]='1.006'

mass[0]=1.e+12
mass[1]=9.e+13
mass[2]=3.e+14
mass[4]=3.e+11
mass[5]=3.e+12

for (( i=0; i<$N_Out; i++ ))
do

mass=${mass[$i]}

for (( j=2; j<$N_Mod; j++ ))
do

mod=${model[$j]}

for (( k=0; k<$N_Red; k++ ))
do

zeta=${z[$k]}

echo $all_ex $N_proc $N_file $mod $zeta $mass
$all_ex $N_proc $N_file $mod $zeta $mass

done
done
done
