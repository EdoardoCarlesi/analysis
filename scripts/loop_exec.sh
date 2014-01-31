#!/bin/bash

home=${HOME}/Analysis/scripts/
out=${HOME}/Analysis/output/
all_ex=$home/all_execute.sh

model[0]='cde000'
model[1]='cde033'
model[2]='cde066'
model[3]='cde099'
model[4]='lcdm'

mass[0]=1.e+10
mass[1]=3.e+11
mass[2]=1.e+12
mass[3]=7.e+13
mass[4]=1.e+14
mass[5]=5.e+12

bins[0]=20
bins[1]=15
bins[2]=15
bins[3]=10
bins[4]=9
bins[5]=15

output[0]='1e10/'
output[1]='3e11/'
output[2]='1e12/'
output[3]='7e13/'
output[4]='1e14/'
output[5]='5e12/'

for (( i=3; i<4; i++ ))
do

mass=${mass[$i]}
out_dir=$out${output[$i]}
bin=${bins[$i]}

for (( j=0; j<5; j++ ))
do

mod=${model[$j]}

$all_ex 5 $mod $mass $bin
echo $all_ex 5 $mod $mass $bin

done

cd $out
mkdir ${output[$i]}
mkdir ${output[$i]}/data/
mv $out*dat $out_dir/data/

done
