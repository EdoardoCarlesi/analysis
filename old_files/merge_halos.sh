#!/bin/bash
# Combine into a single catalog the haloes found
# during the analysis by the different CPUs
# using AHF 

echo 'Merge halo catalogues files into a single snapshot'

if [ $1 == 1 ] ; then

totSnaps=$2
totCPUs=$3
model=$4
box=$5
dirBase=${HOME}/data/dedm/$box/$model/ahf/

ahf_merge=${HOME}'/amiga_old/bin/ahf_merge'
baseName="snapshot_0" 
baseNameOut=${HOME}'/data/dedm/'$box/$model'/MERGED/merged_snapshot_'

else

totSnaps=30
totCPUs=2
particles=128
model=cde3
box=50
dirBase=${HOME}/data/dedm/$particles/$box/$model/ahf/

ahf_merge='/home/edoardo/ahftools/tools/ahf_merge'
baseName="snapshot_0" 
baseNameOut=${HOME}'/data/dedm/'/$particles/$box/$model'/MERGED/merged_snapshot_'

fi

cd $dirBase

for((i=0; i<totSnaps; i++))

do

snapNumber=$i

if [ $i -gt 9 ] ; then

fileName=$baseName$snapNumber
num=$i

else

fileName=$baseName"0"$snapNumber
num=0$i

fi

for((k=0; k<totCPUs; k++))

do

if [ $k -gt 9 ] ; then

fileNameA=$fileName".00"$k

else 

fileNameA=$fileName".000"$k

fi

done

# we generate a file containing all the "names" which will have to be used 

ls $fileNameA* > new
sed 's/.AHF_profiles//' <new >new0
sed 's/.AHF_halos//' <new0 >new1
sed 's/.AHF_particles//' <new1 >new2
sed "s/"$fileNameA."//" <new2 >new3
sed "s/\./ /" <new3 >new4
sed 's/z//' <new4 >new5

exec < new5
read a1 a2

echo $ahf_merge $fileName 4 z$a1.$a2 $totCPUs $baseNameOut"_"$num"_z"$a1
$ahf_merge $fileName 4 z$a1.$a2 $totCPUs $baseNameOut"_"$num"_z"$a1
rm -rf new0 new new1 new2 new3 new4 new5

done 

exit 0
