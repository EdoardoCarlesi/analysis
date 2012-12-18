#!/bin/bash
# Script to generate Pk_s using amigaPk #
#########################################

swap=$1
amigadir=$2
filesdir=$3
snap_zero=$4
tot_snaps=$5
grid=$6
list_file=$7


if [ $swap -eq 1 ] ; then
amigaPk=$amigadir'/bin/amigaPk-swap'
else 
#amigaPk=$amigadir'/bin/amigaPk'
amigaPk=$amigadir'/bin/simuPk'
fi

cd $filesdir

for (( i=$snap_zero; i<=$tot_snaps; i++))

do

if [ $i -gt 9 ] ; then

snap="snapshot_0$i"

else 

snap="snapshot_00$i"

fi

if [ -e $snap ] ; then
$amigaPk $snap $grid 1
else
echo 'File ' $snap ' does not exist in ' $filesdir
echo ' '
fi

done

echo "P(k) done!"

exit 0
