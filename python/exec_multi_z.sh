#!/bin/bash

multi_z='multi_z_analysis.py'
model='lcdm'

sys='comodo'
#sys='castor'
#sys='laptop'

if [ $sys == "laptop" ]
then
folder='/home/edoardo/Dropbox/Uni/UnsortedDataForPlots/HiResDeDm/all/data/'
fi

if [ $sys == "comodo" ]
then
folder='/home/edoardo/Analysis/output/'
fi

suff[0]='avg_profiles'
suff[1]='all_halo'
suff[2]='numerical_mass'

ftype='-z'
redshift=z[0-9].[0-9][0-9][0-9]

list_f='../temp/files.list'
list_z='../temp/zetas.list'

ls $folder$model*$ftype*$suffix*'.dat' > $list_f

grep -o $redshift $list_f | sed 's/z// ' > $list_z

more $list_f
more $list_z

suffix=${suff[0]}
output='../output/'$model'-'$suffix'_z.dat'
./$multi_z $output $list_z $list_f $suffix $line_skip
