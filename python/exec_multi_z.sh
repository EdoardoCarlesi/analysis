#!/bin/bash

multi_z='multi_z_analysis.py'
model='lcdm'
folder='/home/edoardo/Dropbox/Uni/UnsortedDataForPlots/HiResDeDm/all/data/'

ftype='-z'
suffix='avg_profiles'

output='../output/'$model'-'$suffix'_z.dat'
redshift=z[0-9].[0-9][0-9][0-9]

list_f='../temp/files.list'
list_z='../temp/zetas.list'

ls $folder$model*$ftype*$suffix*'.dat' > $list_f

grep -o $redshift $list_f | sed 's/z// ' > $list_z

./$multi_z $output $list_z $list_f $suffix
