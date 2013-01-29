#/bin/bash

echo 'MPI checking if settings are consistent'

mkfile=$1/Makefile.config
tmp=$1/Makefile.config.tmp

# Mpi is set
if [ $2 -eq 1 ]
then

if grep '#WITH_MPI=' $mkfile
then
sed 's/#WITH_MPI=/WITH_MPI=/' <$mkfile >$tmp
mv $tmp $mkfile
fi 

else # Mpi option is not chosen

if grep '^WITH_MPI=' $mkfile
then
sed 's/WITH_MPI=/#WITH_MPI=/' <$mkfile >$tmp
mv $tmp $mkfile
fi

fi
