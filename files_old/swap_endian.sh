#!/bin/bash
# Script to convert Gadget files' endianness
#

echo "Convert endianness"

cd ../swapgadget
make clean; make

echo "LCDM (1), DEDM (2) or VDE(3)?"
read model

echo "Initial Snapshot: "
read is

echo "To Snapshot: "
read ifa

# Home directory /home/edoardo
homedir=${directories[0]}

# Where the gadget snapshots are stored
datadir=${directories[1]}

# Where the Analysis program is found 
analysisdir=${directories[2]}

# Subdirectories for the different models
lcdmdir=$datadir${directories[3]}
dedmdir=$datadir${directories[4]}
vdedir=$datadir${directories[5]}

if [ $model -eq 1  ] ; then
cd $lcdmdir
fi

if [ $model -eq 2 ] ; then
cd $dedmdir
fi

if [ $model -eq 3 ] ; then
cd $vdedir
fi

for (( i=is; i<=ifa; i++))

do

if [ $i -gt 9 ] ; then

snap="snapshot_0$i"

else 

snap="snapshot_00$i"

fi

../bin/gadget_swap -snap $snap
done

exit 0
