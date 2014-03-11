#!/bin/bash

tmp1=snapshot.1.tmp
tmp2=snapshot.2.tmp

cat $* > $tmp1
sort -b -nr -k 5 $tmp1 > $tmp2
rm -rf $tmp1 
