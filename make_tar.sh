#!/bin/bash

make clean
mkdir Analysis
mkdir Analysis/temp
mkdir Analysis/output
mkdir Analysis/bin
cp -r src Analysis/
cp -r scripts Analysis/
cp Makefile* Analysis/
tar -cjf analysis.tar.gz Analysis
