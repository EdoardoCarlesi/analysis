#!/bin/bash

NP=2

cd ${HOME}/Analysis/src/; make clean; make test
cd ..

mpiexec -n $NP ./bin/test $NP
