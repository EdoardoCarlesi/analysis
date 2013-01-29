#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "density.h"
#include "grid.h"

#ifdef WITH_MPI
#include <mpi.h>
#include "../libparallel/general.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


struct density Density;
