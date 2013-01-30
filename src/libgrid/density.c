#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../libmath/mathtools.h"
#include "../general_variables.h"
#include "../general_functions.h"

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


void init_density()
{
	int i=0;
	
	Density.N = Grid.N3;
	Density.r_n = (double *) calloc(Grid.N3, sizeof(double));
	Density.mass = (double *) calloc(Grid.N3, sizeof(double));
	Density.rho = (double *) calloc(Grid.N3, sizeof(double));

		for(i=0; i<Grid.N3; i++)
		{	
			Density.r_n[i] = (0.5 + i) * Grid.cell_volume;
			Density.mass[i] = Node[i].M;
			Density.rho[i] = Node[i].M / Grid.cell_volume;	
		}
}



void free_density()
{
	free(Density.r_n);
	free(Density.mass);
	free(Density.rho);
}



void find_density_maxima_and_minima()
{
	double MAX[5], MIN[5];
	int i=0, posMax[5], posMin[5];	

	maxima(Density.rho, Grid.N3, MAX, posMax, 5);
	minima(Density.rho, Grid.N3, MIN, posMin, 5);

		for(i=0; i<5; i++)
			fprintf(stdout, "MAX[%d] with M=%e at node %d, X=%lf, Y=%lf, Z=%lf\n", 
			i, MAX[i], posMax[i], Node[posMax[i]].X_cm[0], Node[posMax[i]].X_cm[1], Node[posMax[i]].X_cm[2]
				);
		
	INFO_MSG("\n");	

		for(i=0; i<5; i++)
			fprintf(stdout, "MIN[%d] with M=%e at node %d; X=%lf, Y=%lf, Z=%lf\n", 
			i, MIN[i], posMin[i], Node[posMin[i]].X_cm[0], Node[posMin[i]].X_cm[1], Node[posMin[i]].X_cm[2]
				);
}
