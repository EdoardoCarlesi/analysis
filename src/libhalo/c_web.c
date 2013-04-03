#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../libmath/math.h"
#include "../libcosmo/cosmo.h"
#include "../general_def.h"

#include "halo.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef WITH_MPI
#include "../libparallel/general.h"
#endif


/*
 * Declare functions
 */
#define Node(N, x, y, z) (x + N*y + N*N*z)

struct c_web
{
	double lambda[3];
	int type;
} *VWeb, *TWeb;


void compute_node_type(int);



/*
 * Initialize functions
 */ 
void read_v_web()
{
	int skip, i, NWeb;
	double a;
	char dummy[512];
	FILE *fweb = fopen(Urls.c_web_file, "r");

	NWeb = Settings.c_web_size;

	VWeb = malloc(pow3(NWeb) * sizeof(struct c_web));

	i = 0;
	skip = 1;

	while(i < pow3(NWeb))
	{
		fgets(dummy, 512, fweb);

		if(i >= skip) 
		{			  //1  2  3  4  5  6  7  8  9  10 
			sscanf(dummy, "%f %f %f %f %f %f %f %f %f %f", 
			&a, &a, &a, &a, VWeb[i-1].lambda[0], VWeb[i-1].lambda[1], VWeb[i-1].lambda[2]);

			compute_node_type(i-1);
			
			i++;
		}
	}
}



void assign_haloes_to_web()
{
	int i=0, j=0, m=0, index, ix, iy, iz, nHaloes, NWeb;
	double L, GridUnit;

	NWeb = Settings.c_web_size;
	L = Settings.box_size;
	GridUnit = NWeb / L;

	nHaloes = Settings.n_haloes;

		for(j=0; j<nHaloes; j++)
		{
			if(halo_condition(j) == 1)
			{
				ix = (int) GridUnit * Haloes[j].X[0]; 
				iy = (int) GridUnit * Haloes[j].X[1]; 
				iz = (int) GridUnit * Haloes[j].X[2]; 
	
					// Init all node types to zero
				for(i=0; i<4; i++)
					Haloes[j].web_type[i] = 0;

				index = VWeb[Node(NWeb, ix, iy, iz)].type;
				Haloes[j].web_type[index] = 1;
			}
		}

	free(VWeb);
}



void compute_node_type(int i)
{
	int j=0;
	double l1, l2, l3;

	l1 = VWeb[i].lambda[0];
	l2 = VWeb[i].lambda[1];
	l3 = VWeb[i].lambda[2];

		if(l1 > 0) j++;
		if(l2 > 0) j++;
		if(l3 > 0) j++;

	VWeb[i].type = j;
}
