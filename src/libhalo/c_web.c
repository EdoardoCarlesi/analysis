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
	float lambda[3];
	int type;
} *VWeb, *TWeb;


int compute_node_type(int);



/*
 * Initialize functions
 */ 
void read_v_web()
{
	int t, skip, i, NWeb;
	int type[4];
	double a;
	char dummy[512];
	FILE *fweb = fopen(Urls.c_web_file, "r");

	INFO_MSG("Reading V Web file");
	INFO_MSG(Urls.c_web_file);

	NWeb = Settings.c_web_size;

	VWeb = malloc(pow3(NWeb) * sizeof(struct c_web));

	for(i=0; i<4; i++)
		type[i] = 0;

	i = 0;
	skip = 1;

	for (i=0; i<pow3(NWeb)+1; i++)
	{
		fgets(dummy, 512, fweb);

		if(i >= skip) 
		{	              //1  2  3  4  5  6  7  8  9  10 
			sscanf(dummy, "%f %f %f %f %f %f %f %f %f %f", 
				       &a,&a,&a,&a,&a,&a,&a, 
				       &VWeb[i-1].lambda[0], &VWeb[i-1].lambda[1], &VWeb[i-1].lambda[2]);

			t = compute_node_type(i-1);
			type[t]++;
		}
	}
	
	for(i=0; i<4; i++)
	{		
		fprintf(stdout, "type=%d, fraction=%d\n", i, type[i]);
		fprintf(stdout, "type=%d, fraction=%f\n", i, (float)type[i] / (float) pow3(NWeb));
	}

	fclose(fweb);
}



void assign_haloes_to_web()
{
	int i=0, j=0, m=0, index, ix, iy, iz, nHaloes, NWeb;
	double L, GridUnit;

	INFO_MSG("Assigning haloes to V Web nodes");

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

				fprintf(stderr, "(%d , %d) type=%d, x=%f, y=%f z=%f; ix=%d iy=%d iz=%d\n",
			j,Node(NWeb, ix, iy, iz),index,Haloes[j].X[0], Haloes[j].X[1], Haloes[j].X[2], ix, iy, iz);
			}
		}

	free(VWeb);
}



int compute_node_type(int i)
{
	int j=0;
	double l0, l1, l2, l3;
	
	l0 = Settings.l_web;
	l1 = VWeb[i].lambda[0];
	l2 = VWeb[i].lambda[1];
	l3 = VWeb[i].lambda[2];


	if(l1 > l2 && l2 > l3)
	{
	} else {
		fprintf(stderr, "incorrect ordering l1=%f, l2=%f, l3=%f\n", l1, l2, l3);
	}

/*
		if(l1 > l0) j++;
		if(l2 > l0) j++;
		if(l3 > l0) j++;
*/		

		if(l3 > l0) j=3;
		if(l2 > l0 && l3 < l0) j=2;
		if(l1 > l0 && l2 < l0) j=1;
		if(l1 < l0) j=0;
		//else j=0;

	VWeb[i].type = j;

	fprintf(stderr, "%d) NodeType=%d\n", i, j);
	
	return j;
}
