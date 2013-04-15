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
#define BIN_SIZE 12

struct c_web
{
	int type;
	float dens;
	float temp;
	float lambda[3];
} *VWeb, *TWeb;


struct web_info
{
	unsigned int N[5];	 // Total number per type of node, 0 is the total nodes
	unsigned int *ids[4];
	float volume[4]; // Volume fraction in node type
	float mass[4];	 // Mass fraction in node type

		// Overdensity values
	float delta[5][BIN_SIZE];
	float delta_dm[5][BIN_SIZE];
	float delta_gas[5][BIN_SIZE];
		// Temperature value at a given overdensity
	float T_vs_delta[5][BIN_SIZE];
		// Correlation coefficient between T and dm/gas/tot mass
	float alpha_T[5][BIN_SIZE];
	float p_alpha_T[5][BIN_SIZE];

} WebInfoDm, WebInfoGas, WebInfoTot;


int compute_node_type(struct c_web *);



/*
 * Initialize functions
 */ 
void read_v_web()
{
	int t=0, skip=0, i=0, NWeb=0;
	int type[4];
	float mass[4];
	float a, mass_norm;
	char dummy[512];
	FILE *fweb = fopen(Urls.c_web_file, "r");

	INFO_MSG("Reading V Web file");
	INFO_MSG(Urls.c_web_file);

	NWeb = Settings.c_web_size;

	VWeb = malloc(pow3(NWeb) * sizeof(struct c_web));
	fprintf(stderr, "Allocating %.2f MB for VWeb storage\n", (float)pow3(NWeb)*sizeof(struct c_web)/1024/1024);

	for(i=0; i<4; i++)
	{
		type[i] = 0;
		mass[i] = 0.0;
	}

	for(i=0; i<5; i++)
		WebInfoDm.ids[i] = malloc(sizeof(unsigned int));

	i = 0;
	skip = 1;

	for (i=0; i<pow3(NWeb)+skip; i++)
	{
		fgets(dummy, 512, fweb);

		if(i >= skip) 
		{	              //1  2  3  4  5  6  7 
			sscanf(dummy, "%f %f %f %f %f %f %f", 
					// Read into dummy variables
				       	&a,&a,&a,
					&VWeb[i-1].dens,
				       	&VWeb[i-1].lambda[0], &VWeb[i-1].lambda[1], &VWeb[i-1].lambda[2]
				);

			//fprintf(stderr, "dens=%f\n", VWeb[i-1].dens);
			//fprintf(stderr, "lamb1=%f\n", VWeb[i-1].lambda[0]);
			//fprintf(stderr, "lamb2=%f\n", VWeb[i-1].lambda[1]);
			t = compute_node_type(&VWeb[i-1]);
			mass[t]+=VWeb[i-1].dens;
			type[t]++;

			WebInfoDm.ids[t] = realloc(WebInfoDm.ids[t], type[t] * sizeof(unsigned int));
			WebInfoDm.ids[t][type[t]-1] = i;
		}
	}
	
	mass_norm = 0.0;

	for(i=0; i<4; i++)
		mass_norm += mass[i];

		WebInfoDm.N[0] = pow3(NWeb);

	for(i=0; i<4; i++)
	{	
		WebInfoDm.N[i+1] = type[i];	
		WebInfoDm.volume[i] = (float)type[i] / (float) pow3(NWeb);
		WebInfoDm.mass[i] = mass[i] / mass_norm;
		//fprintf(stdout, "type=%d, volume fraction=%f\n", i, WebInfoDm.volume[i]);
		//fprintf(stdout, "type=%d, entries =%f\n", i, type[i]);
		//fprintf(stdout, "type=%d, mass fraction=%f\n", i, WebInfoDm.mass[i]);
		//fprintf(stdout, "type=%d, mass =%f\n", i, mass[i]);
	}

	fclose(fweb);
}



void read_t_web()
{
	int t, skip, i, NWeb;
	int type[4];
	float mass[4];
	float mass_norm, a;
	char dummy[512];
	FILE *fweb = fopen(Urls.c_web_gas_file, "r");

	INFO_MSG("Reading T Web file");
	INFO_MSG(Urls.c_web_gas_file);

	NWeb = Settings.c_web_size;
	TWeb = malloc(pow3(NWeb) * sizeof(struct c_web));
	fprintf(stderr, "Allocating %.2f MB for TWeb storage\n", (float)pow3(NWeb)*sizeof(struct c_web)/1024/1024);

	for(i=0; i<4; i++)
	{
		type[i] = 0;
		mass[i] = 0.0; 
	}

	i = 0;
	skip = 1;

	for (i=0; i<pow3(NWeb)+skip; i++)
	{
		fgets(dummy, 512, fweb);

		if(i >= skip) 
		{	              //1  2  3  4  5  6  7  8
			sscanf(dummy, "%f %f %f %f %f %f %f %f", 
					// Read into dummy variables
				        &a,&a,&a,
					&TWeb[i-1].dens, 
				        &TWeb[i-1].lambda[0], &TWeb[i-1].lambda[1], &TWeb[i-1].lambda[2], 
				        &TWeb[i-1].temp
				);

			t = compute_node_type(&TWeb[i-1]);
			type[t]++;
			mass[t]+=TWeb[i-1].dens;
			TWeb[i-1].type=t;
		}
	}
		
	mass_norm = 0.0;

	for(i=0; i<4; i++)
		mass_norm += mass[i];

		WebInfoGas.N[0] = pow3(NWeb);

	for(i=0; i<4; i++)
	{		
		WebInfoGas.N[i+1] = type[i];	
		WebInfoGas.volume[i] = (float)type[i] / (float) pow3(NWeb);
		WebInfoGas.mass[i] = mass[i]/mass_norm;
	//	fprintf(stdout, "type=%d, volume fraction=%f\n", i, WebInfoGas.volume[i]);
	//	fprintf(stdout, "type=%d, entries=%f\n", i, type[i]);
	//	fprintf(stdout, "type=%d, mass fraction=%f\n", i, WebInfoGas.mass[i]);
	//	fprintf(stdout, "type=%d, mass =%f\n", i, mass[i]);
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
	GridUnit = (float) NWeb / L;

	nHaloes = Settings.n_haloes;

	for(i=0; i<4; i++)
		Settings.n_cweb_type[i] = 0;

//#		pragma omp parallel for			\
		private(j,ix,iy,iz,index)		\
		shared(nHaloes,GridUnit,Haloes,VWeb)
		for(j=0; j<nHaloes; j++)
		{
			if(halo_condition(j) == 1)
			{
				ix = (int) (GridUnit * Haloes[j].X[0]); 
				iy = (int) (GridUnit * Haloes[j].X[1]); 
				iz = (int) (GridUnit * Haloes[j].X[2]); 
	
					// Init all node types to minus one
				for(i=0; i<4; i++)
					Haloes[j].web_type[i] = -1;
					
					// Find the nearest node's type assuming that
					// the nodes in the Web file are already ordered
				index = VWeb[Node(NWeb, ix, iy, iz)].type;
				Haloes[j].web_type[index] = 1;
				Settings.n_cweb_type[index]++;
			//fprintf(stderr, "(%d , %d) type=%d, x=%f, y=%f z=%f; ix=%d iy=%d iz=%d\n",
			//j,Node(NWeb, ix, iy, iz),index,Haloes[j].X[0], Haloes[j].X[1], Haloes[j].X[2], ix, iy, iz);
			}
		}

	free(VWeb);
}



int compute_node_type(struct c_web *Web)
{
	int j=0;
	double l0, l1, l2, l3;
	
	l0 = Settings.l_web;
	l1 = Web->lambda[0];
	l2 = Web->lambda[1];
	l3 = Web->lambda[2];

	if(!(l1 > l2 && l2 > l3) && (l1 != 0.0))
		fprintf(stderr, 
	"Problem computing node type, incorrect ordering for the eigenvalues l1=%f, l2=%f, l3=%f\n", 
			l1, l2, l3);

		if(l3 > l0) j=3;
		if(l2 > l0 && l3 < l0) j=2;
		if(l1 > l0 && l2 < l0) j=1;
		if(l1 < l0) j=0;

	Web->type = j;

	//fprintf(stderr, "%d) NodeType=%d\n", i, j);
	
	return j;
}



void sort_web_statistics()
{
	int NWeb, NType, Ns, i=0, j=0, k=0, l=0;
	double w_gas, w_dm;
	double dmMax, dmMin, *dm_bin, *dm_bin_e; 
	double *dm_bin_T, *dm_bin_alpha, *dm_bin_gas;	
	double *delta_dm, *alpha_dm, *p_alpha_dm;
	double *delta_gas, *T_gas, *alpha_gas, *p_alpha_gas;
	double *delta_tot, *alpha_tot, *p_alpha_tot;

	NWeb = Settings.c_web_size;
	w_gas = 0.046 / 0.27;
	w_dm = 1 - w_gas;
	Ns = 5;

	// Now sort statistics per node type
	// i = 0 is the global statistics
	for(i=0; i<5; i++)
	{
		k = 0;
		l = 0;
		NType = WebInfoDm.N[i];
		
		fprintf(stderr, "Found %d nodes of type %d.\n", NType, i);
		dm_bin = malloc((BIN_SIZE+1)*sizeof(double));
		dm_bin_e = malloc((BIN_SIZE)*sizeof(double));
		dm_bin_T = malloc((BIN_SIZE)*sizeof(double));
		dm_bin_alpha = malloc((BIN_SIZE)*sizeof(double));
		dm_bin_gas = malloc((BIN_SIZE)*sizeof(double));

		fprintf(stderr, "Allocating %.2f MB for type analysis...", (double) Ns*NType*sizeof(double)/1024/1024);
		alpha_dm  = malloc(NType * sizeof(double));
		delta_dm  = malloc(NType * sizeof(double));
		delta_gas = malloc(NType * sizeof(double));
		delta_tot = malloc(NType * sizeof(double));
		T_gas = malloc(NType * sizeof(double));
		fprintf(stderr, "Done.\n");	 

#	pragma omp parallel for					\
	private(k, j) 						\
	shared(i, alpha_dm, delta_dm, delta_gas, delta_tot) 	\
	shared(TWeb, VWeb, NWeb, w_dm, w_gas, T_gas, WebInfoDm)	
	for(j=0; j<NType; j++)
	{
		if(i == 0)
		{
			delta_dm[j] =  VWeb[j].dens;	
			delta_gas[j] = TWeb[j].dens;	
			delta_tot[j] = w_gas * TWeb[j].dens + w_dm * VWeb[j].dens;	
			T_gas[j] = convert_u_to_T(TWeb[j].temp);
	
			if(delta_dm[j]!=0.0)
				alpha_dm[j] = T_gas[j] / delta_dm[j];

		}
		else
		{
			k = WebInfoDm.ids[i-1][j];

			delta_dm[j] = VWeb[k].dens;	
			delta_gas[j] = TWeb[k].dens;	
			delta_tot[j] = w_gas * TWeb[k].dens + w_dm * VWeb[k].dens;	
			T_gas[j] = convert_u_to_T(TWeb[k].temp);

			if(delta_dm[j]!=0.0)
				alpha_dm[j] = T_gas[j] / delta_dm[j];	

		}
	}	

			// Now bin into histograms
			dmMin = 0.999 * nonzero_minimum(delta_dm, NType);
			dmMax = 1.001 * maximum(delta_dm, NType);
	
		//fprintf(stderr, "min=%f, max=%f, dm[1]=%f\n", dmMin, dmMax, delta_dm[1]);

			dm_bin = log_stepper(dmMin, dmMax, BIN_SIZE+1);

			average_bin(delta_dm, T_gas, dm_bin, dm_bin_T, dm_bin_e, BIN_SIZE+1, NType);
			average_bin(delta_dm, alpha_dm, dm_bin, dm_bin_alpha, dm_bin_e, BIN_SIZE+1, NType);
			average_bin(delta_dm, delta_gas, dm_bin, dm_bin_gas, dm_bin_e, BIN_SIZE+1, NType);

		for(j=0; j<BIN_SIZE; j++)
		{
			//fprintf(stderr, "%d) dm[%d]=%f, T_bin[%d]=%f\n", i, j, dm_bin[j], j, dm_bin_T[j]);
			//fprintf(stderr, "%d) dm[%d]=%f, alpha_bin[%d]=%f\n", i, j, dm_bin[j], j, dm_bin_alpha[j]);
			//fprintf(stderr, "%d) dm[%d]=%f, gas_bin[%d]=%f\n", i, j, dm_bin[j], j, dm_bin_gas[j]);

			// Overdensity values
			WebInfoDm.delta[i][j] = 0.5 * (dm_bin[j] + dm_bin[j+1]);
			WebInfoDm.delta_gas[i][j] = dm_bin_gas[j];
			WebInfoDm.T_vs_delta[i][j] = dm_bin_T[j];
		}

		free(delta_gas);
		free(delta_tot);
		free(delta_dm);
		free(alpha_dm);
		free(T_gas);

		free(dm_bin);
		free(dm_bin_e);
		free(dm_bin_T);
		free(dm_bin_alpha);
	}
}



void print_web_statistics()
{
	int count=1, i=0, j=0;
	char out_url[200];
 	sprintf(out_url, "%s%.2f%s", Urls.output_prefix, Settings.l_web, "_cosmic_web_properties.dat");
	FILE *file_out = fopen(out_url, "w");	

	INFO_MSG("Print cosmic web statisics to file:");
	INFO_MSG(out_url);

		// Print some general info first
		fprintf(file_out, "# lambda_th: %f\n", Settings.l_web);

		fprintf(file_out, "# DM mass,    void: %f\t", WebInfoDm.mass[0]);
		fprintf(file_out, "sheet: %f\t", WebInfoDm.mass[1]);
		fprintf(file_out, "filament: %f\t", WebInfoDm.mass[2]);
		fprintf(file_out, "node: %f\n", WebInfoDm.mass[3]);

		fprintf(file_out, "# DM volume,  void: %f\t", WebInfoDm.volume[0]);
		fprintf(file_out, "sheet: %f\t", WebInfoDm.volume[1]);
		fprintf(file_out, "filament: %f\t", WebInfoDm.volume[2]);
		fprintf(file_out, "node: %f\n", WebInfoDm.volume[3]);

		fprintf(file_out, "# Gas mass,   void: %f\t", WebInfoGas.mass[0]);
		fprintf(file_out, "sheet: %f\t", WebInfoGas.mass[1]);
		fprintf(file_out, "filament: %f\t", WebInfoGas.mass[2]);
		fprintf(file_out, "node: %f\n", WebInfoGas.mass[3]);

		fprintf(file_out, "# Gas volume, void: %f\t", WebInfoGas.volume[0]);
		fprintf(file_out, "sheet: %f\t", WebInfoGas.volume[1]);
		fprintf(file_out, "filament: %f\t", WebInfoGas.volume[2]);
		fprintf(file_out, "node: %f\n", WebInfoGas.volume[3]);

		fprintf(file_out,"\n#");
		
		// Print header file
		for(j=0; j<5; j++)
		{
			fprintf(file_out,"dDM type=%d (%d)\t", j, count++);
			fprintf(file_out,"Gas type=%d (%d)\t", j, count++);
			fprintf(file_out,"T   type=%d (%d)\t", j, count++);
		}
				
			fprintf(file_out,"\n");

			// Print it all
			for(i=0; i<BIN_SIZE; i++)
			{
				for(j=0; j<5; j++)
				{
					fprintf(file_out, "%f\t", WebInfoDm.delta[j][i]);
					fprintf(file_out, "%f\t", WebInfoDm.delta_gas[j][i]);
					fprintf(file_out, "%f\t", WebInfoDm.T_vs_delta[j][i]);
				}

				fprintf(file_out,"\n");
			}

		INFO_MSG("Done");	

	fclose(file_out);
}
