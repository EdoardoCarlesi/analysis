#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "../libio/io.h"
#include "../general_def.h"

#include "halo.h"

#ifdef WITH_MPI
#include "../libparallel/general.h"
#endif


/*
 * Declare functions
 */
int HALO_INDEX;

int find_host_index(uint64_t);

void stdout_halo_status(int);


/*
 * Initialize functions
 */
void find_substructure()
{
	int i=0, host=0;

	INFO_MSG("Sorting halo substructure");

	SubStructure.N_host = -1;
	SubStructure.N_sub = -1;

	SubStructure.sub = calloc(1, sizeof(struct sub_halo));
	SubStructure.host = calloc(1, sizeof(struct host_halo));

	for(i=0; i<Settings.n_haloes; i++)
	{
		if(Haloes[i].n_satellites > 0)
		{
			SubStructure.N_host++;
			SubStructure.host = realloc(SubStructure.host, (SubStructure.N_host+1)*sizeof(struct host_halo));
			SubStructure.host[SubStructure.N_host].id = Haloes[i].id;
			SubStructure.host[SubStructure.N_host].index = i;

			SubStructure.host[SubStructure.N_host].n_sub = -1;
			SubStructure.host[SubStructure.N_host].sub_index = calloc(1, sizeof(int));

		//	fprintf(stderr, "%d) Host id=%llu, index=%d, nsub=%d\n", 
		//		i, Haloes[i].id, i, Haloes[i].n_satellites);

		}
		
		if(Haloes[i].host > 0)
		{
			host = find_host_index(Haloes[i].host);

			SubStructure.N_sub++;
			SubStructure.sub = realloc(SubStructure.sub, (SubStructure.N_sub+1)*sizeof(struct sub_halo));
			SubStructure.sub[SubStructure.N_sub].id = Haloes[i].id;
			SubStructure.sub[SubStructure.N_sub].index = i;

			SubStructure.sub[SubStructure.N_sub].host_id = Haloes[i].host;
			SubStructure.sub[SubStructure.N_sub].host_index = host;

			SubStructure.host[host].n_sub++;
			SubStructure.host[host].sub_index = 
			realloc(SubStructure.host[host].sub_index, (SubStructure.host[host].n_sub+1) * sizeof(int));
			SubStructure.host[host].sub_index[SubStructure.host[host].n_sub] = i;

	//		fprintf(stderr, "%d) Sub id=%llu, index=%d, host_id=%llu, host_index=%d halo_index=%d\n", 
	//			i, Haloes[i].id, i, Haloes[i].host, host, SubStructure.host[host].index);
		}
	
	}

		fprintf(stdout, "\nFound %d haloes with a total of %d subhaloes\n", 
			SubStructure.N_host, SubStructure.N_sub);

}



int find_host_index(uint64_t host_id)
{
	int i=0, host_index=0;

		for(i=0; i<SubStructure.N_host; i++)
			if(SubStructure.host[i].id == host_id)
				host_index = i;

	return host_index;
}


void compute_halo_and_subhalo_statistics()
{
	fprintf(stdout, "Reading halo url[%d]: %s\n", HALO_INDEX, Urls.urls[HALO_INDEX]);

	Settings.use_cat = HALO_INDEX;
	set_halo_url();

		read_halo_file();

		compute_halo_properties();
		compute_subhalo_properties();

	stdout_halo_status(HALO_INDEX);

}



void stdout_halo_status(int j)
{
	fprintf(stdout, "%lf",   HaloProperties[j].z);   
	fprintf(stdout, "\t%lf", HaloProperties[j].c_0);   
	fprintf(stdout, "\t%lf", HaloProperties[j].l_0);   
	fprintf(stdout, "\t%lf", HaloProperties[j].s0);   
	fprintf(stdout, "\t%lf", HaloProperties[j].t0);   
	fprintf(stdout, "\t%lf", HaloProperties[j].avgSub);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].c_0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].l_0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].s0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].t0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].costh0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].cosphi0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].vel_0);   
	fprintf(stdout, "\t%e",  SubHaloProperties[j].avgMass);   
	fprintf(stdout, "\n");
}



void initialize_halo_storage()
{
	int k=0, j=0, nTot=0;

	nTot = Urls.nCatalogueFiles;

	HaloProperties = (struct halo_properties *) calloc(nTot, sizeof(struct halo_properties));
	SubHaloProperties = (struct halo_properties *) calloc(nTot, sizeof(struct halo_properties));

		for(j=0; j<nTot; j++)
		{
			k = GrowthFac.npts - j - 1;
			HaloProperties[j].z = GrowthFac.z[k];
			SubHaloProperties[j].z = GrowthFac.z[k];
		}
}



void initialize_halo_properties_structure()
{
	int rBins, nBins;

	rBins=Settings.r_bins; 
	nBins=Settings.n_bins;

			// Halo axis alignment
		HaloProperties[HALO_INDEX].r_bins=rBins-1; rBins--;
		HaloProperties[HALO_INDEX].R = (double*) calloc(rBins, sizeof(double));
		HaloProperties[HALO_INDEX].Th_p = (double*) calloc(rBins, sizeof(double));
		HaloProperties[HALO_INDEX].Th_c = (double*) calloc(rBins, sizeof(double));
		HaloProperties[HALO_INDEX].N_pairs = (int*) calloc(rBins, sizeof(double));
			// Other halo properties
		HaloProperties[HALO_INDEX].n_bins=nBins-1; nBins--;
		HaloProperties[HALO_INDEX].c = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].c_avg = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].p_c = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].l = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].p_l = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].err_p_l = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].shape = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].p_shape = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].n_shape = (int*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].triax = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].p_triax = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].n_triax = (int*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].mass = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].radVel = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].err_radVel = (double*) calloc(nBins, sizeof(double));

		SubHaloProperties[HALO_INDEX].n_bins=nBins-1; nBins--;
		SubHaloProperties[HALO_INDEX].costh = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].costh_count = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].cosphi = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].cosphi_count = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].l = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].p_l = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].err_p_l = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].shape = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].p_shape = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_shape = (int*) calloc(nBins, sizeof(int));
		SubHaloProperties[HALO_INDEX].triax = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].p_triax = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_triax = (int*) calloc(nBins, sizeof(int));
		SubHaloProperties[HALO_INDEX].r_sub = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_r_sub = (int*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].cum_n_r_sub = (int*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].r_sub_subset = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_r_sub_subset = (int*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].cum_n_r_sub_subset = (int*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].vel_sub = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_vel_sub = (int*) calloc(nBins, sizeof(int));
		SubHaloProperties[HALO_INDEX].p_vel_sub = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].mass_sub = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_sub = (int*) calloc(nBins, sizeof(int));
		SubHaloProperties[HALO_INDEX].cum_n_sub = (int*) calloc(nBins, sizeof(int));
		SubHaloProperties[HALO_INDEX].ecc = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_ecc = (int*) calloc(nBins, sizeof(int));
#ifdef GAS
		HaloProperties[HALO_INDEX].gas_T = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas_u = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas_fraction = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas_T = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas_u = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas_fraction = (double*) calloc(nBins, sizeof(double));
#endif
}



void free_halo_properties()
{
		free(SubHaloProperties[HALO_INDEX].costh); 
		free(SubHaloProperties[HALO_INDEX].costh_count); 
		free(SubHaloProperties[HALO_INDEX].cosphi); 
		free(SubHaloProperties[HALO_INDEX].cosphi_count);
		free(SubHaloProperties[HALO_INDEX].l); 		
		free(SubHaloProperties[HALO_INDEX].p_l); 	
		free(SubHaloProperties[HALO_INDEX].err_p_l); 
		free(SubHaloProperties[HALO_INDEX].shape); 
		free(SubHaloProperties[HALO_INDEX].p_shape); 
		free(SubHaloProperties[HALO_INDEX].n_shape); 
		free(SubHaloProperties[HALO_INDEX].triax); 
		free(SubHaloProperties[HALO_INDEX].p_triax); 
		free(SubHaloProperties[HALO_INDEX].n_triax); 
		free(SubHaloProperties[HALO_INDEX].r_sub); 
		free(SubHaloProperties[HALO_INDEX].n_r_sub); 
		free(SubHaloProperties[HALO_INDEX].cum_n_r_sub); 
		free(SubHaloProperties[HALO_INDEX].r_sub_subset); 
		free(SubHaloProperties[HALO_INDEX].n_r_sub_subset); 
		free(SubHaloProperties[HALO_INDEX].cum_n_r_sub_subset); 
		free(SubHaloProperties[HALO_INDEX].vel_sub); 
		free(SubHaloProperties[HALO_INDEX].p_vel_sub); 
		free(SubHaloProperties[HALO_INDEX].n_vel_sub); 
		free(SubHaloProperties[HALO_INDEX].mass_sub); 
		free(SubHaloProperties[HALO_INDEX].n_sub); 
		free(SubHaloProperties[HALO_INDEX].cum_n_sub); 
		free(SubHaloProperties[HALO_INDEX].ecc); 
		free(SubHaloProperties[HALO_INDEX].n_ecc);
	
		free(HaloProperties[HALO_INDEX].c);
		free(HaloProperties[HALO_INDEX].p_c);
		free(HaloProperties[HALO_INDEX].l);
		free(HaloProperties[HALO_INDEX].p_l);
		free(HaloProperties[HALO_INDEX].err_p_l);
		free(HaloProperties[HALO_INDEX].p_shape);
		free(HaloProperties[HALO_INDEX].n_shape);
		free(HaloProperties[HALO_INDEX].shape);
		free(HaloProperties[HALO_INDEX].triax);
		free(HaloProperties[HALO_INDEX].p_triax);
		free(HaloProperties[HALO_INDEX].n_triax);
		free(HaloProperties[HALO_INDEX].mass);
		free(HaloProperties[HALO_INDEX].radVel);
		free(HaloProperties[HALO_INDEX].err_radVel);
#ifdef GAS
		free(HaloProperties[HALO_INDEX].gas_T);
		free(HaloProperties[HALO_INDEX].gas_u);
		free(HaloProperties[HALO_INDEX].gas_fraction);
		free(SubHaloProperties[HALO_INDEX].gas_T);
		free(SubHaloProperties[HALO_INDEX].gas_u);
		free(SubHaloProperties[HALO_INDEX].gas_fraction);
#endif 
		free(Haloes);
}



void free_halo_profiles()
{
	int i=0;
	struct halo *HALO;
	struct general_settings *SETTINGS;	

#ifdef WITH_MPI
		HALO = pHaloes[ThisTask];
		SETTINGS = &pSettings[ThisTask];
#else
		HALO = Haloes;
		SETTINGS = &Settings;
#endif
		for(i=0; i<SETTINGS->n_haloes; i++)
		{
			free(HALO[i].radius);
			free(HALO[i].rho);
			free(HALO[i].err);
#ifdef WITH_GAS
			free(HALO[i].m_gas);
			free(HALO[i].u_gas);
#endif
		}
}
