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

void stdout_halo_status(void);


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

			if(halo_condition(i) == 1)
				Settings.n_sub_threshold++;

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



void stdout_halo_status()
{
	fprintf(stdout, "%lf",   HaloProperties[HALO_INDEX].z);   
	fprintf(stdout, "\t%lf", HaloProperties[HALO_INDEX].l_0);   
	fprintf(stdout, "\t%lf", HaloProperties[HALO_INDEX].halo.s0);   
	fprintf(stdout, "\t%lf", HaloProperties[HALO_INDEX].halo.t0);   
	fprintf(stdout, "\t%lf", HaloProperties[HALO_INDEX].avgSub);   
	fprintf(stdout, "\t%lf", SubHaloProperties[HALO_INDEX].l_0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[HALO_INDEX].halo.s0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[HALO_INDEX].halo.t0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[HALO_INDEX].costh0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[HALO_INDEX].cosphi0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[HALO_INDEX].vel_0);   
	fprintf(stdout, "\t%e",  SubHaloProperties[HALO_INDEX].avgMass);   
	fprintf(stdout, "\n");
}



void initialize_halo_storage()
{
	int k=0, j=0, nTot=0;

	nTot = Urls.nCatalogueFiles;

	HaloProperties = (struct halo_properties *) calloc(nTot, sizeof(struct halo_properties));
	SubHaloProperties = (struct halo_properties *) calloc(nTot, sizeof(struct halo_properties));
}



void initialize_halo_properties_structure()
{
	int rBins, nBins;

	rBins=Settings.r_bins-1; 
	nBins=Settings.n_bins-1;

			// Halo axis alignment
		HaloProperties[HALO_INDEX].r_bins=rBins; 
		HaloProperties[HALO_INDEX].R = (double*) calloc(rBins, sizeof(double));
		HaloProperties[HALO_INDEX].Th_p = (double*) calloc(rBins, sizeof(double));
		HaloProperties[HALO_INDEX].Th_c = (double*) calloc(rBins, sizeof(double));
		HaloProperties[HALO_INDEX].N_pairs = (int*) calloc(rBins, sizeof(double));

			// General halo properties
		HaloProperties[HALO_INDEX].n_bins=nBins; 
		HaloProperties[HALO_INDEX].c = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].p_c = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].halo.l = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].halo.p_l = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].halo.s = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].halo.p_s = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].halo.t = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].halo.p_t = (double*) calloc(nBins, sizeof(double));

		HaloProperties[HALO_INDEX].p_fit_nfw.chi = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].p_fit_nfw.p_chi = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].p_fit_nfw.per = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].p_fit_nfw.p_per = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].p_fit_nfw.gof = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].p_fit_nfw.p_gof = (double*) calloc(nBins, sizeof(double));

		HaloProperties[HALO_INDEX].mass = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].n_entry = (int*) calloc(nBins, sizeof(int));
		HaloProperties[HALO_INDEX].vel = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].conc = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].halo.shape = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].halo.triax = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].halo.lambda = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].halo.virial = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].fit_nfw.chi = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].fit_nfw.per = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].fit_nfw.gof = (double*) calloc(nBins, sizeof(double));

			// General subhalo properties
		SubHaloProperties[HALO_INDEX].n_bins=nBins; 
		SubHaloProperties[HALO_INDEX].c = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].p_c = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].halo.l = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].halo.p_l = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].halo.s = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].halo.p_s = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].halo.t = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].halo.p_t = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].mass = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_entry = (int*) calloc(nBins, sizeof(int));
		SubHaloProperties[HALO_INDEX].vel = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].conc = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].halo.shape = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].halo.triax = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].halo.lambda = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].halo.virial = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.shape = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.triax = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.lambda = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].fit_nfw.chi = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].fit_nfw.per = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].fit_nfw.gof = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].p_fit_nfw.chi = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].p_fit_nfw.p_chi = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].p_fit_nfw.per = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].p_fit_nfw.p_per = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].p_fit_nfw.gof = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].p_fit_nfw.p_gof = (double*) calloc(nBins, sizeof(double));

			// Sub halo properties only
		SubHaloProperties[HALO_INDEX].r_sub = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_r_sub = (int*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].cum_n_r_sub = (int*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].r_sub_subset = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_r_sub_subset = (int*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].cum_n_r_sub_subset = (int*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].vel_sub = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_vel_sub = (int*) calloc(nBins, sizeof(int));
		SubHaloProperties[HALO_INDEX].p_vel_sub = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_sub = (int*) calloc(nBins, sizeof(int));
		SubHaloProperties[HALO_INDEX].cum_n_sub = (int*) calloc(nBins, sizeof(int));
		SubHaloProperties[HALO_INDEX].ecc = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].n_ecc = (int*) calloc(nBins, sizeof(int));
		SubHaloProperties[HALO_INDEX].costh = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].costh_count = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].cosphi = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].cosphi_count = (double*) calloc(nBins, sizeof(double));
#ifdef GAS
		HaloProperties[HALO_INDEX].T = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].n_T = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].cm = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].p_cm = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas_dm_cth = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].p_gas_dm_cth = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.shape = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.triax = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.lambda = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.beta = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.b = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.p_b = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.l = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.p_l = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.s = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.p_s = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.t = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.p_t = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.virial = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas.ekin = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].dm.virial = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].dm.ekin = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].diff.l = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].diff.p_l = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].diff.s = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].diff.p_s = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].diff.t = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].diff.p_t = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas_T = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas_u = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas_fraction = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas_dm_costh = (double*) calloc(nBins, sizeof(double));
		HaloProperties[HALO_INDEX].gas_diff_cm = (double*) calloc(nBins, sizeof(double));

		SubHaloProperties[HALO_INDEX].cm = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].p_cm = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas_dm_cth = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].p_gas_dm_cth = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.virial = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.ekin = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].dm.virial = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].dm.ekin = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.shape = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.triax = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.lambda = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas_T = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas_u = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas_fraction = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].diff.l = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].diff.p_l = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].diff.s = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].diff.p_s = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].diff.t = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].diff.p_t = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.l = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.p_l = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.s = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.p_s = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.t = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas.p_t = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas_fraction = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas_dm_costh = (double*) calloc(nBins, sizeof(double));
		SubHaloProperties[HALO_INDEX].gas_diff_cm = (double*) calloc(nBins, sizeof(double));
#endif
}



void free_halo_properties()
{
/*
		free(SubHaloProperties[HALO_INDEX].costh); 
		free(SubHaloProperties[HALO_INDEX].costh_count); 
		free(SubHaloProperties[HALO_INDEX].cosphi); 
		free(SubHaloProperties[HALO_INDEX].cosphi_count);
		free(SubHaloProperties[HALO_INDEX].l); 		
		free(SubHaloProperties[HALO_INDEX].p_l); 	
		free(SubHaloProperties[HALO_INDEX].shape); 
		free(SubHaloProperties[HALO_INDEX].p_s); 
		free(SubHaloProperties[HALO_INDEX].triax); 
		free(SubHaloProperties[HALO_INDEX].p_t); 
		free(SubHaloProperties[HALO_INDEX].r_sub); 
		free(SubHaloProperties[HALO_INDEX].n_r_sub); 
		free(SubHaloProperties[HALO_INDEX].cum_n_r_sub); 
		free(SubHaloProperties[HALO_INDEX].r_sub_subset); 
		free(SubHaloProperties[HALO_INDEX].n_r_sub_subset); 
		free(SubHaloProperties[HALO_INDEX].cum_n_r_sub_subset); 
		free(SubHaloProperties[HALO_INDEX].vel_sub); 
		free(SubHaloProperties[HALO_INDEX].p_vel_sub); 
		free(SubHaloProperties[HALO_INDEX].n_vel_sub); 
		free(SubHaloProperties[HALO_INDEX].mass); 
		free(SubHaloProperties[HALO_INDEX].n_sub); 
		free(SubHaloProperties[HALO_INDEX].cum_n_sub); 
		free(SubHaloProperties[HALO_INDEX].ecc); 
		free(SubHaloProperties[HALO_INDEX].n_ecc);
	
		free(HaloProperties[HALO_INDEX].c);
		free(HaloProperties[HALO_INDEX].p_c);
		free(HaloProperties[HALO_INDEX].l);
		free(HaloProperties[HALO_INDEX].p_l);
		free(HaloProperties[HALO_INDEX].p_s);
		free(HaloProperties[HALO_INDEX].shape);
		free(HaloProperties[HALO_INDEX].triax);
		free(HaloProperties[HALO_INDEX].p_t);
		free(HaloProperties[HALO_INDEX].mass);
		free(HaloProperties[HALO_INDEX].vel);
#ifdef GAS
		free(HaloProperties[HALO_INDEX].gas_T);
		free(HaloProperties[HALO_INDEX].gas_u);
		free(HaloProperties[HALO_INDEX].gas_fraction);
		free(SubHaloProperties[HALO_INDEX].gas_T);
		free(SubHaloProperties[HALO_INDEX].gas_u);
		free(SubHaloProperties[HALO_INDEX].gas_fraction);
#endif 
		free(Haloes);
*/
}



void list_halo_sample(int *index)
{
	int i=0, m=0;

	for(i=0; i<Settings.n_haloes; i++)
	{
		if(halo_condition(i) == 1)
		{
			index[m] = i;
			m++;
		}
	}
}



void free_halo_profiles()
{
	int i=0;
	struct halo *HALO;
	struct general_settings *SETTINGS;	

	INFO_MSG("Freeing halo profiles");

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
			free(HALO[i].mass_r);
			free(HALO[i].npart);
#ifdef GAS
			free(HALO[i].gas_only.u);
			free(HALO[i].gas_only.m);
#endif
		}
}
