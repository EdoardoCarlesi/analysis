#include <stdio.h>
#include <stdlib.h>

#include "subhalo_general.h"
#include "../general_variables.h"


void initialize_subhalo_properties_structure()
{
	int rBins, nBins;
	rBins=Settings.r_bins; nBins=Settings.n_bins;

	fprintf(stdout,"\nInitializing structure to store subhalo properties.\n");

		SubHaloZ.n_bins=nBins-1; nBins--;
		SubHaloZ.costh = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.costh_count = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.cosphi = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.cosphi_count = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.l = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.p_l = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.err_p_l = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.shape = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.p_shape = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.n_shape = (int*) calloc(nBins, sizeof(int));
		SubHaloZ.triax = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.p_triax = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.n_triax = (int*) calloc(nBins, sizeof(int));
		SubHaloZ.r_sub = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.n_r_sub = (int*) calloc(nBins, sizeof(double));
		SubHaloZ.cum_n_r_sub = (int*) calloc(nBins, sizeof(double));
		SubHaloZ.r_sub_subset = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.n_r_sub_subset = (int*) calloc(nBins, sizeof(double));
		SubHaloZ.cum_n_r_sub_subset = (int*) calloc(nBins, sizeof(double));
		SubHaloZ.vel_sub = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.n_vel_sub = (int*) calloc(nBins, sizeof(int));
		SubHaloZ.p_vel_sub = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.mass_sub = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.n_sub = (int*) calloc(nBins, sizeof(int));
		SubHaloZ.cum_n_sub = (int*) calloc(nBins, sizeof(int));
		SubHaloZ.ecc = (double*) calloc(nBins, sizeof(double));
		SubHaloZ.n_ecc = (int*) calloc(nBins, sizeof(int));
	
	fprintf(stdout,"\n");
}



void free_subhalo_properties()
{
	fprintf(stdout, "\nFreeing subhalo properties structure.\n");

		free(SubHaloZ.costh); 
		free(SubHaloZ.costh_count); 
		free(SubHaloZ.cosphi); 
		free(SubHaloZ.cosphi_count);
		free(SubHaloZ.l); 		
		free(SubHaloZ.p_l); 	
		free(SubHaloZ.err_p_l); 
		free(SubHaloZ.shape); 
		free(SubHaloZ.p_shape); 
		free(SubHaloZ.n_shape); 
		free(SubHaloZ.triax); 
		free(SubHaloZ.p_triax); 
		free(SubHaloZ.n_triax); 
		free(SubHaloZ.r_sub); 
		free(SubHaloZ.n_r_sub); 
		free(SubHaloZ.cum_n_r_sub); 
		free(SubHaloZ.r_sub_subset); 
		free(SubHaloZ.n_r_sub_subset); 
		free(SubHaloZ.cum_n_r_sub_subset); 
		free(SubHaloZ.vel_sub); 
		free(SubHaloZ.p_vel_sub); 
		free(SubHaloZ.n_vel_sub); 
		free(SubHaloZ.mass_sub); 
		free(SubHaloZ.n_sub); 
		free(SubHaloZ.cum_n_sub); 
		free(SubHaloZ.ecc); 
		free(SubHaloZ.n_ecc);
	
	fprintf(stdout,"\n");
}



void copy_subhalo_properties(int host, int n_sub)
{
//int index = find_subhalo_index(host);
// TODO for the moment just sort them in the subhaloes struct, then use find_subhalo_index to sort them properly
	int index=host;

	fprintf(stdout,"\nCopying (sub)halo properties into subhalo structures.\n");

		subhaloes[index].Rvir = haloes[n_sub].Rvir; 
		subhaloes[index].Vmax = haloes[n_sub].Vmax; 
                subhaloes[index].n_part = haloes[n_sub].n_part; 
                subhaloes[index].Xc = haloes[n_sub].Xc; 
                subhaloes[index].Yc = haloes[n_sub].Yc; 
                subhaloes[index].Zc = haloes[n_sub].Zc; 
                subhaloes[index].VXc = haloes[n_sub].VXc; 
                subhaloes[index].VYc = haloes[n_sub].VYc; 
                subhaloes[index].VZc = haloes[n_sub].VZc; 
                subhaloes[index].Eax = haloes[n_sub].Eax; 
                subhaloes[index].Eay = haloes[n_sub].Eay; 
                subhaloes[index].Eaz = haloes[n_sub].Eaz; 
                subhaloes[index].Mvir = haloes[n_sub].Mvir; 
                subhaloes[index].c = haloes[n_sub].c; 
                subhaloes[index].shape = haloes[n_sub].shape; 
                subhaloes[index].aa = haloes[n_sub].aa; 
                subhaloes[index].bb = haloes[n_sub].bb; 
                subhaloes[index].cc = haloes[n_sub].cc; 
                subhaloes[index].triax = haloes[n_sub].triax; 
                subhaloes[index].lambda = haloes[n_sub].lambda; 
                subhaloes[index].c_nfw = haloes[n_sub].c_nfw; 
                subhaloes[index].rs_nfw = haloes[n_sub].rs_nfw; 
                subhaloes[index].host = haloes[n_sub].host; 
                subhaloes[index].id = haloes[n_sub].id; 
                subhaloes[index].Epot = haloes[n_sub].Epot; 
                subhaloes[index].Ekin = haloes[n_sub].Ekin; 
                subhaloes[index].AngMom = haloes[n_sub].AngMom; 
                subhaloes[index].Jcirc = haloes[n_sub].Jcirc; 

	fprintf(stdout,"\n");
}



void avg_subhalo()
{
	int totSub=0, i=0;
	double avg_sub=0;

		for (i=0; i<Settings.n_threshold; i++) 
			totSub += haloes[i].n_satellites;

			avg_sub = (double) totSub /( (double) Settings.n_threshold );

		fprintf(stdout, "\nThe average number of subhaloes for hosts above the %e mass limit is: %lf\n", 
			Settings.mass_min, avg_sub);

		fprintf(stdout, "\nThere are %d subhaloes in total for the hosts above that mass limit.\n", 
			totSub);

	HaloZ.avgSub=avg_sub;
}



void init_subhalo_struct()
{
	int i=0, totSub=0, totSubOverN=0, totHost=Settings.n_threshold; 

		for (i=0; i< totHost; i++) 
			totSub += haloes[i].n_satellites;

		for(i=0; i<Settings.n_haloes; i++)
			if(haloes[i].host > -1 && haloes[i].host < totHost && haloes[i].n_part > Settings.n_min-1) 
				totSubOverN++;

			fprintf(stdout, "\nThere are %d sub haloes in total.\n", 
				totSub);

			fprintf(stdout, "\nThere are %d sub haloes in total with more than %d particles.\n", 
				totSubOverN, Settings.n_min);

	Settings.n_subhaloes = totSub;
	Settings.n_subhaloes_nmin = totSubOverN;
	subhaloes = (struct halo*) calloc(totSub, sizeof(struct halo));
}



void load_subhalo_list()
{
	int i=0, j=0;

		for(i=0; i<Settings.n_haloes; i++)
		{
			if(haloes[i].host < Settings.n_threshold && haloes[i].host > -1)
			{
				copy_subhalo_properties(j, haloes[i].id);
				j++;
			}
		}
}



int find_subhalo_index(int host)
{
		// TODO
		//int index=0; int i=0;
	return host;
}
