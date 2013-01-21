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

		SubHaloes[index].Rvir = Haloes[n_sub].Rvir; 
		SubHaloes[index].Vmax = Haloes[n_sub].Vmax; 
                SubHaloes[index].n_part = Haloes[n_sub].n_part; 
                SubHaloes[index].Xc = Haloes[n_sub].Xc; 
                SubHaloes[index].Yc = Haloes[n_sub].Yc; 
                SubHaloes[index].Zc = Haloes[n_sub].Zc; 
                SubHaloes[index].VXc = Haloes[n_sub].VXc; 
                SubHaloes[index].VYc = Haloes[n_sub].VYc; 
                SubHaloes[index].VZc = Haloes[n_sub].VZc; 
                SubHaloes[index].Eax = Haloes[n_sub].Eax; 
                SubHaloes[index].Eay = Haloes[n_sub].Eay; 
                SubHaloes[index].Eaz = Haloes[n_sub].Eaz; 
                SubHaloes[index].Mvir = Haloes[n_sub].Mvir; 
                SubHaloes[index].c = Haloes[n_sub].c; 
                SubHaloes[index].shape = Haloes[n_sub].shape; 
                SubHaloes[index].aa = Haloes[n_sub].aa; 
                SubHaloes[index].bb = Haloes[n_sub].bb; 
                SubHaloes[index].cc = Haloes[n_sub].cc; 
                SubHaloes[index].triax = Haloes[n_sub].triax; 
                SubHaloes[index].lambda = Haloes[n_sub].lambda; 
                SubHaloes[index].c_nfw = Haloes[n_sub].c_nfw; 
                SubHaloes[index].rs_nfw = Haloes[n_sub].rs_nfw; 
                SubHaloes[index].host = Haloes[n_sub].host; 
                SubHaloes[index].id = Haloes[n_sub].id; 
                SubHaloes[index].Epot = Haloes[n_sub].Epot; 
                SubHaloes[index].Ekin = Haloes[n_sub].Ekin; 
                SubHaloes[index].AngMom = Haloes[n_sub].AngMom; 
                SubHaloes[index].Jcirc = Haloes[n_sub].Jcirc; 

	fprintf(stdout,"\n");
}



void avg_subhalo()
{
	int totSub=0, i=0;
	double avg_sub=0;

		for (i=0; i<Settings.n_threshold; i++) 
			totSub += Haloes[i].n_satellites;

			avg_sub = (double) totSub /( (double) Settings.n_threshold );

		fprintf(stdout, "\nThe average number of SubHaloes for hosts above the %e mass limit is: %lf\n", 
			Settings.mass_min, avg_sub);

		fprintf(stdout, "\nThere are %d SubHaloes in total for the hosts above that mass limit.\n", 
			totSub);

	HaloZ.avgSub=avg_sub;
}



void init_subhalo_struct()
{
	int i=0, totSub=0, totSubOverN=0, totHost=Settings.n_threshold; 

		for (i=0; i< totHost; i++) 
			totSub += Haloes[i].n_satellites;

		for(i=0; i<Settings.n_haloes; i++)
			if(Haloes[i].host > -1 && Haloes[i].host < totHost && Haloes[i].n_part > Settings.n_min-1) 
				totSubOverN++;

			fprintf(stdout, "\nThere are %d sub Haloes in total.\n", 
				totSub);

			fprintf(stdout, "\nThere are %d sub Haloes in total with more than %d particles.\n", 
				totSubOverN, Settings.n_min);

	Settings.n_subhaloes = totSub;
	Settings.n_subhaloes_nmin = totSubOverN;
	SubHaloes = (struct halo*) calloc(totSub, sizeof(struct halo));
}



void load_subhalo_list()
{
	int i=0, j=0;

		for(i=0; i<Settings.n_haloes; i++)
		{
			if(Haloes[i].host > 0)
			{
				copy_subhalo_properties(j, Haloes[i].id);
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
