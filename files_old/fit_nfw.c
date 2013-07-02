#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "libio/halo_io.h"
#include "libgen/vars.h"
#include "libgen/mathtools.h"
#include "libgen/commons.h"
#include "libhalo/nfw.h"
#include "libhalo/halo_properties.h"
#include "libstat/statistics.h"

void calculate_and_sort_chi2(){
	int massCut = Settings.haloes_over_threshold;
	int nBins = NFW.bins;
	double vir = Cosmo.virial; 
	double *a_x; double *a_y; double v;
	double cc, mm, rr, rho0, rs, r_back;
	int n=0; int i=0; int k=0;
	int *vir_haloes; int n_vir=0;

	for(k=0; k<massCut; k++){
	v=haloes[k].th_vir;
	if(sqrt(v*v)<vir) n_vir++;
	}

fprintf(stderr, "Total of %d virialized haloes over %d. \n", n_vir, massCut);
	vir_haloes=(int*) calloc(n_vir, sizeof(int));

	for(k=0; k<massCut; k++){
	v=haloes[k].th_vir;
	if(sqrt(v*v)<vir) {
	vir_haloes[i]=k; i++;
	}
	}
	
//for(k=0; k<haloes[0].n_bins; k++)fprintf(stderr, "r:%lf, rho:%lf \n", haloes[0].radius[k], haloes[0].rho[k]);

Chi2.chi2s = (double*) calloc(n_vir, sizeof(double));
fit_and_store_nfw_parameters_from_list(n_vir, vir_haloes, massCut);

i=0;
for(k=0; k<n_vir; k++){
i=vir_haloes[k];
Chi2.chi2s[k] = haloes[i].chi_nfw;
//fprintf(stderr, "%d) chi2:%lf \n", k, Chi2.chi2s[k]);
}
i=0; 
k=0;

Chi2.chi2s = shellsort(Chi2.chi2s, n_vir);
Chi2.bins = nBins;
Chi2.binned_chi = (double*) calloc(Chi2.bins, sizeof(double));
Chi2.outcomes = (int*) calloc(Chi2.bins, sizeof(int));
Chi2.binned_chi = lin_stepper(Chi2.chi2s[0], Chi2.chi2s[massCut-1], nBins);

FILE *nfwfit = fopen(Urls_internal.output_prefix, "w");
//FILE *nfwfit = fopen("nfw_chi_2_distribution-vdez0_1.dat", "w");
//for(n=0; n<massCut; n++) fprintf(stderr, "%d, chi: %lf \n", n, Chi2.chi2s[n]);
fprintf(stderr, "\nBinning chi2...");
i=0; n=0; k=0;

Chi2.outcomes = lin_bin(Chi2.chi2s, Chi2.binned_chi, Chi2.bins, n_vir, Chi2.outcomes);
fprintf(stderr, "Binned.\n");

for(n=0; n<Chi2.bins; n++) {
fprintf(nfwfit, "%lf %d %lf\n", Chi2.binned_chi[n], Chi2.outcomes[Chi2.bins-n-1], (double) Chi2.outcomes[Chi2.bins-n-1]/(double)massCut);
//fprintf(stderr, "%lf %d %lf\n", 0, 0, 0); //Chi2.binned_chi[n], 0, 0); //Chi2.outcomes[Chi2.bins-n-1], (double) Chi2.outcomes[Chi2.bins-n-1]/(double)massCut);
}
fprintf(stderr, "Done.\n");
fclose(nfwfit);
}

int main(int argc, char **argv){
fprintf(stderr, "Computing the Chi^2 distribution for the given halo sample.\n");

	initialize_internal_variables(argv);

	Cosmo.err=0.3;  Cosmo.OmegaM=1.;

	read_halo_file();

	read_profiles_file();

	calculate_and_sort_chi2();

return 0;
}
