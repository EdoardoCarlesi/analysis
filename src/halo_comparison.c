#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "libio/read_io.h"
#include "libio/halo_io.h"
#include "libhalo/nfw.h"
#include "libmath/statistics.h"

#include "general_variables.h"
#include "general_functions.h"

int main(int argc, char **argv){
// TODO review this !!

	/* 
	char *halo_file_1 = argv[1];
	char *halo_file_2 = argv[2];
	char *halo_cc_list = argv[3];
	char *profiles_file_1 = argv[4];
	char *profiles_file_2 = argv[5];

	int skip = atoi(argv[6]);
	int n_halos_comp = atoi(argv[7]);
	int nBins = atoi(argv[8]);
	double box_size = atof(argv[9]);
	int compare_profiles = atoi(argv[10]);
	//char *halo_cc_ratio = merge_strings(Urls_internal.analysis_dir,"/output/halo_cc_ratio.dat");
	//char *halo_cc_ratio_bin = merge_strings(Urls_internal.analysis_dir,"output/halo_cc_binned_ratio.dat");
	char *halo_cc_ratio = argv[11];
	char *halo_cc_ratio_bin = argv[12];

	MF.n_bins=nBins;
	MF.lines_skip=skip; 
	Settings.box_size=box_size; 
	Cosmo.virial=10.5;
	Cosmo.OmegaM=1.;
	Chi2.chi2s = (double*) calloc(n_halos_comp, sizeof(double));
	double err_1=1.1; double err_2=0.75;
	initialize_internal_variables(argv);
	int *cc_ids=calloc(n_halos_comp,sizeof(int));

	Urls_internal.halo_file = halo_file_1;
	Urls_internal.profiles_file = profiles_file_1;
	Cosmo.err = err_1;
	read_halo_file();
	fprintf(stderr, "Initialized halo file 1. \n");

	cc_ids=read_cross_correlated_haloes(halo_cc_list, n_halos_comp, cc_ids);
	int max = int_maximum(cc_ids, n_halos_comp);
fprintf(stderr, "Max number of haloes over selected virial threshold: %d\n", max);
	MF.halos_over_threshold=max; 
	
	if (compare_profiles==1) {
	read_profiles_file();
	fprintf(stderr, "Initialized profiles file 1. \n");
	Chi2.chi2s = (double*) calloc(n_halos_comp, sizeof(double));
	fit_and_store_nfw_parameters_from_list(n_halos_comp, cc_ids, max);
}
	double *masses_1;
	double *radii_1;
	double *angmom_1;
	double *shape_1;
	double *triax_1;
	double *r2_1;
	double *chi_1;
	double *rs_1;
	double *rho0_1;

	masses_1 = (double*) calloc(n_halos_comp, sizeof(double));
	radii_1 = (double*) calloc(n_halos_comp, sizeof(double));
	angmom_1 = (double* ) calloc(n_halos_comp, sizeof(double));
	triax_1 = (double* ) calloc(n_halos_comp, sizeof(double));
	shape_1 = (double* ) calloc(n_halos_comp, sizeof(double));
	r2_1 = (double* ) calloc(n_halos_comp, sizeof(double));
	chi_1 = (double* ) calloc(n_halos_comp, sizeof(double));
	rs_1 = (double* ) calloc(n_halos_comp, sizeof(double));
	rho0_1 = (double* ) calloc(n_halos_comp, sizeof(double));

	int j=0; int kk=0;
	for(j=0; j<n_halos_comp; j++) {
	kk = cc_ids[j];
	masses_1[j] = haloes[kk].Mvir;
	//fprintf(stderr, "j:%d) kk:%d) masses:%e\n", j, kk, masses_1[j]);
	radii_1[j] = haloes[kk].Rvir;
//	angmom_1[j] = haloes[kk].AngMom;
	angmom_1[j] = haloes[kk].lambda;
	triax_1[j] = haloes[kk].triax;
	shape_1[j] = haloes[kk].c_a;
	r2_1[j] = haloes[kk].r2;
	if (compare_profiles==1) {
	chi_1[j] = haloes[kk].chi_nfw;
	rs_1[j] = haloes[kk].rs_nfw;
	rho0_1[j] = haloes[kk].rho0_nfw;
}}

	Urls_internal.halo_file = halo_file_2;
	Urls_internal.profiles_file = profiles_file_2;

	Cosmo.err = err_2;
	read_halo_file();
	fprintf(stderr, "Initialized halo file 2. \n");
	if (compare_profiles==1) {
	read_profiles_file();
	fprintf(stderr, "Initialized profiles file 2. \n");
	fit_and_store_nfw_parameters(n_halos_comp);
}
	double *masses_2;
	double *rho0_2;
	double *radii_2;
	double *angmom_2;
	double *triax_2;
	double *shape_2;
	double *r2_2;
	double *chi_2;
	double *rs_2;

	fprintf(stderr, "Initialized halo file 1. \n");
	masses_2 = (double *) calloc(n_halos_comp, sizeof(double));
	radii_2 = (double *) calloc(n_halos_comp, sizeof(double));
	angmom_2 = (double *) calloc(n_halos_comp, sizeof(double));
	r2_2 = (double *) calloc(n_halos_comp, sizeof(double));
	chi_2 = (double *) calloc(n_halos_comp, sizeof(double));
	rs_2 = (double *) calloc(n_halos_comp, sizeof(double));
	rho0_2 = (double *) calloc(n_halos_comp, sizeof(double));
	shape_2 = (double *) calloc(n_halos_comp, sizeof(double));
	triax_2 = (double *) calloc(n_halos_comp, sizeof(double));

	fprintf(stderr, "Initialized halo file 1. \n");
	for(j=0; j<n_halos_comp; j++) {
	kk = j; 
	//fprintf(stderr, "j:%d) kk:%d)\n", j, kk);
	masses_2[j] = haloes[kk].Mvir;
	radii_2[j] = haloes[kk].Rvir;
//	angmom_2[j] = haloes[kk].AngMom;
	angmom_2[j] = haloes[kk].lambda;
	r2_2[j] = haloes[kk].r2;
	shape_2[j] = haloes[kk].c_a;
	triax_2[j] = haloes[kk].triax;
	if (compare_profiles==1) {
	chi_2[j] = haloes[kk].chi_nfw;
	rs_2[j] = haloes[kk].rs_nfw;
	rho0_2[j] = haloes[kk].rho0_nfw;
}
}

	int i,k;
	fprintf(stderr, "Initialized halo file 1. \n");
	
	double *m_ratio;
	double *r_ratio;
	double *am_ratio;
	double *c_ratio;
	double *chi_ratio;
	double *rs_ratio;
	double *triax_ratio;
	double *shape_ratio;
	double *rho0_ratio;

	m_ratio=(double *) calloc(n_halos_comp, sizeof(double));	
	r_ratio=(double *) calloc(n_halos_comp, sizeof(double));	
	am_ratio=(double *) calloc(n_halos_comp, sizeof(double));	
	shape_ratio=(double *) calloc(n_halos_comp, sizeof(double));	
	triax_ratio=(double *) calloc(n_halos_comp, sizeof(double));	
	c_ratio=(double *) calloc(n_halos_comp, sizeof(double));	
	chi_ratio=(double *) calloc(n_halos_comp, sizeof(double));	
	rs_ratio=(double *) calloc(n_halos_comp, sizeof(double));	
	rho0_ratio=(double *) calloc(n_halos_comp, sizeof(double));	

	FILE *hccr = fopen(halo_cc_ratio, "w");
	FILE *hccr_bin = fopen(halo_cc_ratio_bin, "w");

// Print header file	
fprintf(hccr, "#massLCDM(1)  m_ratio(2)  r_ratio(3)  lambda_ratio(4)  c_ratio(5)  chi_ratio(6)  rs_ratio(7)  rho0_ratio(8)  shape_ratio(9)  triax_ratio(10)\n");

// masses_1 = VDE, haloes = LCDM

	for(j=0; j<n_halos_comp; j++){
	m_ratio[j] = masses_1[j]/haloes[j].Mvir;if(m_ratio[j] > 100) m_ratio[j] = 1;
	r_ratio[j] = radii_1[j]/haloes[j].Rvir; if(r_ratio[j] > 100) r_ratio[j] = 1;
	//am_ratio[j] = angmom_1[j]/haloes[j].AngMom; if(am_ratio[j] > 100) am_ratio[j] = 1;
	am_ratio[j] = angmom_1[j]/haloes[j].lambda; if(am_ratio[j] > 100) am_ratio[j] = 1;

	if (compare_profiles==1) {
	c_ratio[j] = r_ratio[j]*(haloes[j].r2/r2_1[j]);if(c_ratio[j] > 100) c_ratio[j] = 1;
	chi_ratio[j] = chi_1[j]/haloes[j].chi_nfw; if(chi_ratio[j] > 9.9) chi_ratio[j] = 1;
	if(chi_ratio[j] < 0.11) chi_ratio[j] = 1;
	rs_ratio[j] = rs_1[j]/haloes[j].rs_nfw;if(rs_ratio[j] > 100) rs_ratio[j] = 1;
	rho0_ratio[j] = rho0_1[j]/haloes[j].rho0_nfw;if(rho0_ratio[j] > 100) rho0_ratio[j] = 1;
} else {

	c_ratio[j] = r_ratio[j]*(haloes[j].r2/r2_1[j]);if(c_ratio[j] > 100) c_ratio[j] = 1;
}
	shape_ratio[j] = shape_1[j]/haloes[j].c_a;if(shape_ratio[j] > 100) shape_ratio[j] = 1;
	triax_ratio[j] = triax_1[j]/haloes[j].triax;if(triax_ratio[j] > 100) triax_ratio[j] = 1;
	//fprintf(stderr, "R1: %lf, R2: %lf, R3: %lf \n", m_ratio[j], r_ratio[j], am_ratio[j]);
	//fprintf(stderr, "chi1: %lf, chi2: %lf, rs1: %lf chi_2: %lf\n", chi_1[j], haloes[j].chi_nfw, rs_1[j], chi_2[j] );
	fprintf(hccr, "%e   %lf   %lf   %lf   %lf   %lf   %lf  %lf  %lf  %lf\n", masses_2[j]/haloes[j].Mvir, m_ratio[j], r_ratio[j], am_ratio[j], c_ratio[j], chi_ratio[j], rs_ratio[j], rho0_ratio[j], shape_ratio[j], triax_ratio[j]);
	}

	fprintf(stderr, "Binnig the ratios.\n");

	//double mMax=1.e14; //masses_1[0]; 
	double mMax=masses_2[0]; 
	double mMin=masses_2[n_halos_comp-1];
	//double mMax=5.e+15; double mMin=1.e+14;
	//double mMax=masses_1[0]; double mMin=masses_1[max-1];
	double *steps; 
	steps = (double *) calloc(nBins, sizeof(double));
	steps = stepper(mMin, mMax, nBins);

	double *m_ratio_bin;
	double *r_ratio_bin;
	double *am_ratio_bin;
	double *void_err;
	double *c_ratio_bin;
	double *chi_ratio_bin;
	double *rs_ratio_bin;
	double *rho0_ratio_bin;
	double *shape_ratio_bin;
	double *triax_ratio_bin;

	void_err=(double *) calloc(nBins, sizeof(double));	

	m_ratio_bin=(double *) calloc(nBins, sizeof(double));	
	m_ratio_bin=average_bin(masses_2, m_ratio, void_err, steps, m_ratio_bin, nBins, n_halos_comp);

	r_ratio_bin=(double *) calloc(nBins, sizeof(double));	
	r_ratio_bin=average_bin(masses_2, r_ratio, void_err, steps, r_ratio_bin, nBins, n_halos_comp);

	am_ratio_bin=(double *) calloc(nBins, sizeof(double));	
	am_ratio_bin=average_bin(masses_2, am_ratio, void_err, steps, am_ratio_bin, nBins, n_halos_comp);

	c_ratio_bin=(double *) calloc(nBins, sizeof(double));	
	c_ratio_bin=average_bin(masses_2, c_ratio, void_err, steps, c_ratio_bin, nBins, n_halos_comp);

	chi_ratio_bin=(double *) calloc(nBins, sizeof(double));

        rs_ratio_bin=(double *) calloc(nBins, sizeof(double));	

	rho0_ratio_bin=(double *) calloc(nBins, sizeof(double));	

	if (compare_profiles==1) {
	chi_ratio_bin=average_bin(masses_2, chi_ratio, void_err, steps, chi_ratio_bin, nBins, n_halos_comp);
	rs_ratio_bin=average_bin(masses_2, rs_ratio, void_err, steps, rs_ratio_bin, nBins, n_halos_comp);
	rho0_ratio_bin=average_bin(masses_2, rho0_ratio, void_err, steps, rho0_ratio_bin, nBins, n_halos_comp); }
	
	shape_ratio_bin=(double *) calloc(nBins, sizeof(double));	
	shape_ratio_bin=average_bin(masses_2, shape_ratio, void_err, steps, shape_ratio_bin, nBins, n_halos_comp);
	
	triax_ratio_bin=(double *) calloc(nBins, sizeof(double));	
	triax_ratio_bin=average_bin(masses_2, triax_ratio, void_err, steps, triax_ratio_bin, nBins, n_halos_comp);

	//for(i=0; i<nBins; i++) fprintf(stderr, "M:%e   r:%lf \n", steps[i], r_ratio_bin[i]);
	
fprintf(hccr_bin, "#mass(1)  m_ratio(2)  r_ratio(3)  lambda_ratio(4)  c_ratio(5)  chi_ratio(6)  rs_ratio(7)  rho0_ratio(8)  shape_ratio(9)  triax_ratio(10)\n");
	
	for(i=1; i<nBins; i++){ 
	fprintf(hccr_bin, "%e  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n", 
steps[i-1], m_ratio_bin[i], r_ratio_bin[i], am_ratio_bin[i], 
c_ratio_bin[i], chi_ratio_bin[i], rs_ratio_bin[i], rho0_ratio_bin[i],
shape_ratio_bin[i], triax_ratio_bin[i]);
	}
	// Save the last bin
	fprintf(hccr_bin, "%e  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n", 
steps[nBins-1], m_ratio_bin[nBins-1], r_ratio_bin[nBins-1], am_ratio_bin[nBins-1], 
c_ratio_bin[nBins-1], chi_ratio_bin[nBins-1], rs_ratio_bin[nBins-1], rho0_ratio_bin[nBins-1],
shape_ratio_bin[i], triax_ratio_bin[i]);
	fprintf(stderr, "Done.\n");

	free(m_ratio);
	free(am_ratio);
	free(r_ratio);
	free(c_ratio);
	free(chi_ratio);
	free(rs_ratio);
	free(rho0_ratio);
	free(shape_ratio);
	free(triax_ratio);
	fclose(hccr);

*/
return 0;
}
