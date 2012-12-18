#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "libio/halo_io.h"
#include "libio/read_io.h"

#include "general_variables.h"
#include "general_functions.h"

int main(int argc, char **argv){
	
	// TODO rewrite this

/*	int skip = atoi(argv[1]);
	int n_haloes=atoi(argv[2]);
	int model_1=atoi(argv[3]);
	int model_2=atoi(argv[4]);
	int nBins=atoi(argv[5]);
	double vir=atof(argv[6]);
	double mass=atof(argv[7]);
	char *analysis_dir=argv[8];
	char *halo_files_1=argv[9];
	char *halo_files_2=argv[10];
	char *cc_file=argv[11];
	char *out = argv[12]; 
	char *out_bin = argv[13];

	Urls_internal.analysis_dir=analysis_dir;
	Urls_internal.actual_dir=halo_files_1;

	Settings.model_type=model_1;
	MF.lines_skip=skip;
	MF.n_bins=nBins;
	Cosmo.virial=vir;
	Cosmo.OmegaM=1.;

	int threshold_1=0; int threshold_2=0; int max=0; 
	int gen_mt=0;

	//char *out = merge_strings(analysis_dir,"z_formation_lcdmvde05.dat");
	//char *out_bin = merge_strings(analysis_dir,"z_formation_binned_lcdmvde05.dat");

	double *zf_ratio; int *cc_ids;
	double *z_form_1; double *mass_1; 
	double *z_form_2; double *mass_2; 
	double *z_form_1b; double *mass_1b;
	double *z_form_2b; double *mass_2b;
	z_form_1 = (double*) calloc(n_haloes, sizeof(double));
	z_form_2 = (double*) calloc(n_haloes, sizeof(double));
	mass_1 = (double*) calloc(n_haloes, sizeof(double));
	mass_2 = (double*) calloc(n_haloes, sizeof(double));
	zf_ratio = (double*) calloc(n_haloes, sizeof(double));
	cc_ids = (int*) calloc(n_haloes, sizeof(int));

	initialize_redshift();

	read_and_store_halo_urls();
	
	Urls_internal.halo_file = FC.urls[FC.numFiles-1];

	read_halo_file(); fprintf(stderr, "Initialized halo file 1. \n");
	
	threshold_1=haloes_over_threshold(mass); 
	MF.halos_over_threshold=threshold_1; 
	fprintf(stderr, "Haloes over mass threshold: %d\n", threshold_1);

	initialize_merger_tree();

	if(gen_mt==1) {
generate_merger_tree();
} else {
	read_merger_tree_files();

	Urls_internal.halo_file = FC.urls[FC.numFiles-1];

	read_halo_file();

	MF.halos_over_threshold=threshold_1; 

	print_merger_tree();

	int i=0; int j=0; int k=0;

	for(i=0; i<n_haloes; i++){
	j = i; //cc_ids[i];
z_form_1[i] = haloes[j].z_form;
mass_1[i] = haloes[j].Mvir;
}

	z_form_1b = (double*) calloc(threshold_1, sizeof(double));
	mass_1b = (double*) calloc(threshold_1, sizeof(double));

for(k=0; k<threshold_1; k++){
z_form_1b[k]=haloes[k].z_form;
mass_1b[k]=haloes[k].Mvir;
}

	Urls_internal.actual_dir=halo_files_2;

	Settings.model_type=model_2;

	read_and_store_halo_urls();

	Urls_internal.halo_file = FC.urls[FC.numFiles-1];

	read_halo_file(); fprintf(stderr, "Initialized halo file 2. \n");

	cc_ids=read_cross_correlated_haloes(cc_file, n_haloes, cc_ids);

	max=int_maximum(cc_ids, n_haloes);
        threshold_2=haloes_over_threshold(mass);
	fprintf(stderr, "Haloes over mass threshold: %d, max: %d\n", threshold_2, max);
		
if(max>threshold_2){
	MF.halos_over_threshold=max; 
} else {
	MF.halos_over_threshold=threshold_2; 
}
	initialize_merger_tree();

	if(gen_mt==1) generate_merger_tree();

	read_merger_tree_files();

	Urls_internal.halo_file = FC.urls[FC.numFiles-1];

	read_halo_file();

if(max>threshold_2){
	MF.halos_over_threshold=max; 
} else {
	MF.halos_over_threshold=threshold_2; 
}

	print_merger_tree();

fprintf(stderr, "Done merger tree, printing files.\n");

	FILE *out_z=NULL;  out_z = fopen(out, "w");
	FILE *out_z_bin=NULL; out_z_bin = fopen(out_bin, "w");

fprintf(out_z, "#mass_1(1)   mass(2)   z_form_1(3)   z_form_2(4)   ratio(5)\n");
fprintf(out_z_bin,"#mass_1(1)   mass_2(2)   binned_1(3)   binned_2(4)   mass_bin(5)  ratio(6)\n");

	z_form_2b = (double*) calloc(threshold_2, sizeof(double));
	mass_2b = (double*) calloc(threshold_2, sizeof(double));
for(i=0; i<threshold_2; i++){
mass_2b[i]=haloes[i].Mvir;
z_form_2b[i]=haloes[i].z_form;
}

for(i=0; i<n_haloes; i++){
//fprintf(stderr, "i:%d, n:%d \n", i, n_haloes);
	j = cc_ids[i];
z_form_2[i]=haloes[j].z_form;
mass_2[i]=haloes[j].Mvir;
zf_ratio[i] = z_form_2[i]/z_form_1[i];
//fprintf(stderr, "i:%d) j:%d) z_f_1:%lf, z_f_2:%lf, r:%lf \n", i, j, z_form_1[i], z_form_2[i], zf_ratio[i]);
fprintf(out_z,"%e  %e  %lf  %lf  %lf\n", mass_1[i], mass_2[i], z_form_1[i], z_form_2[i], zf_ratio[i]);
}

double *z_binned_ratio; double *mass_bins;
double *mass_bins_1; double *mass_bins_2; 
double *err_1; double *err_2;
double *z_binned_1; double *z_binned_2;

z_binned_ratio = (double*) calloc(nBins,sizeof(double));
z_binned_1 = (double*) calloc(nBins,sizeof(double));
z_binned_2 = (double*) calloc(nBins,sizeof(double));

mass_bins = (double*) calloc(nBins,sizeof(double));
mass_bins_1 = (double*) calloc(nBins,sizeof(double));
mass_bins_2 = (double*) calloc(nBins,sizeof(double));

err_1 = (double*) calloc(nBins,sizeof(double));
err_2 = (double*) calloc(nBins,sizeof(double));

double mMax=mass_1[0]; double mMin=mass_1[n_haloes-1];
double mMax1=mass_1b[0]; double mMin1=mass_1b[threshold_1-1];
double mMax2=mass_2b[0]; double mMin2=mass_2b[threshold_2-1];

mass_bins=stepper(mMin,mMax,nBins);
mass_bins_1=stepper(mMin1,mMax1,nBins);
mass_bins_2=stepper(mMin2,mMax2,nBins);

z_binned_ratio=average_bin(mass_1,zf_ratio,err_1,mass_bins,z_binned_ratio,nBins,n_haloes);
z_binned_1=average_bin(mass_1b,z_form_1b,err_1,mass_bins_1,z_binned_1,nBins,threshold_1);
z_binned_2=average_bin(mass_2b,z_form_2b,err_2,mass_bins_2,z_binned_2,nBins,threshold_2);

fprintf(stderr, "Printing files... \n");

	for(i=0; i<nBins-1; i++){
//fprintf(stderr,"%e  %lf  %lf\n", mass_bins[i], z_binned_1[i], z_binned_2[i]);
//fprintf(stderr,"%e  %lf\n", mass_bins[i], z_binned_ratio[i]);
fprintf(out_z_bin,"%e  %e  %lf  %lf  %e  %lf\n",mass_bins_1[i],mass_bins_2[i],z_binned_1[i],z_binned_2[i],mass_bins[i],z_binned_ratio[i]);
}
i=nBins-1;
fprintf(out_z_bin,"%e  %e  %lf  %lf  %e  %lf\n",mass_bins_1[i],mass_bins_2[i],z_binned_1[i],z_binned_2[i],mass_bins[i],z_binned_ratio[i-1]);

fclose(out_z_bin);
fclose(out_z);
}
*/
return 0;
}
