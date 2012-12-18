#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "libhalo/nfw.h"
#include "libio/read_io.h"
#include "libmath/statistics.h"
#include "libmath/mathtools.h"
#include "general_variables.h"
#include "general_functions.h"

void initialize_internal_variables(char **argv){
	fprintf(stderr, "Initializing internal variables...\n");

	Urls_internal.analysis_dir = argv[1];
	Urls_internal.a_outputs = argv[2];
	Urls_internal.pk_file = argv[3];
	Urls_internal.halo_file = argv[4];
	Urls_internal.profiles_file = argv[5];
	Urls_internal.pk_root = argv[6];
	Urls_internal.snaps_dir = argv[7];
	Urls_internal.halo_dir = argv[8];

		int URLS=8;

	Settings.box_size = atof(argv[URLS+1]); 
	Settings.nP_1D = atoi(argv[URLS+2]); 
	Settings.n_bins=atoi(argv[URLS+3]);
	Settings.pk_skip=atoi(argv[URLS+4]);
	Settings.halo_skip=atoi(argv[URLS+5]);
	Settings.fit = atoi(argv[URLS+6]);
	Settings.zStart = atof(argv[URLS+7]);
	Settings.thMass=atof(argv[URLS+8]);
	Settings.Mmin = atof(argv[URLS+9]);
	Settings.Mmax = atof(argv[URLS+10]);
	Settings.Rmin=atof(argv[URLS+11]);
	Settings.Rmax=atof(argv[URLS+12]);
	Settings.r_bins=atoi(argv[URLS+13]);
	Settings.thNum=(int) atof(argv[URLS+14]);
	Settings.n_min=atof(argv[URLS+15]);

		int SET=15;

	Cosmo.h=atof(argv[SET+URLS+1]);
	Cosmo.s8=atof(argv[SET+URLS+2]);
	Cosmo.OmegaM=atof(argv[SET+URLS+3]);
	Cosmo.OmegaL=atof(argv[SET+URLS+4]);
	Cosmo.delta_c=atof(argv[SET+URLS+5]);
	Cosmo.spin=atof(argv[SET+URLS+6]);
	Cosmo.virial=atof(argv[SET+URLS+7]);
	GF.scale_k=atof(argv[SET+URLS+8]);
	ND.zMax=atof(argv[SET+URLS+9]);

		int COSMO=9;

	// Extra analysis dependent parameters - set as the last one
	Urls_internal.output_prefix = argv[SET+COSMO+URLS+1];
	Urls_internal.pk_list = argv[SET+COSMO+URLS+2];
	Urls_internal.halo_list = argv[SET+COSMO+URLS+3];
	Urls_internal.profile_list = argv[SET+COSMO+URLS+4];
	Urls_internal.subhalo_list = argv[SET+COSMO+URLS+5];
	FC.numFiles = atoi(argv[SET+COSMO+URLS+6]);

#ifdef PRINT_INFO
	int kk=0; 
		fprintf(stderr, "settings: %d\n", FC.numFiles);

		for(kk=1; kk<URLS+SET+COSMO+7; kk++) 
			fprintf(stderr, "argv[%d]: %s \n", kk, argv[kk]);
#endif

	// Setting extra useful variables
	Cosmo.H_0=Cosmo.h*100;
	AMF.Mmin=Settings.Mmin;
	AMF.Mmax=Settings.Mmax;
	Settings.Gn=6.672e-8;
}

void print_counter(int freq){
		if(Settings.tick==freq) {
	fprintf(stderr,"."); 
	Settings.tick=0;
		}
	else {
	Settings.tick++;
	}
}

void normalize_to_one(char * url_in, char * url_out){

		FILE *file_in = fopen(url_in, "r");
		FILE *file_out = fopen(url_out, "w");
		int inv=0, k=0, size = get_lines(file_in, url_in);
		double norm=1., *x, *y;
		char line[256];

		x   = (double*) calloc(size,sizeof(double));
		y   = (double*) calloc(size,sizeof(double));

		for(k=0; k<size; k++){
		fgets(line,256,file_in);
		sscanf(line,"%lf %lf", &x[k], &y[k]);
		}

		if (inv==1) {
			norm = y[size-1];
		fprintf(stderr, "norm: %lf \n", y[size-1]);
		} else {
			norm = y[0];
		fprintf(stderr, "norm: %lf \n", y[0]);
		}

		for(k=0; k<size; k++)
			fprintf(file_out, "%e  %e\n", x[k], y[k]/norm);

	fclose(file_in);
	fclose(file_out);
}

char* merge_strings(char* str_1, char* str_2){

	int dim = strlen(str_1) + strlen(str_2) + 2;
	char *str_3;
	str_3 = (char *) calloc(dim, sizeof(char));

	sprintf(str_3, "%s%s", str_1, str_2);	

return str_3;
};

void compute_files_ratio(char *url1, char *url2, char *url_out, int bins, int skip, int logStep){

	FILE *f1 = fopen(url1, "r");
	FILE *f2 = fopen(url2, "r");
	FILE *output = fopen(url_out, "w");
	char dummyline[2048];
	int n1 = get_lines(f1, url1) - skip, n2 = get_lines(f2, url2) - skip, i=0, j=0;
	
	double *vec_x1, *vec_y1, *vec_x2, *vec_y2;
	double *arr_x1, *arr_y1, *arr_x2, *arr_y2;
	double max, max1, max2, min, min1, min2;
	double x, ratio, y1=0, y2=0, eps = 0.001, *step;

	vec_x1 = (double *) malloc(n1 * sizeof(double));
    	vec_y1 = (double *) malloc(n1 * sizeof(double));
	vec_x2 = (double *) malloc(n2 * sizeof(double));
	vec_y2 = (double *) malloc(n2 * sizeof(double));
	arr_x1 = (double *) malloc(n1 * sizeof(double));
    	arr_y1 = (double *) malloc(n1 * sizeof(double));
	arr_x2 = (double *) malloc(n2 * sizeof(double));
	arr_y2 = (double *) malloc(n2 * sizeof(double));

	for(i=0;i<n1+skip;i++){
			fgets(dummyline, 2048, f1);
	if(i>skip-eps)
			{
				j=i-skip;
			sscanf(dummyline, "%lf %lf", &vec_x1[j], &vec_y1[j]);
		}
	}


	i=0;

			while(!feof(f2)){
				fgets(dummyline, 2048, f2);
			if(i>skip-eps){
				j=i-skip;
					sscanf(dummyline, "%lf %lf", &vec_x2[j], &vec_y2[j]);
				}
					i++;
						}

		fclose(f1);
		fclose(f2);

	/* Check if the arrays have to be inverted - GSL consistency */
	if(vec_x1[n1-1] < vec_x1[0]) {
			fprintf(stdout, "\nInvert array 1 \n");

		for(i=0; i<n1; i++)
			{
	arr_x1[i]=vec_x1[n1-i-1];
	arr_y1[i]=vec_y1[n1-i-1];
		}

			} else {
	fprintf(stdout, "\nDo not invert array 1\n");

		for(i=0; i<n1; i++)
			{
	arr_x1[i]=vec_x1[i];
	arr_y1[i]=vec_y1[i];
		}

	}

	if(vec_x2[n2-1] < vec_x2[0]) {
	fprintf(stdout, "\nInvert array 2\n");

		for(i=0; i<n2; i++)
				{
			arr_x2[i]=vec_x2[n2-i-1];
			arr_y2[i]=vec_y2[n2-i-1];
		}
	} else {
	fprintf(stdout, "\nDo not invert array 2 \n");
	
		for(i=0; i<n2; i++)
				{
	arr_x2[i]=vec_x2[i];
	arr_y2[i]=vec_y2[i];
		}
	}

	max1=arr_x1[n1-1];
	max2=arr_x2[n2-1];
	
		if(max1 > max2)
			{		
		max = max2;
			} else {
		max = max1;
		}

	min1=arr_x1[0];
	min2=arr_x2[0];

		if(min1 < min2) 
		{
		min = min2;
			} else {
		min = min1;
		}

	if(logStep==1)
			{
		step = log_stepper(min,max,bins);
			} else {
		step = lin_stepper(min,max,bins);
	}

		fprintf(stdout, "\nMAX: %e Min:%e\n", max, min);

	for(i=0; i<bins; i++)
		{
			x = step[i];
				y1 = get_interpolated_value(arr_x1, arr_y1, n1, x);
				y2 = get_interpolated_value(arr_x2, arr_y2, n2, x);
			ratio = y2/y1;
		fprintf(output, "%e %lf \n", x, ratio);
	}

	free(vec_x1);
	free(vec_x2);
	free(vec_y1);
	free(vec_y2);
	free(arr_x1);
	free(arr_x2);
	free(arr_y1);
	free(arr_y2);

fclose(output);
}
