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

#ifdef WITH_MPI
#include <mpi.h>
#include "libparallel/general.h"
#endif



void initialize_internal_variables(char **argv){
	fprintf(stdout, "Initializing internal variables...\n");
		
	int count=1;

	// Basic working file urls
	Urls.a_outputs = argv[count++];
	Urls.halo_file = argv[count++];
	Urls.profiles_file = argv[count++];
	Urls.pk_file = argv[count++];

	// Basic simulation and analysis settings
	Settings.box_size = atof(argv[count++]); 
	Settings.n_part_1D = atoi(argv[count++]); 
	Settings.n_bins = atoi(argv[count++]);
	Settings.n_bins_th = atoi(argv[count++]);
	Settings.r_bins = atoi(argv[count++]);
	Settings.pk_skip = atoi(argv[count++]);
	Settings.halo_skip = atoi(argv[count++]);
	Settings.cat_number = atoi(argv[count++]);
	Settings.fit = atoi(argv[count++]);
	Settings.zStart = atof(argv[count++]);
	Settings.mass_min = atof(argv[count++]);
	Settings.Mmin = atof(argv[count++]);
	Settings.Mmax = atof(argv[count++]);
	Settings.Rmin = atof(argv[count++]);
	Settings.Rmax = atof(argv[count++]);
	Settings.n_min = atof(argv[count++]);
	Settings.use_n_min = atof(argv[count++]);
	Settings.n_haloes_to_use = atoi(argv[count++]); 
	
	// Cosmological parameters
	Cosmo.h = atof(argv[count++]);
	Cosmo.sigma8 = atof(argv[count++]);
	Cosmo.OmegaM = atof(argv[count++]);
	Cosmo.OmegaL = atof(argv[count++]);
	Cosmo.delta_c = atof(argv[count++]);
	Cosmo.spin = atof(argv[count++]);
	Cosmo.virial = atof(argv[count++]);
	GrowthFac.scale_k = atof(argv[count++]);
	NumDen.zMax = atof(argv[count++]);

	// Other urls related parameters
	Urls.output_prefix = argv[count++];
	Urls.halo_list = argv[count++];
	Urls.profile_list = argv[count++];
	Urls.subhalo_list = argv[count++];
	Urls.pk_list = argv[count++];
	Urls.nCatalogueFiles = atoi(argv[count++]);

//	fprintf(stdout, "Prefix=%s\n", Urls.output_prefix);

#ifdef PRINT_INFO
	int kk=0;
 
	for(kk=1; kk<count; kk++) 
		fprintf(stdout, "argv[%d]: %s \n", kk, argv[kk]);
#endif

	default_init();

#ifdef _OPENMP
	omp_set_num_threads(4);
#endif
}



void default_init()
{
	// Setting extra useful variables
	Cosmo.H_0=Cosmo.h*100;
	Cosmo.G=6.672e-8;
	ThMassFunc.Mmin=Settings.Mmin;
	ThMassFunc.Mmax=Settings.Mmax;
	MassFunc.bins  = Settings.n_bins; 
	ThMassFunc.bins = Settings.n_bins_th;
	Settings.use_cat=Urls.nCatalogueFiles-Settings.cat_number;
	
	// Init some commonly used structures to default values
	GrowthFac.z = (double *) calloc(1, sizeof(double));
	GrowthFac.a = (double *) calloc(1, sizeof(double));
	GrowthFac.gf= (double *) calloc(1, sizeof(double));
	Pks = (struct power_spectrum *) calloc(1, sizeof(struct power_spectrum));
}



void print_counter(int freq){

	if(Settings.tick==freq) 
	{
		fprintf(stdout,"."); 
		Settings.tick=0;
	} else {
		Settings.tick++;
	}
}



void normalize_to_one(char * url_in, char * url_out)
{
		int inv=0, k=0, size=0;
		double norm=1.; 
		double *x=NULL, *y=NULL;
		char line[256];

		FILE *file_in = fopen(url_in, "r");
		FILE *file_out = fopen(url_out, "w");

		size = get_lines(file_in, url_in);

		x   = (double*) calloc(size,sizeof(double));
		y   = (double*) calloc(size,sizeof(double));

		for(k=0; k<size; k++)
		{
			fgets(line,256,file_in);
			sscanf(line,"%lf %lf", &x[k], &y[k]);
		}

		if (inv==1) 
		{
			norm = y[size-1];
			fprintf(stdout, "norm: %lf \n", y[size-1]);
		} else {
			norm = y[0];
			fprintf(stdout, "norm: %lf \n", y[0]);
		}

		for(k=0; k<size; k++)
			fprintf(file_out, "%e  %e\n", x[k], y[k]/norm);

	fclose(file_in);
	fclose(file_out);
}



char* merge_strings(char* str_1, char* str_2)
{
	int dim=0; 
	char *str_3=NULL;

		dim = strlen(str_1) + strlen(str_2) + 2;
		str_3 = (char *) calloc(dim, sizeof(char));

	sprintf(str_3, "%s%s", str_1, str_2);	

	return str_3;
}



void compute_files_ratio(char *url1, char *url2, char *url_out, int bins, int skip, int logStep)
{
	int n1=0, n2=0, i=0, j=0;
	
	double *vec_x1, *vec_y1, *vec_x2, *vec_y2;
	double *arr_x1, *arr_y1, *arr_x2, *arr_y2;
	double x=0, ratio=0, y1=0, y2=0, eps=0.001;
	double max, max1, max2, min, min1, min2;
	double *step=NULL;
	char dummyline[2048];

	FILE *f1 = fopen(url1, "r");
	FILE *f2 = fopen(url2, "r");
	FILE *output = fopen(url_out, "w");

	n1 = get_lines(f1, url1) - skip; 
	n2 = get_lines(f2, url2) - skip; 

	vec_x1 = (double *) malloc(n1 * sizeof(double));
    	vec_y1 = (double *) malloc(n1 * sizeof(double));
	vec_x2 = (double *) malloc(n2 * sizeof(double));
	vec_y2 = (double *) malloc(n2 * sizeof(double));
	arr_x1 = (double *) malloc(n1 * sizeof(double));
    	arr_y1 = (double *) malloc(n1 * sizeof(double));
	arr_x2 = (double *) malloc(n2 * sizeof(double));
	arr_y2 = (double *) malloc(n2 * sizeof(double));

	for(i=0;i<n1+skip;i++)
	{
			fgets(dummyline, 2048, f1);

		if(i>skip-eps)
		{
			j=i-skip;
			sscanf(dummyline, "%lf %lf", &vec_x1[j], &vec_y1[j]);
		}
	}

		i=0;

			while(!feof(f2))
			{
				fgets(dummyline, 2048, f2);

				if(i>skip-eps)
				{
					j=i-skip;
					sscanf(dummyline, "%lf %lf", &vec_x2[j], &vec_y2[j]);
				}
					i++;
			}

		fclose(f1);
		fclose(f2);

	/* Check if the arrays have to be inverted - GSL consistency */
	if(vec_x1[n1-1] < vec_x1[0]) 
	{
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

	if(vec_x2[n2-1] < vec_x2[0]) 
	{
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
