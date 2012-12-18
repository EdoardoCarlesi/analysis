#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "read_io.h"
#include "../general_variables.h"
#include "../libcosmo/cosmological_relations.h"


int get_lines(FILE *f, char *url)
{
	int n=-1;
	char line[2048];
	
		if(f==NULL) 
			fprintf(stderr,"\nError in get_lines(). File %s could not be opened.\n", url);

		while(!feof(f))
		{
			fgets(line,2048,f);
			n++;
		}

		rewind(f);

	return n;
}



void read_redshift_file()
{
	int size=0, i=0;
	double aa;
	char dummyline[50];
	FILE *redshifts_file = NULL;

		redshifts_file = fopen(Urls_internal.a_outputs, "r");

		if(redshifts_file == NULL) 
			fprintf(stderr, "\nError opening redshift file %s\n", Urls_internal.a_outputs);
		else 
			fprintf(stderr, "Redshift file opened correctly:%s\n", Urls_internal.a_outputs);

			size = get_lines(redshifts_file, Urls_internal.a_outputs);

			Settings.n_pk_files = size;
			GF.npts = size;
			GF.z = (double *) calloc(size, sizeof(double));
			GF.a = (double *) calloc(size, sizeof(double));
			GF.gf_z = (double *) calloc(size, sizeof(double));
			GF.gf_over_a_z = (double *) calloc(size, sizeof(double));


				for(i=0; i<size; i++) 
				{
					fgets(dummyline, 50, redshifts_file);
					sscanf(dummyline, "%lf", &aa);
					GF.a[i] = aa;
					GF.z[i] = 1./aa - 1.;
#ifdef PRINT_INFO
			fprintf(stderr, "%d) a:%lf z:%lf\n", i, aa, GF.z[i]);
#endif
		}
}		



void read_hubble()
{
	int i=0, n=0; 
	char dummyline[512];
	FILE *hub;

	fprintf(stderr, "\nread_hubble(). Using file: %s \n", Urls_internal.hubble_file);
	hub = fopen(Urls_internal.hubble_file, "r");
	n = get_lines(hub, Urls_internal.hubble_file);

	Cosmo.npts = n;
	Cosmo.a_hub = (double *) calloc(n, sizeof(double));
	Cosmo.z_hub = (double *) calloc(n, sizeof(double));
	Cosmo.Hubble = (double *) calloc(n, sizeof(double));

		for(i=0; i<n; i++)
		{
			fgets(dummyline, 512, hub);
			sscanf(dummyline, "%lf %lf", &Cosmo.a_hub[i], &Cosmo.Hubble[i]);
			Cosmo.z_hub[i] = 1./Cosmo.a_hub[i] - 1.;
		}

	fclose(hub);
}


		/* Initialize cosmology: read w(z) and H(z) from a file */
void read_cosmology(char *curl)
{
	int dim=0, k=0;
	char line[512];
	FILE *cfile = NULL;
	
	fprintf(stderr, "\nread_cosmology(). Initializing H(z) and w(z) from file:%s\n", curl);

	cfile=fopen(curl,"r");
	if(cfile==NULL) 
		fprintf(stderr, "Could not open file: %s.\n", curl);

		dim = get_lines(cfile, curl) - 1;

		Cosmo.npts = dim;
		Cosmo.z_hub = (double*) calloc(dim, sizeof(double));
		Cosmo.a_hub = (double*) calloc(dim, sizeof(double));
		Cosmo.Hubble = (double*) calloc(dim, sizeof(double));
		Cosmo.w = (double*) calloc(dim, sizeof(double));

			for(k=0; k<dim; k++)
			{
				fgets(line, 512, cfile);
				sscanf(line, "%lf  %lf  %lf", &Cosmo.z_hub[k], &Cosmo.Hubble[k], &Cosmo.w[k]);
				Cosmo.a_hub[k] = 1/(1+Cosmo.z_hub[k]);
			}

	fclose(cfile);
}

		/* A routine reading a mass function from a given file */
void read_numerical_mass_function(char *URL)
{
		int dim=0, k=0;
		char lines[50];
		FILE* nmf = fopen(URL, "r");

		fprintf(stderr, "\nread_numerical_mass_function() from file:%s\n", URL);

		dim = get_lines(nmf, URL);
		MF.bins  =  dim;
		MF.n = (double*) calloc(dim, sizeof(double));
		MF.num_masses  = (double*) calloc(dim, sizeof(double));

			for(k=0; k<dim; k++)
			{
				fgets(lines, 512, nmf);
				sscanf(lines, "%lf %lf", &MF.num_masses[k], &MF.n[k]);
			}
	fclose(nmf);
}
