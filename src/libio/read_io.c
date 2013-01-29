#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "read_io.h"
#include "../general_variables.h"
#include "../general_functions.h"
#include "../libcosmo/cosmological_relations.h"


int get_lines(FILE *f, char *url)
{
	int n=-1;
	char line[2048];
	
		if(f==NULL) 
			ERROR("File not found when trying to get the number of lines", url);

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

		redshifts_file = fopen(Urls.a_outputs, "r");

		if(redshifts_file == NULL) 
			ERROR("File not found", Urls.a_outputs);
		else 
			fprintf(stdout, "Redshift file opened correctly:%s\n", Urls.a_outputs);

			size = get_lines(redshifts_file, Urls.a_outputs);

			Urls.nPkFiles = size;
			GrowthFac.npts = size;
			GrowthFac.z = (double *) calloc(size, sizeof(double));
			GrowthFac.a = (double *) calloc(size, sizeof(double));
			GrowthFac.gf = (double *) calloc(size, sizeof(double));

				for(i=0; i<size; i++) 
				{
					fgets(dummyline, 50, redshifts_file);
					sscanf(dummyline, "%lf", &aa);
					GrowthFac.a[size-i-1] = aa;
					GrowthFac.z[size-i-1] = 1./aa - 1.;
#ifdef PRINT_INFO
					fprintf(stdout, "%d) a:%lf z:%lf\n", size-i-1, aa, GrowthFac.z[size-i-1]);
#endif
				}
}		



void read_hubble()
{
	int i=0, n=0; 
	char dummyline[512];
	FILE *hub;

	fprintf(stdout, "\nread_hubble(). Using file: %s \n", Urls.hubble_file);
	hub = fopen(Urls.hubble_file, "r");
	n = get_lines(hub, Urls.hubble_file);

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
	
	fprintf(stdout, "\nread_cosmology(). Initializing H(z) and w(z) from file:%s\n", curl);

	cfile=fopen(curl,"r");

		if(cfile==NULL) 
			ERROR("File not found", curl);

		dim = get_lines(cfile, curl) - 1;

		Cosmo.npts = dim;
		Cosmo.z_hub = (double*) calloc(dim, sizeof(double));
		Cosmo.a_hub = (double*) calloc(dim, sizeof(double));
		Cosmo.Hubble = (double*) calloc(dim, sizeof(double));
		Cosmo.w = (double*) calloc(dim, sizeof(double));

			for(k=0; k<dim; k++)
			{
				fgets(line, 512, cfile);
				sscanf(line, "%lf  %lf  %lf", &Cosmo.z_hub[k], 
					&Cosmo.Hubble[k], &Cosmo.w[k]);
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

		fprintf(stdout, "\nread_numerical_mass_function() from file:%s\n", URL);

		dim = get_lines(nmf, URL);
		MassFunc.bins  =  dim;
		MassFunc.n = (double*) calloc(dim, sizeof(double));
		MassFunc.mass  = (double*) calloc(dim, sizeof(double));

			for(k=0; k<dim; k++)
			{
				fgets(lines, 512, nmf);
				sscanf(lines, "%lf %lf", &MassFunc.mass[k], 
					&MassFunc.n[k]);
			}
	fclose(nmf);
}
