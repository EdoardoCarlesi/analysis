#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "io.h"
#include "../general_def.h"
#include "../libcosmo/cosmo.h"



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
