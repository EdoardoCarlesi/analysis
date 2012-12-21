#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "power_io.h"
#include "../libio/read_io.h"
#include "../libmath/mathtools.h"
#include "../general_variables.h"
#include "../general_functions.h"


void fill_pk_file(int num)
{
	FILE *f = fopen(Pks[num].url, "r");
	double kk=0, ppkk=0;
	char dummyline[200];
	int i=0, skip=Settings.pk_skip;

		for(i=0; i< Pks[num].n_pk_entries + skip; i++) 
			{
				fgets(dummyline, 200, f);
			if(i>=skip)
			{
				sscanf(dummyline, "%lf  %lf", &kk, &ppkk);
				Pks[num].k[i-skip] = kk;
				Pks[num].pk[i-skip] = ppkk;
			}
		}
}



void read_pk_snapshots()
{
	char dummyline[120]; 
	int n=0, nPkSnaps=0;
	FILE *pk=NULL;

#ifdef PRINT_INFO
	fprintf(stderr, "\nget_pk_snapshot()\n");
#endif

		if(Settings.use_one_pk == 1)
		{
			fprintf(stderr, "\nReading one P(k) from: %s\n", Urls_internal.pk_file);
				Pks[0].url = Urls_internal.pk_file;
			
		} else {

			fprintf(stderr, "\nReading pk file list from: %s\n", Urls_internal.pk_list);
				pk = fopen(Urls_internal.pk_list,"r");
		if(pk==NULL) 
				fprintf(stderr, "Could not find pk.list file: %s\n", Urls_internal.pk_list);

		// Check that redshifts and snapshot number match
	nPkSnaps = get_lines(pk, Urls_internal.pk_list);
	rewind(pk);

	if(nPkSnaps != Settings.n_pk_files) 
	{
		fprintf(stderr, "\nNumber of redshifts (%d) and number of snapshots (%d) does not match!\n", 
		Settings.n_pk_files, nPkSnaps);
	}

			for(n=0; n<nPkSnaps; n++)	
			{
				fgets(dummyline,120,pk);
				dummyline[strlen(dummyline)-1]='\0';
				Pks[n].url = (char *) calloc(strlen(dummyline)+1, sizeof(char));	
				strcpy(Pks[n].url,dummyline);
#ifdef PRINT_INFO
	fprintf(stdout, "%s\n", dummyline);
#endif
			}		
		} // if there is more than one P(k)
}



void init_pks()
{
	int m=0, numPkFiles=0, dimPkFile=0;  
	FILE *f;  
 
		if(Settings.use_one_pk == 1) 
			{
				numPkFiles = 1;

				} else {
					numPkFiles = Settings.n_pk_files;
				}

			Pks = (struct power_spectrum *) calloc(numPkFiles, sizeof(struct power_spectrum));

			read_pk_snapshots();

			for(m=0; m<numPkFiles; m++)
			{
				f = fopen(Pks[m].url, "r");
				dimPkFile = get_lines(f, Pks[m].url) - Settings.pk_skip;
				Pks[m].n_pk_entries = dimPkFile;

#ifdef PRINT_INFO
		fprintf(stderr, "Number of P(k) files:%d, dimension of file[%d]:%d \n", numPkFiles, m+1, dimPkFile);
#endif

			Pks[m].k    = (double *) calloc(dimPkFile, sizeof(double));
			Pks[m].pk   = (double *) calloc(dimPkFile, sizeof(double));

				if(Settings.use_one_pk > 1)
				{
					Pks[m].z = GF.z[m];
					Pks[m].a = GF.a[m];
				}

		fill_pk_file(m);
	}
}
