#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../libio/io.h"
#include "../libmath/math.h"
#include "../general_def.h"

#include "io.h"


/*
 * Define functions
 */
void fill_pk_file(int);


/*
 * Init functions
 */
void fill_pk_file(int num)
{
	FILE *f = fopen(Urls.urls_pks[num], "r");
	double kk=0, ppkk=0;
	char dummyline[200];
	int i=0, skip=Settings.pk_skip;

		for(i=0; i< Pks[num].npts + skip; i++) 
			{
				fgets(dummyline, 200, f);
			if(i>=skip)
			{
				sscanf(dummyline, "%lf  %lf", &kk, &ppkk);
#ifdef USE_UNIT_MPC
				kk *= 1.e3;
				ppkk *= 1.e-9;
#endif
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
	fprintf(stdout, "\nget_pk_snapshot()\n");
#endif

		if(Settings.use_one_pk == 1)
		{
			fprintf(stdout, "\nReading one P(k) from: %s\n", Urls.pk_file);
			Urls.urls_pks[0] = (char *) calloc(strlen(Urls.pk_file), sizeof(char));
			strcpy(Urls.urls_pks[0],Urls.pk_file);
			
		} else {

			fprintf(stdout, "\nReading pk file list from: %s\n", Urls.pk_list);
				pk = fopen(Urls.pk_list,"r");
		if(pk==NULL) 
				ERROR("File not found", Urls.pk_list);

		// Check that redshifts and snapshot number match
	nPkSnaps = get_lines(pk, Urls.pk_list);
	rewind(pk);

	if(nPkSnaps != Urls.nPkFiles) 
	{
		fprintf(stdout, "\nNumber of redshifts (%d) and number of snapshots (%d) does not match!\n", 
			Urls.nPkFiles, nPkSnaps);
	}

			for(n=0; n<nPkSnaps; n++)	
			{
				fgets(dummyline,120,pk);
				Urls.urls_pks[n] = (char *) calloc(strlen(dummyline)+1, sizeof(char));	
				strcpy(Urls.urls_pks[n],dummyline);
				Urls.urls_pks[n][strlen(dummyline)-1]='\0';
#ifdef PRINT_INFO
	fprintf(stdout, "%s\n", dummyline);
#endif
			}		
		} // if there is more than one P(k)
}



void init_pks()
{
	int m=0, numPkFiles=0, dimPkFile=0;  
	FILE *f=NULL;  
 
		if(Settings.use_one_pk == 1) 
		{
			numPkFiles = 1;	
		} else {

			numPkFiles = Urls.nPkFiles;
		}

		Pks = (struct power_spectrum *) calloc(numPkFiles, sizeof(struct power_spectrum));
		Urls.urls_pks = (char **) calloc(numPkFiles, sizeof(char *));
			
		read_pk_snapshots();

		for(m=0; m<numPkFiles; m++)
		{
			f = fopen(Urls.urls_pks[m], "r");
			dimPkFile = get_lines(f, Urls.urls_pks[m]) - Settings.pk_skip;
			Pks[m].npts = dimPkFile;

#ifdef PRINT_INFO
		fprintf(stdout, "Number of P(k) files:%d, dimension of file[%d]:%d \n", numPkFiles, m+1, dimPkFile);
#endif

				Pks[m].k    = (double *) calloc(dimPkFile, sizeof(double));
				Pks[m].pk   = (double *) calloc(dimPkFile, sizeof(double));

				if(Settings.use_one_pk > 1)
				{
					Pks[m].z = GrowthFac.z[m];
					Pks[m].a = GrowthFac.a[m];
				}

		fill_pk_file(m);
	}
}
