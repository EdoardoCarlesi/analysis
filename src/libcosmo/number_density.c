#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tinker.h"
#include "mass_function.h"
#include "number_density.h"
#include "cosmological_relations.h"

#include "../libio/power_io.h"
#include "../libio/halo_io.h"
#include "../libmath/mathtools.h"
#include "../general_functions.h"
#include "../general_variables.h"

void get_n_M_z(double M)
{
	int k=0, nTot=0, nPks=0;
	
	NumDen.npts = snaps_within_z_range();
	nTot = NumDen.npts;
	nPks = Urls.nPkFiles;	

	NumDen.z = (double*) calloc(nTot, sizeof(double));
	NumDen.n_th = (double*) calloc(nTot, sizeof(double));
	NumDen.n_num = (double*) calloc(nTot, sizeof(double));

	for(k=0; k<NumDen.npts; k++)
	{
		NumDen.z[k] = Pks[nPks-k-1].z;

		NumDen.n_th[k] = 
			get_interpolated_value(MassFuncZ[k].mass, 
				MassFuncZ[k].n, MassFuncZ[k].bins, M);

		NumDen.n_num[k] = 
			get_interpolated_value(MassFuncZ[k].mass, 
				MassFuncZ[k].n, MassFuncZ[k].bins, M);
	}
}



int snaps_within_z_range()
{
	int k=0;

	do 
	{ 
		k++;

	} while(NumDen.z[k]<NumDen.zMax);

	return k;
}



double integrand_tin_n_M_z(double z, void *p)
{
	return comoving_vol(z,p)*get_interpolated_value(NumDen.z, 
		NumDen.n_th, NumDen.npts, z);;
}



double integrand_num_n_M_z(double z, void *p)
{
	return comoving_vol(z,p)*get_interpolated_value(NumDen.z, 
		NumDen.n_num, NumDen.npts, z);;
}



double* integrate_number_density(double z0, double z1)
{
	double result0=0, error0, result1=0, error1, *result; 

	result = (double*) calloc(2,sizeof(double));

		gsl_function H, I;

		gsl_integration_workspace *w0 = gsl_integration_workspace_alloc(1000);
		gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(1000);

			H.function=&integrand_tin_n_M_z;
			I.function=&integrand_num_n_M_z;
			H.params=&z0;
			I.params=&z0;

		gsl_integration_qags(&H, z0, z1, 0, 1e-4, 1000, w0, &result0, &error0);
		gsl_integration_qags(&I, z0, z1, 0, 1e-4, 1000, w1, &result1, &error1);

		gsl_integration_workspace_free (w0);
		gsl_integration_workspace_free (w1);
			
		result[0] = result0;
		result[1] = result1;
	
	return result;
}



void compute_number_density()
{
	int index=0, m=0, nPk=0, nBin=0, nFc=0; 
	double zz=0;

	fprintf(stderr, "\ncompute_number_density().\n");

	nPk = Urls.nPkFiles;
	nFc = Urls.nCatalogueFiles;
	nBin = MassFunc.bins;
	
	init_pks();
	
	initialize_mass_function_datastore();

		do {
			index = nPk-m-1;
			zz=Pks[index].z;
			Settings.zStart = zz;

			fprintf(stderr,"\nstep %d, analyzing snapshot at redshift z: %lf\n", m, zz);
	
			Urls.halo_file = (char *) calloc(strlen(Urls.urls[nFc-m-1]), sizeof(char));
 			strcpy(Urls.halo_file, Urls.urls[nFc-m-1]);
			read_halo_file();

#ifndef TH_ONLY
			compute_numerical_mass_function();
#endif
			compute_theoretical_mass_function();

			ThMassFunc.Mmin=0.1*MassFunc.mass[0];
			ThMassFunc.Mmax=10*MassFunc.mass[nBin-1];
			store_mf(m);

			fprintf(stderr, "Done z: %lf. zMax: %lf \n", zz, NumDen.zMax);
			m++;

		} while(zz<NumDen.zMax);
}


