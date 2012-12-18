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
		int k=0;
		ND.npts = snaps_within_z_range();
		ND.z = (double*) calloc(ND.npts, sizeof(double));
		ND.n_tin = (double*) calloc(ND.npts, sizeof(double));
		ND.n_num = (double*) calloc(ND.npts, sizeof(double));

		for(k=0; k<ND.npts; k++){
			ND.z[k] = Pks[Settings.n_pk_files-k-1].z;
			ND.n_tin[k] = get_interpolated_value(mf[k].th_masses, mf[k].tin, mf[k].bins, M);
			ND.n_num[k] = get_interpolated_value(mf[k].num_masses, mf[k].n, mf[k].bins, M);
		}
}



int snaps_within_z_range()
{
	int k=0;
		do { 
			k++;

			} while(ND.z[k]<ND.zMax);
	return k;
}



double integrand_tin_n_M_z(double z, void *p)
{
	return comoving_vol(z,p)*get_interpolated_value(ND.z, ND.n_tin, ND.npts, z);;
}



double integrand_num_n_M_z(double z, void *p)
{
	return comoving_vol(z,p)*get_interpolated_value(ND.z, ND.n_num, ND.npts, z);;
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
	int index=0, m=0, nPk=0, nBin=0, nFc=0; double zz;

	fprintf(stderr, "\ncompute_number_density().\n");

	nPk = Pks[0].n_pk_entries;
	nFc = FC.numFiles;
	nBin= MF.bins;
	
	init_pks();
	
	initialize_mass_function_datastore();

		do {
			index = nPk-m-1;
			zz=Pks[index].z;
			Settings.zStart = zz;

			fprintf(stderr,"\nstep %d, analyzing snapshot at redshift z: %lf\n", m, zz);

 			Urls_internal.halo_file = FC.urls[nFc-m-1];
			read_halo_file();

#ifndef TH_ONLY
			compute_numerical_mass_function();
#endif
			compute_theoretical_mass_function(index);

			AMF.Mmin=0.1*MF.num_masses[0];
			AMF.Mmax=10*MF.num_masses[nBin-1];
			store_mf(m);

			fprintf(stderr, "Done z: %lf. zMax: %lf \n", zz, ND.zMax);
			m++;

		} while(zz<ND.zMax);
}


