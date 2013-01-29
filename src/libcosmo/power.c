#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>

#include "power.h"
#include "mass_function.h"
#include "../general_variables.h"
#include "../libio/power_io.h"
#include "../libmath/mathtools.h"


	/* Interpolated P(k) */
double power_k(double k, int index)
{
	return get_interpolated_value(Pks[index].k, Pks[index].pk, 
		Pks[index].npts, k);
}



void normalize_all_power_spectra_to_sigma8()
{
	int i=0, j=0, Npk=0, Nentr=0;
	Npk = Urls.nPkFiles;
	
		for(i=0; i<Npk; i++)
		{	
			Nentr = Pks[i].npts;

			for(j=0; j<Nentr; j++)
				Pks[i].pk[j] *= Cosmo.norm_sigma8;
		}

}



void compute_growth_factor()
{
	fprintf(stdout,"compute_growth_factor().\n");

	int i=0, j=0, pk0=0, pkNorm=0, dim_eff=0; 
	double norm=1., pk_norm=1., delta=1.;
	double a=0, z=0, pk=0, K=0, kk=0; 

		kk = Pks[0].k[0];
		K = GrowthFac.scale_k, 

		dim_eff = Urls.nPkFiles;
		pkNorm  = 0; // This way we normalize GrowthFactor(z=0) = 1

			if(K < kk) 
			{
				fprintf(stdout, 
					"\nChosen k (%lf) is smaller than the numerical range available.\n\
						Setting k=%lf.\n", K, kk);
					K = kk; 
			}
				// pk0 is the highest redshift power spectrum
				// pkNorm is the z at which we want the GrowthFac to be exactly one
				pk_norm = power_k(K, pk0);
				norm = power_k(K, pkNorm);

		fprintf(stdout, "Normalizing to Pk[%d]=%lf at scale k:%lf\n", pk0, pk_norm, K);

		for (j=0; j<dim_eff; j++)
		{
			i = dim_eff-j-1;
			pk = power_k(K, i);
			delta = sqrt(pk/pk_norm);
			delta *= sqrt(pk_norm/norm);
			a = GrowthFac.a[i];
			z = 1./a - 1;

#ifdef PRINT_INFO
		fprintf(stdout, "GrowthFac.z[%d]=%lf, delta:%lf, delta_a:%lf\n", i, z, delta, delta/a);
#endif

			GrowthFac.gf[i] = delta;
		}
}



void fit_correlation_function()
{
//TODO
}



void calculate_power_k_from_xi()
{
// TODO
}



void compute_correlation_function(int index)
{
	int i=0, N_r=0;
	double r_min=0.1, r_max=100;
	double r=0, vol=0, fit=0, r_0=13.; 
	double *R=NULL;

	fprintf(stdout, "Computing Xi(r) between r_min:%lf and r_max:%lf\n", r_min, r_max);

	Xi.npts = Settings.n_bins;
	N_r = Xi.npts;

	Xi.r = (double *) calloc(N_r, sizeof(double));
	Xi.xi_r = (double *) calloc(N_r, sizeof(double));
	Xi.xi_fit = (double *) calloc(N_r, sizeof(double));

		R = log_stepper(r_min, r_max, N_r);

		for(i=0; i<N_r; i++)
		{
			r = R[i];
			vol = 4*PI*r*r*r*0.3333;
			vol = 1.0;
			fit = pow(r/r_0, -1.8);
			Xi.r[i] = r;
			Xi.xi_r[i] = correlation_r(r,index)/vol;
			Xi.xi_fit[i] = fit;
#ifdef PRINT_INFO
	fprintf(stdout, "At r:%lf, Xi(r):%lf; fitted value:%lf \n", Xi.r[i], Xi.xi_r[i], fit);
#endif
		}
}



double correlation_r(double r, int index)
{
	int WORKSP = 5000, nPkFiles=0; 
	double result=0., error=0., k_min=0., k_max=0.; 
	double par[2]; 
	
	nPkFiles = Urls.nPkFiles;

	par[0] = r; 
	par[1] = (double) index;

	k_min = Pks[index].k[0];
	k_max = Pks[index].k[nPkFiles-1];

#ifdef PRINT_INFO
	fprintf(stdout, "correlation_r(), r:%lf\n", r);
#endif

		gsl_integration_workspace *w = gsl_integration_workspace_alloc(WORKSP);
		gsl_function F;
		F.function = &correlation_integral;
		F.params = &par;
		gsl_integration_qags(&F, k_min, k_max, 0, 1e-6, WORKSP, w, &result, &error);

		result*=(1./(2*PI*PI));

		gsl_integration_workspace_free (w);

	return result;
}



double correlation_integral(double k, void *p)
{
	int i=0;
	double r, ind, *par;

		par = (double*) p;
		r   = par[0];
		ind = par[1];
		i = (int) ind;

	return power_k(k, i) * (sin(k*r)/k*r);
}
