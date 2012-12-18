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
	return get_interpolated_value(Pks[index].k, Pks[index].pk, Pks[index].n_pk_entries, k);
}



void normalize_all_power_spectra_to_sigma8()
{
	int i=0, j=0, Npk=0, Nentr=0;
	Npk = Settings.n_pk_files;
	
		for(i=0; i<Npk; i++)
		{	
			Nentr = Pks[i].n_pk_entries;
			for(j=0; j<Nentr; j++)
				Pks[i].pk[j] *= Cosmo.n_s8;
		}

}



void compute_growth_factor()
{
		fprintf(stderr,"compute_growth_factor().\n");

		int i=0, vecOne=0, dim_eff=0; 
		double norm=1., delta=1., delta_0=1., K, delta_a, a, z, pk, pk_norm, kk; 

			kk = Pks[0].k[0];
			K = GF.scale_k, 
			dim_eff = Settings.n_pk_files;
			vecOne = 2;

			if(K < kk) 
				{
				fprintf(stderr, 
				"\n** WARNING! **\n Chosen k (%lf) is smaller than the numerical range. Setting k=%lf.\n", 
				K, kk);
				K = kk; 
				}

				pk_norm = power_k(K, vecOne);
				norm = power_k(K, dim_eff-1);

		fprintf(stderr, "vecOne:%d, pk:%lf, norm:%lf, dimEff:%d, K:%lf\n", vecOne, pk_norm, norm, dim_eff, K);

		for (i=0; i<dim_eff; i++)
		{
			pk = power_k(K, vecOne);
			delta = sqrt(pk/pk_norm);
			delta *= sqrt(pk_norm/norm);

				if(i==vecOne) delta_0 = GF.a[vecOne]/delta;
	
			delta *= delta_0;
			a=Pks[i].a;
			z = 1./a - 1;
			delta_a = delta/a;

#ifdef PRINT_INFO
		fprintf(stderr, "%d) z:%lf delta:%lf, delta_a:%lf\n", i, z, delta, delta_a);
#endif

			GF.gf_z[i] = delta;
			GF.gf_over_a_z[i] = delta_a;
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
	double r_min = 0.1, r_max = 100;

	fprintf(stderr, "Computing Xi(r) between r_min:%lf and r_max:%lf\n", r_min, r_max);

	double r=0, vol=0, fit=0, r_0=13., *R;
	int i, N_r;

	Xi.n_xi_entries = Settings.n_bins;
	N_r = Xi.n_xi_entries;
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
	fprintf(stderr, "At r:%lf, Xi(r):%lf; fitted value:%lf \n", Xi.r[i], Xi.xi_r[i], fit);
#endif
		}
}



double correlation_r(double r, int index)
{
	int WORKSP = 5000; 
	double result, error, k_min, k_max, par[2]; 

	par[0] = r; 
	par[1] = (double) index;
	k_min = Pks[index].k[0];
	k_max = Pks[index].k[Pks[index].n_pk_entries-1];

#ifdef PRINT_INFO
	fprintf(stderr, "correlation_r(), r:%lf\n", r);
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
