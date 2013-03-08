#ifdef GAS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>

#include "../libmath/math.h"
#include "../general_def.h"

#include "halo.h"

#ifdef WITH_MPI
#include "../libparallel/general.h"
#endif


/*
 *  Function declaration
 */





/*
 *  Function implementations
 */
double polytropic_T(double T, double a, double A)
{	
	// T is in T0 units, returns rho in rho0 units
	// a = 1 / (gamma - 1)
	return A * pow(T, a);
}


void fit_and_store_polytropic_T_parameters(void)
{
	double A=0, a=0, rho0, T0, V0, gof; 
	double *x, *y, *e, *y_th, *params; 
	int bins, skip, j=0;
	struct general_settings *SETTINGS;

#ifdef WITH_MPI
	SETTINGS = &pSettings[ThisTask];
#else
	SETTINGS = &Settings;
#endif

		A = 1.;
		a = 5./3.;

		bins = HALO->n_bins;
		skip = HALO->neg_r_bins;
		V0 = 4./3. * M_PI * pow3(HALO->radius[skip]);
		rho0 = HALO->rho[skip] / (V0 * SETTINGS->rho_b);
		T0 = convert_u_to_T(HALO->gas_only.u[skip]);

		params = (double*) calloc(2, sizeof(double));
		x = (double*) calloc(N, sizeof(double));
		y = (double*) calloc(N, sizeof(double));
		e = (double*) calloc(N, sizeof(double));
			
		params[0] = a;
		params[1] = A;

		for(j=0; j<N; j++)
		{
			x[j] = HALO->radius[j+skip];
			y[j] = HALO->rho[j+skip];
			e[j] = HALO->err[j+skip];
		}

			best_fit_power_law(x, y, e, N);

			HALO->fit.rho0_nfw = rho0;
			HALO->fit.rs_nfw = rs;
			HALO->fit.c_nfw = HALO->Rvir/rs;

			y_th = (double*) calloc(bins-skip,sizeof(double));

		for(j=skip; j<HALO->n_bins; j++)
		{
			y_th[j-skip] = nfw(HALO->radius[j], HALO->fit.rs_nfw, HALO->fit.rho0_nfw);
			//fprintf(stderr, "%d) %f  %f  %f\n", j, HALO->radius[j], HALO->rho[j], y_th[j-skip]);
		}

	// Various estimators for the goodness of fit
	chi = chi_square(y_th,HALO->rho,HALO->err,bins,skip);
	gof = goodness_of_fit(y_th,HALO->rho,bins,skip);
	per = percentage_error(y_th,HALO->rho,bins,skip);
	
	chi /= (double) (bins-skip);

	HALO->fit.chi_nfw = chi;
	HALO->fit.gof_nfw = gof;
	HALO->fit.per_nfw = per;

}
#endif
