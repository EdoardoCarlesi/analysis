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
#include "../libcosmo/cosmo.h"
#include "../libhalo/halo.h"
#include "../general_def.h"

#include "halo.h"

#ifdef WITH_MPI
#include "../libparallel/general.h"
#endif


/*
 *  Function declaration
 */
void fit_polytropic_T(struct halo *HALO);
void sort_f_gas_profile(struct halo *HALO);
void fit_I_X(struct halo *HALO);


/*
 *  Function implementations
 */
void fit_and_store_gas_parameters(void)
{	
	int k=0; 
	int nHaloes;

	struct halo * HALO;

#ifdef WITH_MPI
	HALO = pHaloes[ThisTask];
	nHaloes = pSettings[ThisTask].n_threshold;
	TASK_INFO_MSG(ThisTask, "Fitting gas profiles");
#else 
	HALO = Haloes;
	nHaloes = Settings.n_threshold;
	INFO_MSG("Fitting gas profiles");
#endif

		for(k=0; k<nHaloes; k++)
		{		
			//if(halo_condition(k) == 1)
			{	
				sort_f_gas_profile(&HALO[k]);
				fit_I_X(&HALO[k]);
			}
		}
}



void sort_f_gas_profile(struct halo *HALO)
{
	double h, R, rMax, rMin; 
	double *x, *y, *x_bin, *y_bin, *e_bin; 
	double chi, per, gof;
	int N, bins, skip, j=0;
	struct general_settings *SETTINGS;

#ifdef WITH_MPI
	SETTINGS = &pSettings[ThisTask];
#else
	SETTINGS = &Settings;
#endif

		bins = HALO->n_bins;
		skip = HALO->neg_r_bins;
		N = bins - skip;

		x = (double*) calloc(N, sizeof(double));
		y = (double*) calloc(N, sizeof(double));

		for(j=0; j<N; j++)
		{
			R = HALO->radius[j+skip];
			x[j] = R/HALO->Rvir;
			y[j] = HALO->gas_only.m[j+skip]/HALO->mass_r[j+skip];
		}
	
		x_bin = (double*) calloc(BIN_PROFILE+1, sizeof(double));
		y_bin = (double*) calloc(BIN_PROFILE, sizeof(double));
		e_bin = (double*) calloc(BIN_PROFILE, sizeof(double));

			rMin = Rvir_frac_min;
			rMax = F_MAX * R/HALO->Rvir;
			x_bin = log_stepper(rMin, rMax, BIN_PROFILE+1);

		average_bin(x, y, x_bin, y_bin, e_bin, BIN_PROFILE+1, N);

	for(j=0; j<BIN_PROFILE; j++)
	{
		HALO->f_gas.x[j] = 0.5 * (x_bin[j] + x_bin[j+1]);
		HALO->f_gas.y[j] = y_bin[j];
		//fprintf(stderr, "%d  %f  %f\n", j, HALO->f_gas.x[j], HALO->f_gas.y[j]);
	}
}



void average_gas_profiles(void)
{
	int k=0, i=0, m=0; 
	int nHaloes, nHaloesCut;
	int f_bin, rho_bin, ix_bin;
	double f, f_tot;
	double rho, rho_tot;
	double ix, ix_tot;

	struct halo_properties *HALOPROPERTIES;

	INFO_MSG("Sorting halo radial velocities and concentrations");

	nHaloesCut=n_haloes_per_criterion();
	
	if(Settings.use_sub == 1)
	{
		HALOPROPERTIES = SubHaloProperties;	
	} 
		else 
	{
		HALOPROPERTIES = HaloProperties;
	}

#ifdef WITH_MPI
		nHaloes=Settings.n_haloes; 
#else
		nHaloes=Settings.n_threshold; 
#endif

		for(i=0; i<BIN_PROFILE; i++)
		{
			f_tot = 0;
			f_bin = 0;
	
			rho_tot = 0;
			rho_bin = 0;
	
			ix_tot = 0;
			ix_bin = 0;

			for(k=0; k<nHaloes; k++)
			{
				if(halo_condition(k) == 1)
				{
					f = Haloes[k].f_gas.y[i];
					rho = Haloes[k].rho_gas.y[i];
					ix = Haloes[k].i_x.y[i];

						if(isnan(f) == 0 && f>0)
						{
							f_tot += f;
							f_bin++;
						}

						if(isnan(rho) == 0 && rho>0 && rho < 1./0.)
						{
							rho_tot += rho;
							rho_bin++;
						}
	
						if(isnan(ix) == 0 && ix>0 && ix < 1./0.)
						{
							ix_tot += ix;
							ix_bin++;
						}
				}
		
			}

			HALOPROPERTIES[HALO_INDEX].f_gas.x[i] = Haloes[0].f_gas.x[i];
			HALOPROPERTIES[HALO_INDEX].f_gas.y[i] = f_tot/f_bin;
			HALOPROPERTIES[HALO_INDEX].f_gas.n[i] = f_bin;
			HALOPROPERTIES[HALO_INDEX].rho_gas.x[i] = Haloes[0].rho_gas.x[i];
			HALOPROPERTIES[HALO_INDEX].rho_gas.y[i] = rho_tot/rho_bin;
			HALOPROPERTIES[HALO_INDEX].rho_gas.n[i] = rho_bin;
			HALOPROPERTIES[HALO_INDEX].i_x.x[i] = Haloes[0].i_x.x[i];
			HALOPROPERTIES[HALO_INDEX].i_x.y[i] = ix_tot/ix_bin;
			HALOPROPERTIES[HALO_INDEX].i_x.n[i] = ix_bin;

			if(Haloes[0].rho_gas.x[i]==0)
			{
				HALOPROPERTIES[HALO_INDEX].rho_gas.y[i] = 0; 
				HALOPROPERTIES[HALO_INDEX].i_x.y[i] = 0; 
			}
		}
}



double rhoBeta(double r, double rc, double beta, double rho0)
{
	return rho0 * pow(1 + pow2(r/rc), -3. * beta);
}



double I_X(double r, double rc, double beta, double rho0)
{
		// This is I_x / I_x0
	return (1./rho0) * rhoBeta(r, rc, beta, rho0) * pow(1 + pow2(r/rc), 0.5);
}


	// Fits X-ray surface and beta gas density profile
void fit_I_X(struct halo *HALO)
{
	double rc, A=0, a=0, rho0, ix0; 
	double *r, *x, *y, *e, *y_th, *params; 
	double rMax, rMin, *x_bin, *y_bin, *e_bin;
	double beta, chi, per, gof, R, V;

	int bins, skip, j=0, N=0, M=0;
	struct general_settings *SETTINGS;

#ifdef WITH_MPI
	SETTINGS = &pSettings[ThisTask];
#else
	SETTINGS = &Settings;
#endif

		A = 1.;
		a = -1.;

		bins = HALO->n_bins;
		skip = HALO->neg_r_bins;
		N=1;
		R = HALO->radius[skip];
		rho0 = (HALO->gas_only.m[skip]) / (4./3. * M_PI * pow3(R));

		params = (double*) calloc(2, sizeof(double));
		x = (double*) calloc(1, sizeof(double));
		y = (double*) calloc(1, sizeof(double));
		e = (double*) calloc(1, sizeof(double));
		r = (double*) calloc(1, sizeof(double));
		rc = 0.5 * HALO->r2; // The core radius is half the scale radius

		for(j=0; j<bins; j++)
		{
			R = HALO->radius[j];
			V = 4./3. * M_PI * pow3(R);

			if(R > Rvir_frac_min)
			{
				N++;
				x = (double*) realloc(x, N * sizeof(double));
				y = (double*) realloc(y, N * sizeof(double));
				e = (double*) realloc(e, N * sizeof(double));
				r = (double*) realloc(r, N * sizeof(double));

				r[N-1] = R/HALO->Rvir;
				x[N-1] = 1 + R*R/(rc*rc);
				y[N-1] = (HALO->gas_only.m[j])/ (4./3. * M_PI * pow3(R)) / rho0;
				e[N-1] = y[N-1]/sqrt(HALO->npart[j]);
			}

		}
			
			params[0] = a;
			params[1] = A;
			
			params = best_fit_power_law(x, y, e, N, params);
	
			a = params[0];
			A = params[1];
			beta = -2./3. * a;	

			y_th = (double*) calloc(N, sizeof(double));

			for(j=0; j<N; j++)
			{
				y_th[j] = A * pow(x[j], a);
			}

			x_bin = (double*) calloc(BIN_PROFILE+1, sizeof(double));
			y_bin = (double*) calloc(BIN_PROFILE, sizeof(double));
			e_bin = (double*) calloc(BIN_PROFILE, sizeof(double));

			rMin = 2 * Rvir_frac_min; 
			rMax = F_MAX * R / HALO->Rvir;

			x_bin = log_stepper(rMin, rMax, BIN_PROFILE+1);

		average_bin(r, y, x_bin, y_bin, e_bin, BIN_PROFILE+1, N);

		// Normalize to 1 in the central region	
		ix0 = 1./sqrt(1 + pow2(rMin * HALO->Rvir / rc));

	for(j=0; j<BIN_PROFILE; j++)
	{
		HALO->rho_gas.x[j] = 0.5 * (x_bin[j] + x_bin[j+1]);
		HALO->rho_gas.y[j] = y_bin[j];
		HALO->i_x.x[j] = 0.5 * (x_bin[j] + x_bin[j+1]);
		HALO->i_x.y[j] = ix0 * y_bin[j] * sqrt(1 + pow2(x_bin[j] * HALO->Rvir / rc));
		//fprintf(stderr, "%d  %f  %f\n", j, x_bin[j+1], y_bin[j]);
	}

	// Various estimators for the goodness of fit
	chi = chi_square(y_th, y, e, N);
	gof = goodness_of_fit(y_th, y, N);
	per = percentage_error(y_th, y, N);
	
	chi /= (double) N;

	HALO->fit_beta.rho0 = rho0;
	HALO->fit_beta.beta = beta;
	HALO->fit_IX.ix0 = ix0;
	HALO->fit_IX.beta = beta;
	HALO->fit_IX.chi = chi;
	HALO->fit_IX.gof = gof;
	HALO->fit_IX.per = per;
/*
*/
}



double polytropic_T(double T, double a, double A)
{	
	// T is in T0 units, returns rho in rho0 units
	// a = 1 / (gamma - 1)
	return A * pow(T, a);
}



void fit_polytropic_T(struct halo *HALO)
{
	double T=0, M=0, Mtot, MT, A=0, a=0, rho0, R, V, T0, V0, R0; 
	double *x, *y, *e, *y_th, *params; 
	double chi, per, gof;
	int bins, skip, j=0, N=0;
	struct general_settings *SETTINGS;
		// FIXME
/*
#ifdef WITH_MPI
	SETTINGS = &pSettings[ThisTask];
#else
	SETTINGS = &Settings;
#endif

		A = 1.;
		a = 5./3.;

		bins = HALO->n_bins;
		skip = HALO->neg_r_bins;
		N = bins - skip;
		R = HALO->radius[skip];
		V0 = 4./3. * M_PI * pow3(R);
		R0 = HALO->gas_only.m[skip] / V0;
		rho0 = R0; ///SETTINGS->rho_b;
		T0 = convert_u_to_T(HALO->gas_only.u[skip])/HALO->npart[skip]; // / V0;
		M = 0;
		T = 0;
		MT = 0;
		Mtot = 0;

		params = (double*) calloc(2, sizeof(double));
		x = (double*) calloc(N, sizeof(double));
		y = (double*) calloc(N, sizeof(double));
		e = (double*) calloc(N, sizeof(double));

		for(j=0; j<N; j++)
		{
			R = HALO->radius[j+skip];
			V = 4./3. * M_PI * pow3(R);

			if(j>0)
				M = HALO->gas_only.m[j+skip-1];

			T = convert_u_to_T(HALO->gas_only.u[j+skip]);
			MT += (HALO->gas_only.m[j+skip]-M)*T;

			x[j] = M/V/rho0;
			y[j] = MT/HALO->gas_only.m[j+skip]/T0;

	//		y[j] = HALO->gas_only.u[j+skip]; ///T0;
			e[j] = y[j]/sqrt(HALO->npart[j+skip]);

			fprintf(stderr, "%e %e\n", x[j], y[j]);
		}

			
			params[0] = a;
			params[1] = A;
			
			params = best_fit_power_law(x, y, e, N, params);
	
			a = params[0];
			A = params[1];
	
			y_th = (double*) calloc(N, sizeof(double));

			double gamma = 1 + 1./a;
			F_PRINT("gamma=", gamma);

		for(j=0; j<N; j++)
		{
			y_th[j] = polytropic_T(y[j], a, A);
		}

	// Various estimators for the goodness of fit
	chi = chi_square(y_th, y, e, N);
	gof = goodness_of_fit(y_th, y, N);
	per = percentage_error(y_th, y, N);
	
	chi /= (double) N;

	HALO->fit_poly.chi = chi;
	HALO->fit_poly.gof = gof;
	HALO->fit_poly.per = per;
*/
}
#endif
