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



/*
 *  Function implementations
 */
void fit_and_store_polytropic_T_parameters(void)
{	
	int k=0; 
	int nHaloes;

	struct halo * HALO;

#ifdef WITH_MPI
	HALO = pHaloes[ThisTask];
	nHaloes = pSettings[ThisTask].n_threshold;
	TASK_INFO_MSG(ThisTask, "fitting halo T profiles");
#else 
	HALO = Haloes;
	nHaloes = Settings.n_threshold;
	INFO_MSG("Fitting halo T profiles");
#endif

		for(k=0; k<nHaloes; k++)
		{		
			//fit_polytropic_T(&HALO[k]);			
			sort_f_gas_profile(&HALO[k]);
		}
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
#ifdef USE_UNIT_MPC
		R *= 1.e-3;	
#endif
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
#ifdef USE_UNIT_MPC
			R *= 1.e-3;	
#endif
			V = 4./3. * M_PI * pow3(R);

			if(j>0)
				M = HALO->gas_only.m[j+skip-1];

			T = convert_u_to_T(HALO->gas_only.u[j+skip]);
			MT += (HALO->gas_only.m[j+skip-1]-M)*T;

			x[j] = M/V/rho0;
			y[j] = MT/M/T0;

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
#ifdef USE_UNIT_MPC
			R *= 1.e-3;	
#endif
			x[j] = R/HALO->Rvir;
			y[j] = HALO->gas_only.m[j+skip]/HALO->mass_r[j+skip];
		}
	
		x_bin = (double*) calloc(BIN_PROFILE+1, sizeof(double));
		y_bin = (double*) calloc(BIN_PROFILE, sizeof(double));
		e_bin = (double*) calloc(BIN_PROFILE, sizeof(double));

			rMin = 0.1;
			rMax = 1.;
			x_bin = log_stepper(rMin, rMax, BIN_PROFILE+1);

		average_bin(x, y, x_bin, y_bin, e_bin, BIN_PROFILE+1, N);

		h = 0.5 * (x_bin[1] - x_bin[0]); 

	for(j=0; j<BIN_PROFILE; j++)
	{
		HALO->f_gas.x[j] = x_bin[j+1];// - h;
		HALO->f_gas.y[j] = y_bin[j];
		//fprintf(stderr, "%d  %f  %f\n", j, HALO->f_gas.x[j], HALO->f_gas.y[j]);
	}
}



void average_gas_fraction_profile(void)
{
	int k=0, i=0, m=0, n_bin=0; 
	int nHaloes, nHaloesCut;
	double f, f_tot, *m_bin;
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
			n_bin = 0;
			m = 0;

			for(k=0; k<nHaloes; k++)
			{
				if(halo_condition(k) == 1)
				{
					f = Haloes[m].f_gas.y[i];
					//F_PRINT("y", f);
						if(isnan(f) == 0 && f>0)
						{
							f_tot += f;
							n_bin++;
						}
					m++;
				}
		
			}
			
			HALOPROPERTIES[HALO_INDEX].f_gas.x[i] = Haloes[0].f_gas.x[i];
			HALOPROPERTIES[HALO_INDEX].f_gas.y[i] = f_tot/n_bin;
		}
}
#endif
