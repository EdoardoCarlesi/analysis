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
#include "../libio/io.h"
#include "../general_def.h"

#include "halo.h"

#ifdef WITH_MPI
#include "../libparallel/general.h"
#endif

#define	m_p (1.672e-24) // cgs
#define	mu (0.588) 
#define	k_b (1.3806e-16) // cgs
#define GAMMA (5./3.)
#define GN (6.672e-8) // cgs
#define SOLAR_MASS  1.989e33
#define CM_PER_MPC  3.085678e24
//#define Tfac 5.38
#define Tfac 11.76

/*
 *  Function declaration
 */
#ifndef NO_PROFILES
void fit_polytropic_T(struct halo *HALO);
void sort_f_gas_profile(struct halo *HALO);
void fit_I_X(struct halo *HALO);

double compute_hydrostatic_mass(struct halo *HALO, int);
#else
void average_gas_profiles()
{
}
#endif

/*
 *  Function implementations
 */
double u2TK(double u)
{
	return u*(GAMMA-1)*m_p*mu*(Tfac/k_b)*1.e+10; // Kelvin
}



double u2TkeV(double u)
{
		// Kelvin to eV factor = 8.6173e-8
	return u2TK(u)*(8.6173e-8); // KeV
}


#ifndef NO_PROFILES
void fit_and_store_gas_parameters(void)
{	
	int k=0; 
	int nHaloes;

	struct halo * HALO;

#ifdef WITH_MPI
	HALO = pHaloes[ThisTask];
	nHaloes = pSettings[ThisTask].n_threshold;
	if(ThisTask == 0)
	TASK_INFO_MSG(ThisTask, "Fitting gas profiles");
#else 
	HALO = Haloes;
	nHaloes = Settings.n_threshold;
	INFO_MSG("Fitting gas profiles");
#endif

		for(k=0; k<nHaloes; k++)
		{		
			{	
				sort_f_gas_profile(&HALO[k]);
				fit_I_X(&HALO[k]);
				fit_polytropic_T(&HALO[k]);
#if PRINT_HALO					
				if(HALO[k].Mvir > Settings.Mprint)
					print_halo_profile(k);
#endif

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
			HALO->gas_only.frac[j+skip] = y[j];
		}
	
			x_bin = (double*) calloc(BIN_PROFILE+1, sizeof(double));
			y_bin = (double*) calloc(BIN_PROFILE, sizeof(double));
			e_bin = (double*) calloc(BIN_PROFILE, sizeof(double));

			rMin = 2 * Rvir_frac_min;
			rMax = F_MAX * 1.01; //R/HALO->Rvir;
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
	int f_bin, rho_bin, ix_bin, temp_bin, pressure_bin, hydro_m_bin;
	double f, f_tot;
	double rho, rho_tot;
	double ix, ix_tot;
	double temp, temp_tot;
	double hydro_m, hydro_m_tot;
	double pressure, pressure_tot;

	INFO_MSG("Sorting halo gas profiles");

	nHaloesCut=n_haloes_per_criterion();

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

			temp_tot = 0;
			temp_bin = 0;

			for(k=0; k<nHaloes; k++)
			{
				if(halo_condition(k) == 1)
				{

					f = Haloes[k].f_gas.y[i];
					rho = Haloes[k].rho_gas.y[i];
					ix = Haloes[k].i_x.y[i];
					temp = Haloes[k].temp.y[i];
					pressure = Haloes[k].pressure.y[i];
					hydro_m = Haloes[k].hydro_m.y[i];

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

						if(isnan(temp) == 0 && temp>0 && temp < 1./0.)
						{
							temp_tot += temp;
							temp_bin++;
						}

						if(isnan(hydro_m) == 0 && hydro_m>0 && hydro_m < 1./0.)
						{
							hydro_m_tot += hydro_m;
							hydro_m_bin++;
						}
						
						if(isnan(pressure) == 0 && pressure>0 && pressure < 1./0.)
						{
							pressure_tot += pressure;
							pressure_bin++;
						}
				}
		
			}

			HaloProperties[HALO_INDEX].f_gas.x[i] = Haloes[0].f_gas.x[i];
			HaloProperties[HALO_INDEX].f_gas.y[i] = f_tot/f_bin;
			HaloProperties[HALO_INDEX].f_gas.n[i] = f_bin;
			HaloProperties[HALO_INDEX].rho_gas.x[i] = Haloes[0].rho_gas.x[i];
			HaloProperties[HALO_INDEX].rho_gas.y[i] = rho_tot/rho_bin;
			HaloProperties[HALO_INDEX].rho_gas.n[i] = rho_bin;
			HaloProperties[HALO_INDEX].i_x.x[i] = Haloes[0].i_x.x[i];
			HaloProperties[HALO_INDEX].i_x.y[i] = ix_tot/ix_bin;
			HaloProperties[HALO_INDEX].i_x.n[i] = ix_bin;
			HaloProperties[HALO_INDEX].temp.x[i] = Haloes[0].temp.x[i];
			HaloProperties[HALO_INDEX].temp.y[i] = temp_tot/temp_bin;
			HaloProperties[HALO_INDEX].temp.n[i] = temp_bin;
			HaloProperties[HALO_INDEX].pressure.x[i] = Haloes[0].pressure.x[i];
			HaloProperties[HALO_INDEX].pressure.y[i] = pressure_tot/pressure_bin;
			HaloProperties[HALO_INDEX].pressure.n[i] = pressure_bin;
			HaloProperties[HALO_INDEX].hydro_m.x[i] = Haloes[0].hydro_m.x[i];
			HaloProperties[HALO_INDEX].hydro_m.y[i] = hydro_m_tot/hydro_m_bin;
			HaloProperties[HALO_INDEX].hydro_m.n[i] = hydro_m_bin;

			if(Haloes[0].rho_gas.x[i]==0)
			{
				HaloProperties[HALO_INDEX].rho_gas.y[i] = 0; 
				HaloProperties[HALO_INDEX].i_x.y[i] = 0; 
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
	double rc, Np=0, u=0, T=0, A=0, a=0, rho, rho0, ix0, r500, *inv_r, *inv_rho, *r_hydro, *m_hydro;
	double *r, *x, *y, *m, *p, *e, *y_th, *params; 
	double rMax, rMin, *x_bin, *y_bin, *e_bin, *m_bin, *p_bin;
	double beta, chi, per, gof, R, V;
	double T_ew, T_ed, T_sl, T_sd, rho_p;

	int bins, skip, i=0, j=0, N=0, M=0;
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
		m = (double*) calloc(1, sizeof(double));
		p = (double*) calloc(1, sizeof(double));
		e = (double*) calloc(1, sizeof(double));
		r = (double*) calloc(1, sizeof(double));
		rc = 0.5 * HALO->r2; // The core radius is half the scale radius
		ix0 = 1./sqrt(1 + pow2(HALO->radius[skip] / rc))/rho0;

		for(j=0; j<bins; j++)
		{
			R = HALO->radius[j];
			V = 4./3. * M_PI * pow3(R);
			rho = (HALO->gas_only.m[j]) / V; //(4./3. * M_PI * pow3(R));

			if(j >= skip)
			{
				// Mass weighted T
				Np = HALO->gas_only.m[j]/SETTINGS->gasMass;
				rho_p = (float) rho / rho0;
				T += u2TK(HALO->gas_only.u[j]);
				HALO->gas_only.T[j] = T/Np; //u2T(HALO->gas_only.u[j]);
				HALO->gas_only.hydro_m[j] = compute_hydrostatic_mass(HALO,j);
				HALO->gas_only.pressure[j] = (float) pow(rho/rho0, 5./3.);
				//HALO->gas_only.rho[j] = (float) rho / rho0;
				HALO->gas_only.rho[j] = rho_p; 
				HALO->gas_only.i_x[j] = (float) ix0 * rho * sqrt(1 + pow2(R / rc));
				u += HALO->gas_only.u[j];			
				// Emission weighted T
				T_ew  += (T/Np) * sqrt(T/Np);// * rho;
				T_ed  += sqrt(T/Np);// * rho;
				// Spectroscopic like T
				T_sl  += pow(T/Np, 0.25);// * rho;
				T_sd  += pow(T/Np, -0.75);// * rho;
			}

			if(R > Rvir_frac_min)
			{
				N++;
				x = (double*) realloc(x, N * sizeof(double));
				y = (double*) realloc(y, N * sizeof(double));
				m = (double*) realloc(m, N * sizeof(double));
				p = (double*) realloc(p, N * sizeof(double));
				e = (double*) realloc(e, N * sizeof(double));
				r = (double*) realloc(r, N * sizeof(double));

				r[N-1] = R/HALO->Rvir;
				x[N-1] = 1 + (R*R/(rc*rc));
				y[N-1] = rho / rho0; //(HALO->gas_only.m[j])/ (4./3. * M_PI * pow3(R)) / rho0;
				m[N-1] = ((HALO->gas_only.hydro_m[j]-HALO->mass_r[j])/HALO->mass_r[j]); 
				p[N-1] = HALO->gas_only.pressure[j];
				e[N-1] = y[N-1]/sqrt(HALO->npart[j]);

		//		fprintf(stderr, "hydro=%e, mass=%e, m=%lf\n", HALO->gas_only.hydro_m[j], 
		//			HALO->mass_r[j], m[N-1]);
			}

		}
		
#ifdef HYDRO_500
			int SIZE = 6;
			inv_rho = calloc (SIZE, sizeof(double));
			inv_r   = calloc (SIZE, sizeof(double));

			for(i=0; i<SIZE; i++)
			{
				inv_rho[i] = HALO->rho[bins-i-1];
				//inv_r[i]   = HALO->radius[bins-i-1];
				inv_r[i]   = HALO->gas_only.hydro_m[bins-i-1];
				//fprintf(stderr, "%d) inv_rho=%f, r=%f\n", i, inv_rho[i], inv_r[i]);
			}

			//r500 = get_interpolated_value(inv_rho, inv_r, SIZE, 500);
			//HALO->gas_only.M_hydro = get_interpolated_value(x, HALO->gas_only.hydro_m, bins, r500);
			//HALO->gas_only.M_hydro = r500;
			HALO->gas_only.M_hydro = HALO->gas_only.hydro_m[bins-5];
#else
			HALO->gas_only.M_hydro = HALO->gas_only.hydro_m[bins-1];
#endif
	
			params[0] = a;
			params[1] = A;

			HALO->gas_only.T_mw = u2TkeV(u) / HALO->gas.N;	
			HALO->gas_only.T_ew = T_ew / T_ed / 1e7; 
			HALO->gas_only.T_sl = T_sl / T_sd / 1e7; 
		
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
			m_bin = (double*) calloc(BIN_PROFILE, sizeof(double));
			p_bin = (double*) calloc(BIN_PROFILE, sizeof(double));
			e_bin = (double*) calloc(BIN_PROFILE, sizeof(double));

			rMin = 2 * Rvir_frac_min; 
			rMax = F_MAX * R / HALO->Rvir;

			x_bin = log_stepper(rMin, rMax, BIN_PROFILE+1);

		average_bin(r, y, x_bin, y_bin, e_bin, BIN_PROFILE+1, N);
		average_bin(r, p, x_bin, p_bin, e_bin, BIN_PROFILE+1, N);
		average_bin(r, m, x_bin, m_bin, e_bin, BIN_PROFILE+1, N);

		// Normalize to 1 in the central region	
		R = 0.5*(x_bin[0]+x_bin[1]);
		ix0 = 1./sqrt(1 + pow2(R * HALO->Rvir / rc))/y_bin[0];

	for(j=0; j<BIN_PROFILE; j++)
	{
		R = 0.5 * (x_bin[j] + x_bin[j+1]);
		HALO->rho_gas.x[j] = R;
		HALO->rho_gas.y[j] = y_bin[j];
		HALO->pressure.x[j] = R;
		HALO->pressure.y[j] = pow(y_bin[j], 1.66666);
		HALO->hydro_m.x[j] = R;
		HALO->hydro_m.y[j] = m_bin[j];
		HALO->i_x.x[j] = R; 
		HALO->i_x.y[j] = ix0 * y_bin[j] * sqrt(1 + pow2(R * HALO->Rvir / rc));
	//	fprintf(stderr, "%d  %f  mbin = %f\n", j, x_bin[j+1], m_bin[j]);
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
}


double compute_hydrostatic_mass(struct halo *HALO, int index)
{
	struct general_settings *SETTINGS;
	long double r1, r2, rho1, rho2, T1, T2, V1, V2, N1, N2;
	double d_ln_r, d_ln_rho, d_ln_T;
	double Mr, r, T;

#ifdef WITH_MPI
	SETTINGS = &pSettings[ThisTask];
#else
	SETTINGS = &Settings;
#endif

		r1 = HALO->radius[index-1];
		r2 = HALO->radius[index];

		V1 = 4./3. * M_PI * pow3(r1);
		V2 = 4./3. * M_PI * pow3(r2);

		rho1 = HALO->gas_only.m[index-1]/V1;
		rho2 = HALO->gas_only.m[index]/V2;

		// Mass weighted temperature
		T1 = HALO->gas_only.T[index-1];
		T2 = HALO->gas_only.T[index];

		r = 0.5 * (r1+r2) * CM_PER_MPC;
		T = 0.5 * (T2+T1);

		// d ln(x) = delta x / x
		d_ln_rho = (rho2-rho1)/(rho2+rho1);
		d_ln_r = (r2-r1)/(r2+r1);
		d_ln_T = (T2-T1)/(T2+T1);

		Mr = - (k_b*T*r)/(mu*m_p*GN) * (d_ln_rho/d_ln_r + d_ln_T/d_ln_r) / SOLAR_MASS;
	
	return Mr;
}



double polytropic_T(double T, double a, double A)
{	
	// T must be T0 units, returns rho in rho0 units
	// a = 1 / (gamma - 1)
	return A * pow(T, a);
}



void fit_polytropic_T(struct halo *HALO)
{
	double T1=0, T2=0, T=0, M=0, Mtot, MT, A=0, a=0, rho0, R, V, T0, V0, R0; 
	double *r, *x, *y, *e, *y_th, *params, *t; 
	double *x_bin, *y_bin, *e_bin;
	double rMin, rMax, chi, per, gof, gamma_ply, u_cum, n_part;
	int bins, skip, n_eff, j=0, N=0;
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
		t = calloc(bins, sizeof(double));
		r = calloc(bins, sizeof(double));
		u_cum = 0;

		for (j=1; j<bins; j++)
		{
			u_cum = HALO->gas_only.u[j];
			n_part = (HALO->gas_only.m[j] - HALO->gas_only.m[j-1]) / SETTINGS->gasMass;
			HALO->gas_only.T[j] = u2TK(u_cum);// / n_part;
			t[j] = HALO->gas_only.T[j] / n_part;
			r[j] = HALO->radius[j] / HALO->Rvir;
		//	fprintf(stderr, "%d) npart=%f, r=%f, T=%e\n", n_part, r[j], t[j]);
		}

		for (j=1; j<bins; j++)
		{
			u_cum = HALO->gas_only.u[j];
			n_part = (HALO->gas_only.m[j] - HALO->gas_only.m[j-1]) / SETTINGS->gasMass;
			HALO->gas_only.T[j] = u2TK(u_cum);// / n_part;
			t[j] = HALO->gas_only.T[j] / n_part;
			r[j] = HALO->radius[j] / HALO->Rvir;
		//	fprintf(stderr, "%d) npart=%f, r=%f, T=%e\n", n_part, r[j], t[j]);
		}

/* //TODO fix this mess
		for (j=0; j<bins-1; j++)
		{
			T1 = u2TK(HALO->gas_only.u[j-1]);
			T = u2TK(HALO->gas_only.u[j]);
			T2 = u2TK(HALO->gas_only.u[j+1]);

			if(T > T1 && T > T2)
			{
				n_eff = j;
			}
		}
		
		//n_eff = skip;

		HALO->gas_only.T_0 = HALO->gas_only.T[n_eff]; 
		T0 = HALO->gas_only.T_0;
		R = HALO->radius[n_eff];
		M = HALO->gas_only.m[n_eff];
		V = 4./3. * M_PI * pow3(R);
		rho0 = M/V;

		N = bins - n_eff;

	//	fprintf(stderr, "%d, T1=%e T=%e T2=%e neff=%d, N=%d\n", j, T1, T, T2, n_eff, N);

	// Fit only with enough bins
	if(N>2)
	{
		params = (double*) calloc(2, sizeof(double));
		r = (double*) calloc(N, sizeof(double));
		x = (double*) calloc(N, sizeof(double));
		y = (double*) calloc(N, sizeof(double));
		e = (double*) calloc(N, sizeof(double));

		for(j=n_eff; j<bins; j++)
		{
			R = HALO->radius[j];
			T = HALO->gas_only.T[j];
			M = HALO->gas_only.m[j];
			V = 4./3. * M_PI * pow3(R);

			T = HALO->gas_only.T[j];

			x[j-n_eff] = T/T0;
			r[j-n_eff] = R/HALO->Rvir;
			y[j-n_eff] = (M/V)/rho0;
			e[j-n_eff] = y[j-n_eff]/sqrt((M-HALO->gas_only.m[j-1])/SETTINGS->gasMass);

			//fprintf(stderr, "T/T0=%e rho/rho0=%e e=%e\n", x[j-n_eff], y[j-n_eff], e[j-n_eff]);
			//fprintf(stderr, "rho0=%e M=%e V=%e\n", rho0, M, V);
		}
			
			params[0] = a;
			params[1] = A;
			
			params = best_fit_power_law(x, y, e, N, params);
	
			a = params[0];
			A = params[1];
	
			y_th = (double*) calloc(N, sizeof(double));

	//	double ga = 1 + 1./a;
	//	F_PRINT("********************gamma=", ga);
	//	F_PRINT("********************A=", A);

		for(j=0; j<N; j++)
		{
			y_th[j] = polytropic_T(y[j], a, A);
		}

		// Various estimators for the goodness of fit
		chi = chi_square(y_th, y, e, N);
		chi /= (double) N;
		//gof = goodness_of_fit(y_th, y, N); chi = gof;
		//per = percentage_error(y_th, y, N);
		gamma_ply = 1. + 1./a;

		// Sanity check - remove largely different profiles (inner parts might have been fitted instead of outer)
		if((gamma_ply > 2) || (gamma_ply < 0.5) || (abs(1-A) > 0.4) || (gamma_ply == 1.6) || (chi > 1e+2))
		{
			chi = -1;
			gamma_ply = -1;
		}
		else
		{
		//	fprintf(stderr, "%d) A=%f gamma=%f chi=%f\n", N, A, gamma_ply, chi);
		}

	} // If has less bins
	else
	{
		chi = -1;
		gamma_ply = -1;
	}

		HALO->fit_poly.chi = chi;
		HALO->fit_poly.gamma = gamma_ply;
		HALO->fit_poly.A = A;
*/
		x_bin = (double*) calloc(BIN_PROFILE+1, sizeof(double));
		y_bin = (double*) calloc(BIN_PROFILE, sizeof(double));
		e_bin = (double*) calloc(BIN_PROFILE, sizeof(double));
		rMin = 2 * Rvir_frac_min;
		rMax = F_MAX * 1.01; //R/HALO->Rvir;
		x_bin = log_stepper(rMin, rMax, BIN_PROFILE+1);
	
		// x is T/T0
		//average_bin(r, x, x_bin, y_bin, e_bin, BIN_PROFILE+1, N);
		average_bin(r, t, x_bin, y_bin, e_bin, BIN_PROFILE+1, bins);

		T0 = average(y_bin, BIN_PROFILE);

	for(j=0; j<BIN_PROFILE; j++)
	{
		//HALO->temp.x[j] = 0.5 * (x_bin[j] + x_bin[j+1]);
		HALO->temp.y[j] = y_bin[j]/T0;
	//	fprintf(stderr, "%d  %f\n", j, y_bin[j]/T0);
	}

		//HALO->fit_poly.gof = gof;
		//HALO->fit_poly.per = per;

}
#endif // NO_PROFILES
#endif
