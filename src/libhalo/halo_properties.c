#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../general_variables.h"
#include "../general_functions.h"
#include "../libmath/mathtools.h"
#include "../libmath/log_norm.h"
#include "../libmath/power_law.h"
#include "../libcosmo/mass_function.h"
#include "../libcosmo/cosmological_relations.h"

#include "halo_properties.h"
#include "subhalo_general.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef WITH_MPI
#include "../libparallel/general.h"
#endif


void initialize_halo_properties_structure()
{
	int rBins, nBins;

	rBins=Settings.r_bins; 
	nBins=Settings.n_bins;

			// Halo axis alignment
		HaloZ.r_bins=rBins-1; rBins--;
		HaloZ.R = (double*) calloc(rBins, sizeof(double));
		HaloZ.Th_p = (double*) calloc(rBins, sizeof(double));
		HaloZ.Th_c = (double*) calloc(rBins, sizeof(double));
		HaloZ.N_pairs = (int*) calloc(rBins, sizeof(double));
			// Other halo properties
		HaloZ.n_bins=nBins-1; nBins--;
		HaloZ.c = (double*) calloc(nBins, sizeof(double));
		HaloZ.c_avg = (double*) calloc(nBins, sizeof(double));
		HaloZ.p_c = (double*) calloc(nBins, sizeof(double));
		HaloZ.l = (double*) calloc(nBins, sizeof(double));
		HaloZ.p_l = (double*) calloc(nBins, sizeof(double));
		HaloZ.err_p_l = (double*) calloc(nBins, sizeof(double));
		HaloZ.shape = (double*) calloc(nBins, sizeof(double));
		HaloZ.p_shape = (double*) calloc(nBins, sizeof(double));
		HaloZ.n_shape = (int*) calloc(nBins, sizeof(double));
		HaloZ.triax = (double*) calloc(nBins, sizeof(double));
		HaloZ.p_triax = (double*) calloc(nBins, sizeof(double));
		HaloZ.n_triax = (int*) calloc(nBins, sizeof(double));
		HaloZ.mass = (double*) calloc(nBins, sizeof(double));
		HaloZ.radVel = (double*) calloc(nBins, sizeof(double));
		HaloZ.err_radVel = (double*) calloc(nBins, sizeof(double));

#ifdef GAS
	HaloZ.gas_T = (double*) calloc(nBins, sizeof(double));
	HaloZ.gas_u = (double*) calloc(nBins, sizeof(double));
	HaloZ.gas_fraction = (double*) calloc(nBins, sizeof(double));
#endif
}



void free_halo_properties()
{
		free(HaloZ.c);
		free(HaloZ.p_c);
		free(HaloZ.l);
		free(HaloZ.p_l);
		free(HaloZ.err_p_l);
		free(HaloZ.p_shape);
		free(HaloZ.n_shape);
		free(HaloZ.shape);
		free(HaloZ.triax);
		free(HaloZ.p_triax);
		free(HaloZ.n_triax);
		free(HaloZ.mass);
		free(HaloZ.radVel);
		free(HaloZ.err_radVel);
		free(Haloes);
#ifdef GAS
		free(HaloZ.gas_T);
		free(HaloZ.gas_u);
		free(HaloZ.gas_fraction);
#endif 
}



void free_halo_profiles(int i)
{
	struct halo *HALO;

#ifdef WITH_MPI
		HALO = pHaloes[ThisTask];
#else
		HALO = Haloes;
#endif
		free(HALO[i].radius);
		free(HALO[i].rho);
		free(HALO[i].over_rho);
		free(HALO[i].err);
		free(HALO[i].over_err);
		free(HALO[i].bin);
		free(HALO[i].err_dn);

#ifdef WITH_GAS
	free(HALO[i].m_gas);
	free(HALO[i].u_gas);
#endif
}


	// TODO: add mass cut
void sort_axis_alignement()
{
	int *Nbins, j=0, k=0, i=0, max_haloes, nBins, skip;
	double *radius, *Rbins, *Abins, *Bbins; 
	double Rmin, Rmax, R, A, B;

	fprintf(stdout,"Computing halo major axis alignement angles for %d bins.\n", HaloZ.r_bins);
	
		max_haloes = Settings.n_haloes;
		skip = Settings.halo_skip; 
		nBins = Settings.r_bins; 

		radius = (double *) calloc(nBins, sizeof(double));
		Rbins = (double *) calloc(nBins-1, sizeof(double));
		Abins = (double *) calloc(nBins-1, sizeof(double));
		Bbins = (double *) calloc(nBins-1, sizeof(double));
		Nbins = (int *) calloc(nBins-1, sizeof(int));
	
		Rmin = Settings.Rmin; Rmax = Settings.Rmax;
		radius = log_stepper(Rmin,Rmax,nBins);

			for(i=0; i<nBins-1; i++)
				Rbins[i] = 0.5*(radius[i+1]+radius[i]);

#ifdef _OPENMP
		omp_set_num_threads(OMP_THREADS);
#endif
		
//#		pragma omp parallel for			\
		private(j, k, A, B, R)

			for(j=0; j<max_haloes; j++) 
			{
				for(k=j; k<max_haloes; k++)
				{
					A = 0; B = 0; R = 0;
					R = sqrt(
						pow(Haloes[j].Xc - Haloes[k].Xc,2) +
						pow(Haloes[j].Yc - Haloes[k].Yc,2) +
						pow(Haloes[j].Zc - Haloes[k].Zc,2) );

				if(R > Rmin && R < Rmax)
				{
					A = 
						Haloes[j].Eax*Haloes[k].Eax +
						Haloes[j].Eay*Haloes[k].Eay + 
						Haloes[j].Eaz*Haloes[k].Eaz;
					B = (
						Haloes[j].Eax*(Haloes[j].Xc - Haloes[k].Xc) +
						Haloes[j].Eay*(Haloes[j].Yc - Haloes[k].Yc) + 
						Haloes[j].Eaz*(Haloes[j].Zc - Haloes[k].Zc) ) / R;

				for(i=0; i<nBins-1; i++) 
				{
					if(R>radius[i] && R<radius[i+1])
					{
						Abins[i] += sqrt(A*A);
						Bbins[i] += sqrt(B*B);
						Nbins[i] ++;
					}

				HaloZ.R[i]=Rbins[i]; 
				HaloZ.Th_c[i]=Abins[i]; 
				HaloZ.Th_p[i]=Bbins[i]; 
				HaloZ.N_pairs[i]=Nbins[i];

				}
			}
		}
	}

	free(radius);
	free(Abins);
	free(Bbins);
	free(Rbins);
	free(Nbins);

	fprintf(stdout, "\n");
}



void sort_shape_and_triaxiality()
{
	int i=0, m=0, nBins, nHaloesCut, nHaloes;
	int *array_shape_bin_y, *array_triax_bin_y;
	double *array_shape, *array_triax, *array_shape_bin, *array_triax_bin; 
	double half_t, half_s, sMax, sMin, tMax, tMin, p_s, p_t;

	fprintf(stdout, "\nSorting shape and triaxiality."); 

	nBins=Settings.n_bins;
	nHaloesCut=n_haloes_per_criterion();

#ifdef WITH_MPI
// The loop has to be done using ALL the haloes, since they are not ordered
		nHaloes=Settings.n_haloes; 
#else
		nHaloes=Settings.n_threshold; 
#endif

	Settings.tick=0;
	array_shape = (double*) calloc(nHaloesCut, sizeof(double));	
	array_triax = (double*) calloc(nHaloesCut, sizeof(double));	
	array_shape_bin = (double*) calloc(nBins, sizeof(double));	
	array_triax_bin = (double*) calloc(nBins, sizeof(double));	
	array_shape_bin_y = (int*) calloc(nBins-1, sizeof(int));	
	array_triax_bin_y = (int*) calloc(nBins-1, sizeof(int));	

		for(i=0; i<nHaloes; i++)
		{
			if(halo_condition(i) == 1)
			{
				array_shape[m] = Haloes[i].shape;
				array_triax[m] = Haloes[i].triax;
				m++;
			}
		}

			HaloZ.s0 = average(array_shape, nHaloesCut);
			HaloZ.t0 = average(array_triax, nHaloesCut);

				sMax = maximum(array_shape, nHaloesCut); 
				sMin = minimum(array_shape, nHaloesCut);
				tMax = maximum(array_triax, nHaloesCut); 
				tMin = minimum(array_triax, nHaloesCut);

			array_shape_bin = lin_stepper(sMin, sMax, nBins);
			lin_bin(array_shape, array_shape_bin, nBins, 
				nHaloesCut, array_shape_bin_y);	

			array_triax_bin = lin_stepper(tMin, tMax, nBins);
			lin_bin(array_triax, array_triax_bin, nBins, 
				nHaloesCut, array_triax_bin_y);	

			half_s = 0.5*(array_shape_bin[1]-array_shape_bin[0]);
			half_t = 0.5*(array_triax_bin[1]-array_triax_bin[0]);

			for(i=0; i<nBins-1; i++)	
			{
				p_s = (double) array_shape_bin_y[i]/nHaloesCut;
				p_t = (double) array_triax_bin_y[i]/nHaloesCut;
				HaloZ.shape[i]=array_shape_bin[i]+half_s;
				HaloZ.triax[i]=array_triax_bin[i]+half_t;
				HaloZ.p_shape[i]=p_s;
				HaloZ.p_triax[i]=p_t;
				HaloZ.n_shape[i]=array_shape_bin_y[i];
				HaloZ.n_triax[i]=array_triax_bin_y[i];
			}

	free(array_shape); 
	free(array_triax);
	free(array_shape_bin); 
	free(array_triax_bin);
	free(array_shape_bin_y); 
	free(array_triax_bin_y);
	fprintf(stdout, "\n");
}



void sort_radial_velocity()
{
	int i=0, m=0, nBins, nHaloesCut, nHaloes; 
	double mMax, mMin;
	double *radial_velocity_bin, *radial_velocity_error, *mass, *radial_velocity, *mass_bin; 

	fprintf(stdout, "\nSorting halo radial velocities.\n");

	nBins=Settings.n_bins;
	nHaloesCut=n_haloes_per_criterion();
	
#ifdef WITH_MPI
		nHaloes=Settings.n_haloes; 
#else
		nHaloes=Settings.n_threshold; 
#endif

		mass = (double*) calloc(nHaloesCut, sizeof(double));	
		radial_velocity = (double*) calloc(nHaloesCut, sizeof(double));	

		mass_bin = (double*) calloc(nBins, sizeof(double));	
		radial_velocity_bin = (double*) calloc(nBins-1, sizeof(double));	
		radial_velocity_error = (double*) calloc(nBins-1, sizeof(double));	// TODO error estimation

		for(i=0; i<nHaloes; i++)
		{
			if(halo_condition(i) == 1)
			{
				radial_velocity[m] = Haloes[i].Vmax;
				mass[m] = Haloes[i].Mvir; 
				m++;
			}
		}

			mMax = maximum(mass, nHaloesCut);
			mMin = minimum(mass, nHaloesCut);

			mass_bin = log_stepper(mMin, mMax, nBins);
			average_bin(mass, radial_velocity, mass_bin, radial_velocity_bin, 
				radial_velocity_error, nBins, nHaloesCut);

			for(i=0; i<nBins; i++)	
			{
				HaloZ.radVel[i]=radial_velocity_bin[i];
				HaloZ.err_radVel[i]=radial_velocity_error[i];
			}
	
	free(radial_velocity_error); 
	free(radial_velocity_bin); 
	free(radial_velocity); 
	free(mass_bin);
	free(mass);
	fprintf(stdout, "\n");
}



void sort_lambda()
{
	int nBins, nHaloesCut, nHaloes, i=0, m=0; 
	int *lambda_int_y; 
	double *bin_x, *params, *lambda, *lambda_bin_x, *lambda_err_y, *lambda_double_y;
	double l_0, sig, halfstep, lMax, lMin, delta_l, norm, value;

	fprintf(stdout, "\nSorting spin parameter.\n"); 

	nBins = Settings.n_bins;
	nHaloesCut = n_haloes_per_criterion();

#ifdef WITH_MPI
// The loop has to be done using ALL the haloes, since they are not ordered
		nHaloes=Settings.n_haloes; 
#else
		nHaloes=Settings.n_threshold; 
#endif

	bin_x = (double*) calloc(nBins, sizeof(double));	
	lambda = (double*) calloc(nHaloesCut, sizeof(double));	

	lambda_int_y = (int*) calloc(nBins-1, sizeof(int));	
	lambda_bin_x = (double*) calloc(nBins-1, sizeof(double));	
	lambda_double_y = (double*) calloc(nBins-1, sizeof(double));	
	lambda_err_y = (double*) calloc(nBins-1, sizeof(double));	

	params = (double*) calloc(2, sizeof(double));

		for(i=0; i<nHaloes; i++)
		{
			if(halo_condition(i) == 1)
			{
				lambda[m] = Haloes[i].lambda;
				m++;
			}
		}
	
				lMax = 1.01*maximum(lambda, nHaloesCut);  
				lMin = minimum(lambda, nHaloesCut);
				delta_l = (lMax-lMin)/nBins; 
				norm = 1./(delta_l*nHaloesCut);
	
				bin_x = lin_stepper(lMin, lMax, nBins);
				lin_bin(lambda, bin_x, nBins, nHaloesCut, lambda_int_y);	

			halfstep=(bin_x[1]-bin_x[0])*0.5;

			for(i=0; i<nBins-1; i++)
			{	
				value = (double) lambda_int_y[i];
				lambda_bin_x[i]=bin_x[i]+halfstep;
				lambda_err_y[i]=sqrt(value*norm); 
				lambda_double_y[i]=norm*value; 
			}

			params = best_fit_lognorm(lambda, nHaloesCut, nBins-1, 
					lambda_bin_x, lambda_double_y, lambda_err_y);

			l_0 = params[0]; 
			sig = params[1];
			HaloZ.l_0=l_0;
			HaloZ.l_sig=sig;

		for(i=0; i<nBins-1; i++)
		{		
			HaloZ.l[i]=lambda_bin_x[i];
			HaloZ.p_l[i]=lambda_double_y[i];
			HaloZ.err_p_l[i]=lambda_err_y[i];
		}	

	free(bin_x);
	free(lambda);
	free(lambda_int_y);
	free(lambda_double_y);
	free(lambda_err_y);
	free(lambda_bin_x);
	free(params);
	fprintf(stdout, "\n");
}



void sort_concentration()
{
	int nBins, nHaloes, nHaloesCut, i=0, m=0, *int_c_bin_y;
	double mMax, mMin, *c_err_mass, *params, *params2, *params3;
	double *conc, *mass, *c_avg_mass, *mass_bin, *bin_x, *c_bin_x, *c_bin_y, *c_err_y;
	double value, c_0, halfstep, c_02, sig2, cMax, cMin, norm;
	
	//double max, sig, M_0;

	fprintf(stdout, "\nSorting concentrations.\n");

	nBins = Settings.n_bins;
	nHaloesCut = n_haloes_per_criterion();

#ifdef WITH_MPI
// The loop has to be done using ALL the Haloes, since they are not ordered
		nHaloes=Settings.n_haloes; 
#else
		nHaloes=Settings.n_threshold; 
#endif

	conc = (double*) calloc(nHaloesCut, sizeof(double));	
	mass = (double*) calloc(nHaloesCut, sizeof(double));	

	bin_x = (double*) calloc(nBins, sizeof(double));	
	c_bin_x = (double*) calloc(nBins-1, sizeof(double));	
	c_bin_y = (double*) calloc(nBins-1, sizeof(double));	
	c_err_y = (double*) calloc(nBins-1, sizeof(double));	
	int_c_bin_y = (int*) calloc(nBins-1, sizeof(int));	
	mass_bin = (double*) calloc(nBins, sizeof(double));	
	c_avg_mass = (double*) calloc(nBins-1, sizeof(double));	
	c_err_mass = (double*) calloc(nBins-1, sizeof(double));	

	params = (double*) calloc(2, sizeof(double));
	params2 = (double*) calloc(2, sizeof(double));
	params3 = (double*) calloc(2, sizeof(double));

		for(i=0; i<nHaloes; i++) 
		{
			if(halo_condition(i) == 1)
			{ 
				conc[m] = Haloes[i].c_nfw;
	
					if(conc[m] == -1) 
						conc[m] = Haloes[i].c;

				mass[m] = Haloes[i].Mvir;
			//	fprintf(stdout, "mass: %e, conc: %f\n", mass[m], conc[m]);
				m++;
			}
		}
		
					mMax = 1.01*maximum(mass, nHaloesCut); 
					mMin = minimum(mass, nHaloesCut);
					cMax = 1.01*maximum(conc, nHaloesCut); 
					cMin = minimum(conc, nHaloesCut); 
 		
					bin_x=lin_stepper(cMin,cMax,nBins);
					lin_bin(conc,bin_x,nBins,nHaloesCut,int_c_bin_y);	

				norm=(nBins)/((cMax-cMin)*nHaloesCut);
				halfstep=0.5*(bin_x[1]-bin_x[0]);

			for(i=0; i<nBins-1; i++) 
			{
				value = (double) int_c_bin_y[i];
				c_bin_x[i] = bin_x[i] + halfstep;
				c_bin_y[i] = norm*value;
				c_err_y[i] = sqrt(norm*value);
			}	

			//c_0 = average(conc, nHaloes_conc);
			c_0 = maximum(c_bin_y, nBins-1);
			/*

				params = best_fit_lognorm(conc, nHaloes_conc, nBins-1, 
					c_bin_x, c_bin_y, c_err_y);

				c_0 = params[0]; 
				sig = params[1]; 
				max=maximum(c_bin_y, nBins-1);
			
				HaloZ.c_0=c_0;
				HaloZ.c_sig=sig;
			*/

			for(i=0; i<nBins-1; i++)
			{	
				HaloZ.c[i]=c_bin_x[i];
				HaloZ.p_c[i]=c_bin_y[i];
			}

			for(i=0; i<nHaloesCut; i++) 
				conc[i]/=c_0;

		for(i=0; i<nBins-1; i++)
		{	
			c_bin_x[i]/=c_0;
			c_bin_y[i]/=c_0;
		}

		params2 = best_fit_lognorm(conc, nHaloesCut, nBins-1, 
			c_bin_x, c_bin_y, c_err_y);

		c_02 = params2[0]; 
		sig2 = params2[1];

			for(i=0; i<nHaloesCut; i++) 
				conc[i]*=c_0;

			mass_bin = log_stepper(mMin, mMax, nBins);
			average_bin(mass, conc, mass_bin, c_avg_mass,
				c_err_mass, nBins, nHaloesCut);

			for(i=0; i<nBins-1; i++)
			{
				HaloZ.c_avg[i]=c_avg_mass[i];
			}
	
/*	
			params3[0] = 1.5; 
			params3[1] = pow(10.e-14, params[0]);
			params3=best_fit_power_law(mass_bin, c_avg_mass, c_err_y, nBins-1, params3);

		M_0 = -log(params3[1])/params3[0];
		fprintf(stdout, "M-Conc    a:%lf   M_0 10e+%f SM\n",params3[0],M_0/log(10));
*/
	free(conc); 
	free(mass); 
	free(c_avg_mass); 
	free(c_err_mass); 
	free(int_c_bin_y);
	free(bin_x); 
	free(c_bin_x); 
	free(c_bin_y); 
	free(c_err_y);
	free(params);
	free(params2);
	free(params3);
	fprintf(stdout, "\n");
}


#ifdef GAS
void sort_gas_fraction()
{
	int nBins, nHaloes, nHaloesCut, i=0, m=0;
	double *gas_fraction, *gas_fraction_bin, *gas_fraction_error, *mass_bin, *mass;
	double mMax, mMin;
	
	fprintf(stdout, "\nSorting halo gas fraction.\n");

	Settings.tick=0;
	mass = (double*) calloc(nHaloes, sizeof(double));	
	mass_bin = (double*) calloc(nBins, sizeof(double));	
	gas_fraction = (double*) calloc(nHaloes, sizeof(double));	
	gas_fraction_bin = (double*) calloc(nBins-1, sizeof(double));	
	gas_fraction_error = (double*) calloc(nBins-1, sizeof(double));	

		for(i=0; i<nHaloes; i++)
		{
			mass[i] = Haloes[i].Mvir;
			gas_fraction[i] = Haloes[i].b_fraction;
		}
	
			mass_bin = log_stepper(mMin, mMax, nBins);
			average_bin(mass, gas_fraction, mass_bin, gas_fraction_bin,  
					gas_fraction_error, nBins, nHaloes);

		for(i=0; i<nBins-1; i++)
		{
			HaloZ.mass[i]=0.5*(mass_bin[i]+mass_bin[i+1]);
			HaloZ.gas_fraction[i]=gas_fraction_bin[i];
		}

	free(mass);
	free(mass_bin);
	free(gas_fraction); 
	free(gas_fraction_bin); 
	free(gas_fraction_error); 
	fprintf(stdout, "\n");
}



void sort_and_fit_mass_temperature_relation()
{
 	int nBins=Settings.n_bins, nHaloes=Settings.n_threshold, i=0, m=0, *n_per_mass_bin;
	double *a, *temperature, *temperature_bin, *temperature_error, *mass_bin, *mass;
	double M_0, mMax, mMin;
	
	fprintf(stdout, "\nSorting and fitting mass temperature relation.\n");
	
	Settings.tick=0;
	a = (double*) calloc(2,sizeof(double));
	mass = (double*) calloc(nHaloes, sizeof(double));	
	mass_bin = (double*) calloc(nBins, sizeof(double));	
	temperature = (double*) calloc(nHaloes, sizeof(double));	
	temperature_bin = (double*) calloc(nBins-1, sizeof(double));	
	temperature_error = (double*) calloc(nBins-1, sizeof(double));	

		for(i=0; i<nHaloes; i++)
		{
			mass[i] = Haloes[i].Mvir;
			temperature[i] = Haloes[i].T_gas;
		}

			mass_bin = log_stepper(mMin, mMax, nBins);
			average_bin(mass, temperature, mass_bin, temperature_bin
					temperature_error, nBins, nHaloes);

			for(i=0; i<nBins-1; i++)
			{
				HaloZ.gas_T[i]=temperature_bin[i];
				HaloZ.mass[i]=0.5*(mass_bin[i]+mass_bin[i+1]);
				temperature_error[i] = HaloZ.gas_T[i]/sqrt(MF.n_bin[i]);
			}
	
			a[0] = 1.5; 
			a[1] = pow(10.e-14, a[0]);
	
			a=best_fit_power_law(mass_bin, temperature_bin, temperature_error, nBins-1, a);

		M_0 = -log(a[1])/a[0];
		fprintf(stdout, "M-Tx    a:%lf   M_0 10e+%f SM\n",a[0],M_0/log(10));

	free(a);
	free(mass);
	free(mass_bin);
	free(temperature); 
	free(temperature_bin); 
	free(temperature_error); 
	fprintf(stdout, "\n");
}
#endif // Gas



void compute_halo_properties()
{
		avg_subhalo();

		compute_numerical_mass_function();

			sort_axis_alignement();
			sort_shape_and_triaxiality();
			sort_lambda();
			sort_concentration();

#ifdef GAS
		sort_gas_fraction();
		sort_and_fit_mass_temperature_relation();
#endif
}
