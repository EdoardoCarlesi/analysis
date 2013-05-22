#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../libmath/math.h"
#include "../libcosmo/cosmo.h"
#include "../general_def.h"

#include "halo.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef WITH_MPI
#include "../libparallel/general.h"
#endif

/*
 * Declare functions
 */
void sort_lambda_and_concentration(void);
void sort_shape_and_triaxiality(void);
void sort_nfw_parameters(void);
void sort_mass_relations(void);

#ifdef GAS
void sort_T_mass_function(void);
void sort_gas_relations(void);
void sort_alignment_and_displacement(void);
#endif


/*
 * Initialize functions
 */ 
void sort_axis_alignment()
{
	int m=0, i=0, j=0, k=0, p=0, n=0, max_haloes, nBins, skip;
	int *Nbins, *index;
	double *radius, *Rbins, *Abins, *Bbins; 
	double Rmin, Rmax, R, A, B, sum;

	INFO_MSG("Computing halo major axis alignement angles");

#ifdef USE_SAMPLE	
		max_haloes = Settings.n_threshold; 
#else
		max_haloes = Settings.n_haloes;
#endif
		skip = Settings.halo_skip; 
		nBins = Settings.r_bins; 

		radius = (double *) calloc(nBins, sizeof(double));
		Rbins = (double *) calloc(nBins-1, sizeof(double));
		Abins = (double *) calloc(nBins-1, sizeof(double));
		Bbins = (double *) calloc(nBins-1, sizeof(double));
		Nbins = (int *) calloc(nBins-1, sizeof(int));
		
		index = (int *) calloc(max_haloes, sizeof(int));
	
		list_halo_sample(index);

		Rmin = Settings.Rmin; Rmax = Settings.Rmax;
		radius = log_stepper(Rmin,Rmax,nBins);

			for(i=0; i<nBins-1; i++)
				Rbins[i] = 0.5*(radius[i+1]+radius[i]);

#ifdef USE_SAMPLE
#		pragma omp parallel for			\
		private(m, i, j, k, A, B, R, sum) shared(Haloes, Rmin, Rmax)
			for(j=0; j<max_haloes; j++) 
			{
				p = index[j];
				for(k=j; k<max_haloes; k++)
				{
					n = index[k];
					A = 0; B = 0; R = 0; sum = 0;
					for(m=0; m<3; m++)
						sum += pow2(Haloes[p].X[m] - Haloes[n].X[m]);
	
					R = sqrt(sum);

					if(R > Rmin && R < Rmax)
					{
						for(m=0; m<3; m++)
							A += Haloes[p].Ea[m]*Haloes[n].Ea[m];

						for(m=0; m<3; m++)
							B += Haloes[p].Ea[m]*(Haloes[p].X[m] - Haloes[n].X[m]);

							B /= R;

				for(i=0; i<nBins-1; i++) 
				{
					if(R>radius[i] && R<radius[i+1])
					{
						Abins[i] += sqrt(A*A);
						Bbins[i] += sqrt(B*B);
						Nbins[i]++;
					}

					HaloProperties[HALO_INDEX].R[i]=Rbins[i]; 
					HaloProperties[HALO_INDEX].Th_c[i]=Abins[i]; 
					HaloProperties[HALO_INDEX].Th_p[i]=Bbins[i]; 
					HaloProperties[HALO_INDEX].N_pairs[i]=Nbins[i]; 

				}
			}
		}
	}
#else
#		pragma omp parallel for			\
		private(m, i, j, k, A, B, R, sum) shared(Haloes, Rmin, Rmax)
			for(j=0; j<max_haloes; j++) 
			{
				for(k=j; k<max_haloes; k++)
				{
					A = 0; B = 0; R = 0; sum = 0;
				// Use haloes above a given threshold 
				if(Haloes[j].Mvir > Settings.mass_min && Haloes[k].Mvir > Settings.mass_min)
				{	
					for(m=0; m<3; m++)
						sum += pow2(Haloes[j].X[m] - Haloes[k].X[m]);
	
					R = sqrt(sum);

					if(R > Rmin && R < Rmax)
					{
						for(m=0; m<3; m++)
							A += Haloes[j].Ea[m]*Haloes[k].Ea[m];

						for(m=0; m<3; m++)
							B += Haloes[j].Ea[m]*(Haloes[j].X[m] - Haloes[k].X[m]);

							B /= R;

				for(i=0; i<nBins-1; i++) 
				{
					if(R>radius[i] && R<radius[i+1])
					{
						Abins[i] += sqrt(A*A);
						Bbins[i] += sqrt(B*B);
						Nbins[i]++;
					}

					HaloProperties[HALO_INDEX].R[i]=Rbins[i]; 
					HaloProperties[HALO_INDEX].Th_c[i]=Abins[i]; 
					HaloProperties[HALO_INDEX].Th_p[i]=Bbins[i]; 
					HaloProperties[HALO_INDEX].N_pairs[i]=Nbins[i]; 

					}
				}
			} // End if mass cut
		}
	}
#endif

	free(radius);
	free(Abins);
	free(Bbins);
	free(Rbins);
	free(Nbins);
}



void sort_numerical_mass_function(void)
{
	int nBins=0, nHaloes=0, nHaloesCut=0, i=0, j=0, k=0, l=0; 
	double dn_norm=1., volume=0, mMin=0, mMax=0, vMax=0, vMin=0, halfstep=0, dM=0, dV=0; 
	double gMax, gMin, ngMax, ngMin;
	double *mass, *mass_bin; 
	int *n_mass, *cum_n_mass;
	double *vel, *vel_bin; 
	int *n_vel, *cum_n_vel;
	double *gas, *gas_bin; 
	int *n_gas, *cum_n_gas;
	double *nogas, *nogas_bin; 
	int *n_nogas, *cum_n_nogas;
	double totMass, totDm, totGas;

	nBins = Settings.n_bins; 

	nBins = Settings.n_bins; 
	nHaloes = Settings.n_haloes; 

	Settings.totHaloMass = 0; 
	Settings.totDarkMass = 0;
	Settings.totGasMassInHalo = 0;

	if(Settings.use_web == 1)
	{
		nHaloesCut = Settings.n_cweb_type[Settings.use_web_type]; 
		D_PRINT("Sorting mass function for halo number=", nHaloesCut);
	}
		else if(Settings.use_sub == 1)
	{ 
		nHaloesCut = SubStructure.N_sub;
		D_PRINT("Sorting mass function for subhalo number=", nHaloesCut);
	}
		else
	{
		nHaloesCut = Settings.n_haloes; 
		D_PRINT("Sorting mass function for halo number=", nHaloesCut);
	}
	
		Settings.tick=0;
	
		mass = (double*) calloc(nHaloesCut, sizeof(double));
		mass_bin = (double*) calloc(nBins, sizeof(double));
		n_mass = (int*) calloc(nBins-1, sizeof(int));
		cum_n_mass = (int*) calloc(nBins-1, sizeof(int));

		vel = (double*) calloc(nHaloesCut, sizeof(double));
		vel_bin = (double*) calloc(nBins, sizeof(double));
		n_vel = (int*) calloc(nBins-1, sizeof(int));
		cum_n_vel = (int*) calloc(nBins-1, sizeof(int));

		gas = (double*) calloc(nHaloesCut, sizeof(double));
		gas_bin = (double*) calloc(nBins, sizeof(double));
		n_gas = (int*) calloc(nBins-1, sizeof(int));
		cum_n_gas = (int*) calloc(nBins-1, sizeof(int));

		nogas = (double*) calloc(nHaloesCut, sizeof(double));
		nogas_bin = (double*) calloc(nBins, sizeof(double));
		n_nogas = (int*) calloc(nBins-1, sizeof(int));
		cum_n_nogas = (int*) calloc(nBins-1, sizeof(int));

		alloc_mass_function(&MassFunc[MF_INDEX], nBins);
		alloc_mass_function(&VelFunc[MF_INDEX], nBins);
		alloc_mass_function(&GasFunc[MF_INDEX], nBins);
		alloc_mass_function(&NoGasFunc[MF_INDEX], nBins);

		for(i=0; i<nHaloes; i++)
		{	
			if(Settings.use_sub == 1) 
			{
				if(Haloes[i].host > 0)
				{
					Settings.totHaloMass += Haloes[i].Mvir;
					mass[j] = Haloes[i].Mvir;
					vel[j]  = Haloes[i].Vmax;
					j++;

					if(Haloes[i].gas.M > 0)
					{	
						Settings.totGasMassInHalo += Haloes[i].gas.M;
						gas[k] = Haloes[i].gas.M;
						k++;
					} 
					else
					{
						Settings.totDarkMass += Haloes[i].Mvir;
						nogas[l] = Haloes[i].Mvir;
						l++;
					}
				}
			}
			else if (Settings.use_sub == 0) 
			{
				if(Settings.use_web == 1) 
				{
					if(Haloes[i].web_type[Settings.use_web_type] == 1)
					{
						Settings.totHaloMass += Haloes[i].Mvir;
						mass[j] = Haloes[i].Mvir;
						vel[j]  = Haloes[i].Vmax;
						j++;
					
						if(Haloes[i].gas.M > 0)
						{
							Settings.totGasMassInHalo += Haloes[i].gas.M;
							gas[k] = Haloes[i].gas.M;
							k++;
						} 
						else
						{
							Settings.totDarkMass += Haloes[i].Mvir;
							nogas[l] = Haloes[i].Mvir;
							l++;
						}
	
					}
				}
				else 
				{
			
					if(Haloes[i].gas.M > 0)
					{
						Settings.totGasMassInHalo += Haloes[i].gas.M;
						gas[k] = Haloes[i].gas.M;
						k++;
					} 
					else
					{
						Settings.totDarkMass += Haloes[i].Mvir;
						nogas[l] = Haloes[i].Mvir;
						l++;
					}

					Settings.totHaloMass += Haloes[i].Mvir;
					mass[i] = Haloes[i].Mvir;			
					vel[i]  = Haloes[i].Vmax;
				}
			}
		}

		fprintf(stderr, "\nBinned\n");

		totDm = Settings.dmMass * pow3(Settings.n_part_1D);		
		totGas = Settings.gasMass * pow3(Settings.n_part_1D);		
		totMass = totDm + totGas;

		Settings.totGasMassInHalo /= totGas;
		Settings.totHaloMass /= totMass;
		Settings.totDarkMass /= totMass;

		mMin = F_MIN * minimum(mass, nHaloesCut);
		mMax = F_MAX * maximum(mass, nHaloesCut);
		vMin = F_MIN * minimum(vel, nHaloesCut);
		vMax = F_MAX * maximum(vel, nHaloesCut);
	
		gMin = F_MIN * minimum(gas, k);	
		gMax = F_MAX * maximum(gas, k);	
	
		ngMin = F_MIN * minimum(nogas, l);	
		ngMax = F_MAX * maximum(nogas, l);	

		mass_bin = log_stepper(mMin, mMax, nBins);
		vel_bin = log_stepper(vMin, vMax, nBins);
		gas_bin = log_stepper(gMin, gMax, nBins);
		nogas_bin = log_stepper(ngMin, ngMax, nBins);
	
		lin_bin(mass, mass_bin, nBins, nHaloesCut, n_mass);	
		cum_bin(n_mass, cum_n_mass, nBins-1);
	
		lin_bin(vel, vel_bin, nBins, nHaloesCut, n_vel);	
		cum_bin(n_vel, cum_n_vel, nBins-1);
	
		lin_bin(gas, gas_bin, nBins, k, n_gas);	
		cum_bin(n_gas, cum_n_gas, nBins-1);
	
		lin_bin(nogas, nogas_bin, nBins, l, n_nogas);	
		cum_bin(n_nogas, cum_n_nogas, nBins-1);

		volume=Settings.box_size*Settings.box_size*Settings.box_size;


		for(i=0; i<nBins-1; i++)
		{
			dn_norm = 2*halfstep/nHaloesCut;
			dM = mass_bin[i+1]-mass_bin[i];
			dV = vel_bin[i+1]-vel_bin[i];

			MassFunc[MF_INDEX].mass[i]=mass_bin[i];
			MassFunc[MF_INDEX].mass_halfstep[i]=0.5 * (mass_bin[i]+mass_bin[i+1]);
			MassFunc[MF_INDEX].dn[i]=n_mass[i]/(volume*dM);
			MassFunc[MF_INDEX].n[i]=cum_n_mass[i]/volume;
			MassFunc[MF_INDEX].err_dn[i]=sqrt(n_mass[i])/(volume*dM);
			MassFunc[MF_INDEX].err[i]=sqrt(cum_n_mass[i])/volume;

			MassFunc[MF_INDEX].n_bin[i]=n_mass[i];
			MassFunc[MF_INDEX].n_tot[i]=cum_n_mass[i];

			VelFunc[MF_INDEX].mass[i] = 0.5 * (vel_bin[i] + vel_bin[i+1]);
			VelFunc[MF_INDEX].n[i] = cum_n_vel[i]/volume;
			VelFunc[MF_INDEX].dn[i] = n_vel[i]/volume/dV;
			VelFunc[MF_INDEX].n_bin[i] = n_vel[i];
			VelFunc[MF_INDEX].n_tot[i] = cum_n_vel[i];
	
			GasFunc[MF_INDEX].mass[i] = 0.5 * (gas_bin[i] + gas_bin[i+1]);
			GasFunc[MF_INDEX].n[i] = cum_n_gas[i]/volume;
			GasFunc[MF_INDEX].n_tot[i] = cum_n_gas[i];

			NoGasFunc[MF_INDEX].mass[i] = 0.5 * (nogas_bin[i] + nogas_bin[i+1]);
			NoGasFunc[MF_INDEX].n[i] = cum_n_nogas[i]/volume;
			NoGasFunc[MF_INDEX].n_tot[i] = cum_n_nogas[i];
		}
	
	free(cum_n_mass);
	free(mass_bin); 
	free(n_mass); 
	free(mass); 
	free(cum_n_vel);
	free(vel_bin); 
	free(n_vel); 
	free(vel); 
	free(cum_n_gas);
	free(gas_bin); 
	free(n_gas); 
	free(gas); 
	free(cum_n_nogas);
	free(nogas_bin); 
	free(n_nogas); 
	free(nogas); 
}



void sort_shape_and_triaxiality()
{
	int i=0, m=0, nBins, nHaloesCut, nHaloes;
	int *array_shape_bin_y, *array_triax_bin_y;
	double *array_shape, *array_triax, *array_shape_bin, *array_triax_bin; 
	double half_t, half_s, sMax, sMin, tMax, tMin, p_s, p_t;

#ifdef GAS
	int *array_gas_shape_bin_y, *array_gas_triax_bin_y;
	double *array_gas_shape, *array_gas_triax, *array_gas_shape_bin, *array_gas_triax_bin; 
	int *array_diff_shape_bin_y, *array_diff_triax_bin_y;
	double *array_diff_shape, *array_diff_triax, *array_diff_shape_bin, *array_diff_triax_bin; 
	double half_gt, half_gs, gsMax, gsMin, gtMax, gtMin, p_gs, p_gt;
	double half_dt, half_ds, dsMax, dsMin, dtMax, dtMin, p_ds, p_dt;
#endif

	INFO_MSG("Sorting shape and triaxiality");

	nBins = Settings.n_bins;
	nHaloesCut = n_haloes_per_criterion();

#ifdef WITH_MPI
		nHaloes = Settings.n_haloes; 
#else
		nHaloes = Settings.n_threshold; 
#endif

	Settings.tick=0;
	array_shape = (double*) calloc(nHaloesCut, sizeof(double));	
	array_triax = (double*) calloc(nHaloesCut, sizeof(double));	
	array_shape_bin = (double*) calloc(nBins, sizeof(double));	
	array_triax_bin = (double*) calloc(nBins, sizeof(double));	
	array_shape_bin_y = (int*) calloc(nBins-1, sizeof(int));	
	array_triax_bin_y = (int*) calloc(nBins-1, sizeof(int));	
#ifdef GAS
	array_diff_shape = (double*) calloc(nHaloesCut, sizeof(double));	
	array_diff_triax = (double*) calloc(nHaloesCut, sizeof(double));	
	array_diff_shape_bin = (double*) calloc(nBins, sizeof(double));	
	array_diff_triax_bin = (double*) calloc(nBins, sizeof(double));	
	array_diff_shape_bin_y = (int*) calloc(nBins-1, sizeof(int));	
	array_diff_triax_bin_y = (int*) calloc(nBins-1, sizeof(int));	

	array_gas_shape = (double*) calloc(nHaloesCut, sizeof(double));	
	array_gas_triax = (double*) calloc(nHaloesCut, sizeof(double));	
	array_gas_shape_bin = (double*) calloc(nBins, sizeof(double));	
	array_gas_triax_bin = (double*) calloc(nBins, sizeof(double));	
	array_gas_shape_bin_y = (int*) calloc(nBins-1, sizeof(int));	
	array_gas_triax_bin_y = (int*) calloc(nBins-1, sizeof(int));	
#endif
		for(i=0; i<nHaloes; i++)
		{
			if(halo_condition(i) == 1)
			{
				array_shape[m] = Haloes[i].shape;
				array_triax[m] = Haloes[i].triax;
#ifdef GAS
				array_gas_shape[m] = Haloes[i].gas_only.shape;
				array_gas_triax[m] = Haloes[i].gas_only.triax;
				array_diff_shape[m] = Haloes[i].gas_only.diff.shape;
				array_diff_triax[m] = Haloes[i].gas_only.diff.triax;
#endif
				m++;
			}
		}

			HaloProperties[HALO_INDEX].halo.s0 = average(array_shape, nHaloesCut);
			HaloProperties[HALO_INDEX].halo.t0 = average(array_triax, nHaloesCut);
#ifdef GAS
			HaloProperties[HALO_INDEX].diff.s0 = average(array_diff_shape, nHaloesCut);
			HaloProperties[HALO_INDEX].diff.t0 = average(array_diff_triax, nHaloesCut);

			HaloProperties[HALO_INDEX].gas.s0 = average(array_gas_shape, nHaloesCut);
			HaloProperties[HALO_INDEX].gas.t0 = average(array_gas_triax, nHaloesCut);
#endif
			sMax = F_MAX*maximum(array_shape, nHaloesCut); 
			sMin = F_MIN*minimum(array_shape, nHaloesCut);
			tMax = F_MAX*maximum(array_triax, nHaloesCut); 
			tMin = F_MIN*minimum(array_triax, nHaloesCut);
#ifdef GAS
			dsMax = F_MAX*maximum(array_diff_shape, nHaloesCut); 
			dsMin = F_MIN*minimum(array_diff_shape, nHaloesCut);
			dtMax = F_MAX*maximum(array_diff_triax, nHaloesCut); 
			dtMin = F_MIN*minimum(array_diff_triax, nHaloesCut);
	
			gsMax = F_MAX*maximum(array_gas_shape, nHaloesCut); 
			gsMin = F_MIN*minimum(array_gas_shape, nHaloesCut);
			gtMax = F_MAX*maximum(array_gas_triax, nHaloesCut); 
			gtMin = F_MIN*minimum(array_gas_triax, nHaloesCut);
#endif
			array_shape_bin = lin_stepper(sMin, sMax, nBins);
			lin_bin(array_shape, array_shape_bin, nBins, nHaloesCut, array_shape_bin_y);	

			array_triax_bin = lin_stepper(tMin, tMax, nBins);
			lin_bin(array_triax, array_triax_bin, nBins, nHaloesCut, array_triax_bin_y);	

			half_s = 0.5*(array_shape_bin[1]-array_shape_bin[0]);
			half_t = 0.5*(array_triax_bin[1]-array_triax_bin[0]);

#ifdef GAS
			array_gas_shape_bin = lin_stepper(gsMin, gsMax, nBins);
			lin_bin(array_gas_shape, array_gas_shape_bin, nBins, nHaloesCut, array_gas_shape_bin_y);	

			array_gas_triax_bin = lin_stepper(gtMin, gtMax, nBins);
			lin_bin(array_gas_triax, array_gas_triax_bin, nBins, nHaloesCut, array_gas_triax_bin_y);	

			array_diff_shape_bin = lin_stepper(dsMin, dsMax, nBins);
			lin_bin(array_diff_shape, array_diff_shape_bin, nBins, nHaloesCut, array_diff_shape_bin_y);	

			array_diff_triax_bin = lin_stepper(dtMin, dtMax, nBins);
			lin_bin(array_diff_triax, array_diff_triax_bin, nBins, nHaloesCut, array_diff_triax_bin_y);	

			half_ds = 0.5*(array_diff_shape_bin[1]-array_diff_shape_bin[0]);
			half_dt = 0.5*(array_diff_triax_bin[1]-array_diff_triax_bin[0]);

			half_gs = 0.5*(array_gas_shape_bin[1]-array_gas_shape_bin[0]);
			half_gt = 0.5*(array_gas_triax_bin[1]-array_gas_triax_bin[0]);
#endif

			for(i=0; i<nBins-1; i++)	
			{
				p_s = (double) array_shape_bin_y[i]/nHaloesCut;
				p_t = (double) array_triax_bin_y[i]/nHaloesCut;
				HaloProperties[HALO_INDEX].halo.s[i]=array_shape_bin[i]+half_s;
				HaloProperties[HALO_INDEX].halo.t[i]=array_triax_bin[i]+half_t;
				HaloProperties[HALO_INDEX].halo.p_s[i]=p_s;
				HaloProperties[HALO_INDEX].halo.p_t[i]=p_t;
#ifdef GAS
				p_ds = (double) array_diff_shape_bin_y[i]/nHaloesCut;
				p_dt = (double) array_diff_triax_bin_y[i]/nHaloesCut;
				HaloProperties[HALO_INDEX].diff.s[i]=array_diff_shape_bin[i]+half_ds;
				HaloProperties[HALO_INDEX].diff.t[i]=array_diff_triax_bin[i]+half_dt;
				HaloProperties[HALO_INDEX].diff.p_s[i]=p_ds;
				HaloProperties[HALO_INDEX].diff.p_t[i]=p_dt;

				p_gs = (double) array_gas_shape_bin_y[i]/nHaloesCut;
				p_gt = (double) array_gas_triax_bin_y[i]/nHaloesCut;
				HaloProperties[HALO_INDEX].gas.s[i]=array_gas_shape_bin[i]+half_gs;
				HaloProperties[HALO_INDEX].gas.t[i]=array_gas_triax_bin[i]+half_gt;
				HaloProperties[HALO_INDEX].gas.p_s[i]=p_gs;
				HaloProperties[HALO_INDEX].gas.p_t[i]=p_gt;
#endif
			}

	free(array_shape); 
	free(array_triax);
	free(array_shape_bin); 
	free(array_triax_bin);
	free(array_shape_bin_y); 
	free(array_triax_bin_y);
#ifdef GAS
	free(array_gas_shape); 
	free(array_gas_triax);
	free(array_gas_shape_bin); 
	free(array_gas_triax_bin);
	free(array_gas_shape_bin_y); 
	free(array_gas_triax_bin_y);

	free(array_diff_shape); 
	free(array_diff_triax);
	free(array_diff_shape_bin); 
	free(array_diff_triax_bin);
	free(array_diff_shape_bin_y); 
	free(array_diff_triax_bin_y);
#endif
}



void sort_lambda_and_concentration()
{
	int nBins, nHaloesCut, nHaloes, i=0, m=0; 
	int *lambda_int_y; 
	double *l_bin_x, *params, *lambda, *lambda_bin_x, *lambda_err_y, *lambda_double_y;
	double l_0, sig, l_halfstep, lMax, lMin, delta_l, l_norm, l_value;
	int *conc_int_y; 
	double *c_bin_x, *conc, *conc_bin_x, *conc_err_y, *conc_double_y;
	double c_halfstep, cMax, cMin, c_norm, c_value;
	int *avg_sub_int_y; 
	double *a_bin_x, *avg_sub, *avg_sub_bin_x, *avg_sub_err_y, *avg_sub_double_y;
	double a_halfstep, aMax, aMin, a_norm, a_value;
#ifdef GAS
	// TODO add gas lambda only
	//int *lambda_int_y; 
	//double *l_bin_x, *params, *lambda, *lambda_bin_x, *lambda_err_y, *lambda_double_y;
	//double l_0, sig, l_halfstep, lMax, lMin, delta_l, l_norm, l_value;
#endif

	INFO_MSG("Sorting spin parameter"); 

	nBins = Settings.n_bins;
	nHaloesCut = n_haloes_per_criterion();

		HaloProperties = HaloProperties;


#ifdef WITH_MPI
		nHaloes = Settings.n_haloes; 
#else
		nHaloes = Settings.n_threshold; 
#endif

	l_bin_x = (double*) calloc(nBins, sizeof(double));	
	c_bin_x = (double*) calloc(nBins, sizeof(double));	
	a_bin_x = (double*) calloc(nBins, sizeof(double));	

	lambda_int_y = (int*) calloc(nBins-1, sizeof(int));	
	lambda_bin_x = (double*) calloc(nBins-1, sizeof(double));	
	lambda_double_y = (double*) calloc(nBins-1, sizeof(double));	
	lambda_err_y = (double*) calloc(nBins-1, sizeof(double));	

	conc_int_y = (int*) calloc(nBins-1, sizeof(int));	
	conc_bin_x = (double*) calloc(nBins-1, sizeof(double));	
	conc_double_y = (double*) calloc(nBins-1, sizeof(double));	
	conc_err_y = (double*) calloc(nBins-1, sizeof(double));	

	avg_sub_int_y = (int*) calloc(nBins-1, sizeof(int));	
	avg_sub_bin_x = (double*) calloc(nBins-1, sizeof(double));	
	avg_sub_double_y = (double*) calloc(nBins-1, sizeof(double));	
	avg_sub_err_y = (double*) calloc(nBins-1, sizeof(double));	

	params = (double*) calloc(2, sizeof(double));

	lambda = (double*) calloc(nHaloesCut, sizeof(double));	
	conc = (double*) calloc(nHaloesCut, sizeof(double));	
	avg_sub = (double*) calloc(nHaloesCut, sizeof(double));	

		for(i=0; i<nHaloes; i++)
		{
			if(halo_condition(i) == 1)
			{
				lambda[m] = Haloes[i].lambda;
				//conc[m] = Haloes[i].c_nfw;
				conc[m] = Haloes[i].Rvir/Haloes[i].fit_nfw.rs;
				//conc[m] = Haloes[i].Rvir/Haloes[i].fit_nfw.c;
				avg_sub[m] = (double) Haloes[i].n_satellites;

				//if(conc[m] == -1) 
				//	conc[m] = Haloes[i].c;
				m++;
			}
		}

			lMax = F_MAX*maximum(lambda, nHaloesCut);  
			lMin = F_MIN*minimum(lambda, nHaloesCut);
			delta_l = (lMax-lMin)/nBins; 
			l_norm = 1./(delta_l*nHaloesCut);

			cMax = F_MAX*maximum(conc, nHaloesCut);  
			cMin = F_MIN*minimum(conc, nHaloesCut);
			c_norm = 1./(nHaloesCut);

			aMax = F_MAX*maximum(avg_sub, nHaloesCut);  
			aMin = F_MIN*minimum(avg_sub, nHaloesCut);
			a_norm = 1./(nHaloesCut);

			l_bin_x = lin_stepper(lMin, lMax, nBins);
			lin_bin(lambda, l_bin_x, nBins, nHaloesCut, lambda_int_y);	
			l_halfstep=(l_bin_x[1]-l_bin_x[0])*0.5;

			c_bin_x = lin_stepper(cMin, cMax, nBins);
			lin_bin(conc, c_bin_x, nBins, nHaloesCut, conc_int_y);	
			c_halfstep=(c_bin_x[1]-c_bin_x[0])*0.5;

			a_bin_x = lin_stepper(aMin, aMax, nBins);
			lin_bin(avg_sub, a_bin_x, nBins, nHaloesCut, avg_sub_int_y);	
			a_halfstep=(a_bin_x[1]-a_bin_x[0])*0.5;

			for(i=0; i<nBins-1; i++)
			{	
				l_value = (double) lambda_int_y[i];
				lambda_bin_x[i]=l_bin_x[i]+l_halfstep;
				lambda_err_y[i]=sqrt(l_value)*(l_norm); 
				lambda_double_y[i]=l_norm*l_value; 

				c_value = (double) conc_int_y[i];
				conc_bin_x[i]=c_bin_x[i]+c_halfstep;
				conc_err_y[i]=sqrt(c_value)*c_norm; 
				conc_double_y[i]=c_norm*c_value; 

				a_value = (double) avg_sub_int_y[i];
				avg_sub_bin_x[i]=a_bin_x[i]+a_halfstep;
				avg_sub_err_y[i]=sqrt(a_value)*a_norm; 
				avg_sub_double_y[i]=a_norm*a_value; 
			}

		for(i=0; i<nBins-1; i++)
		{		
			HaloProperties[HALO_INDEX].halo.l[i]=lambda_bin_x[i];
			HaloProperties[HALO_INDEX].halo.p_l[i]=lambda_double_y[i];
			HaloProperties[HALO_INDEX].c[i]=conc_bin_x[i];
			HaloProperties[HALO_INDEX].p_c[i]=conc_double_y[i];
			HaloProperties[HALO_INDEX].avg_sub[i]=conc_bin_x[i];
			HaloProperties[HALO_INDEX].p_avg_sub[i]=conc_double_y[i];
		}
	
			INFO_MSG("Fitting spin parameter distribution to a lognorm");
			params = best_fit_lognorm(lambda, nHaloesCut, nBins-1, 
				lambda_bin_x, lambda_double_y, lambda_err_y);

			l_0 = params[0]; 
			sig = params[1];
			HaloProperties[HALO_INDEX].l_0=l_0;
			HaloProperties[HALO_INDEX].l_sig=sig;

	free(conc_double_y);
	free(conc_int_y);
	free(conc_err_y);
	free(conc_bin_x);
	free(conc);
	free(c_bin_x);
	free(lambda_double_y);
	free(lambda_int_y);
	free(lambda_err_y);
	free(lambda_bin_x);
	free(lambda);
	free(l_bin_x);
	free(avg_sub_double_y);
	free(avg_sub_int_y);
	free(avg_sub_err_y);
	free(avg_sub_bin_x);
	free(avg_sub);
	free(a_bin_x);
	free(params);
}



void sort_nfw_parameters()
{
	int nBins, nHaloesCut, nHaloes, i=0, m=0; 
	int *nfw_gof_int_y; 
	double *gbin_x,*nfw_gof, *nfw_gof_bin_x, *nfw_gof_err_y, *nfw_gof_double_y;
	double ghalfstep, gMax, gMin, delta_g, gnorm, gvalue;

	INFO_MSG("Sorting NFW parameters"); 

	nBins = Settings.n_bins;
	nHaloesCut = n_haloes_per_criterion();

#ifdef WITH_MPI
		nHaloes = Settings.n_haloes; 
#else
		nHaloes = Settings.n_threshold; 
#endif

	gbin_x = (double*) calloc(nBins, sizeof(double));	

	nfw_gof_int_y = (int*) calloc(nBins-1, sizeof(int));	
	nfw_gof_bin_x = (double*) calloc(nBins-1, sizeof(double));	
	nfw_gof_double_y = (double*) calloc(nBins-1, sizeof(double));	
	nfw_gof_err_y = (double*) calloc(nBins-1, sizeof(double));	

	nfw_gof = (double*) calloc(nHaloesCut, sizeof(double));	

		for(i=0; i<nHaloes; i++)
		{
			if(halo_condition(i) == 1)
			{
				nfw_gof[m] = Haloes[i].fit_nfw.gof;
				m++;
			}
		}

			gMax = F_MAX*maximum(nfw_gof, nHaloesCut);  
			gMin = F_MIN*minimum(nfw_gof, nHaloesCut);
			delta_g = (gMax-gMin)/nBins; 
			gnorm = 1./(delta_g*nHaloesCut);

			gbin_x = lin_stepper(gMin, gMax, nBins);
			lin_bin(nfw_gof, gbin_x, nBins, nHaloesCut, nfw_gof_int_y);	


			ghalfstep=(gbin_x[1]-gbin_x[0])*0.5;

			for(i=0; i<nBins-1; i++)
			{	
				gvalue = (double) nfw_gof_int_y[i];
				nfw_gof_bin_x[i]=gbin_x[i]+ghalfstep;
				nfw_gof_double_y[i]=gnorm*gvalue; 
			}

		for(i=0; i<nBins-1; i++)
		{		
			HaloProperties[HALO_INDEX].p_fit_nfw.gof[i]=nfw_gof_bin_x[i];
			HaloProperties[HALO_INDEX].p_fit_nfw.p_gof[i]=nfw_gof_double_y[i];
		}	

	free(nfw_gof_err_y);
	free(nfw_gof_bin_x);
	free(nfw_gof);
	free(gbin_x);
}



void sort_mass_relations()
{
	int i=0, m=0, *n_mass, nBins, nHaloesCut, nHaloes; 
	double mMax, mMin;
	double *vel_bin, *vel_err, *mass, *vel, *mass_bin; 
	double *vir_bin, *vir_err, *vir; 
	double *conc_bin, *conc_err, *conc;
	double *avg_sub_bin, *avg_sub_err, *avg_sub;
	double *lambda_bin, *lambda_err, *lambda;
	double *shape_bin, *shape_err, *shape;
	double *triax_bin, *triax_err, *triax;
	double *chi_bin, *chi_err, *chi;
	double *gof_bin, *gof_err, *gof;
	double *per_bin, *per_err, *per;
	double c_0, v_0;
	double *param_conc, *param_vel;

	INFO_MSG("Sorting halo radial velocities and concentrations");

	nBins = Settings.n_bins;

	nHaloesCut = n_haloes_per_criterion();
	
#ifdef WITH_MPI
		nHaloes=Settings.n_haloes; 
#else
		nHaloes=Settings.n_threshold; 
#endif

		lambda = (double*) calloc(nHaloesCut, sizeof(double));	
		triax = (double*) calloc(nHaloesCut, sizeof(double));	
		shape = (double*) calloc(nHaloesCut, sizeof(double));	
		mass = (double*) calloc(nHaloesCut, sizeof(double));	
		conc = (double*) calloc(nHaloesCut, sizeof(double));	
		avg_sub = (double*) calloc(nHaloesCut, sizeof(double));	
		vel = (double*) calloc(nHaloesCut, sizeof(double));	
		vir = (double*) calloc(nHaloesCut, sizeof(double));	
		gof = (double*) calloc(nHaloesCut, sizeof(double));	
		per = (double*) calloc(nHaloesCut, sizeof(double));	
		chi = (double*) calloc(nHaloesCut, sizeof(double));	

		mass_bin = (double*) calloc(nBins, sizeof(double));	
		n_mass = (int*) calloc(nBins-1, sizeof(int));	
		vel_bin = (double*) calloc(nBins-1, sizeof(double));	
		vel_err = (double*) calloc(nBins-1, sizeof(double));	
		vir_bin = (double*) calloc(nBins-1, sizeof(double));	
		vir_err = (double*) calloc(nBins-1, sizeof(double));	
		conc_bin = (double*) calloc(nBins-1, sizeof(double));	
		conc_err = (double*) calloc(nBins-1, sizeof(double));	
		avg_sub_bin = (double*) calloc(nBins-1, sizeof(double));	
		avg_sub_err = (double*) calloc(nBins-1, sizeof(double));	
		lambda_bin = (double*) calloc(nBins-1, sizeof(double));	
		lambda_err = (double*) calloc(nBins-1, sizeof(double));	
		shape_bin = (double*) calloc(nBins-1, sizeof(double));	
		shape_err = (double*) calloc(nBins-1, sizeof(double));	
		triax_bin = (double*) calloc(nBins-1, sizeof(double));	
		triax_err = (double*) calloc(nBins-1, sizeof(double));	
		chi_bin = (double*) calloc(nBins-1, sizeof(double));	
		chi_err = (double*) calloc(nBins-1, sizeof(double));	
		per_bin = (double*) calloc(nBins-1, sizeof(double));	
		per_err = (double*) calloc(nBins-1, sizeof(double));	
		gof_bin = (double*) calloc(nBins-1, sizeof(double));	
		gof_err = (double*) calloc(nBins-1, sizeof(double));	

		param_vel = (double*) calloc(2, sizeof(double)); 
		param_conc = (double*) calloc(2, sizeof(double)); 

		for(i=0; i<nHaloes; i++)
		{
			if(halo_condition(i) == 1 && m < nHaloesCut)
			{
			//	D_PRINT("M=", m);
			//	D_PRINT("i=", i);
				vel[m] = Haloes[i].Vmax;
				vir[m] = Haloes[i].abs_th_vir;
				chi[m] = Haloes[i].fit_nfw.chi;
				gof[m] = Haloes[i].fit_nfw.gof;
				per[m] = Haloes[i].fit_nfw.per;
				mass[m] = Haloes[i].Mvir; 
				//conc[m] = Haloes[i].c_nfw;
				//conc[m] = Haloes[i].Rvir/Haloes[i].fit_nfw.rs;
				//conc[m] = Haloes[i].Rvir/Haloes[i].r2;
				conc[m] = Haloes[i].fit_nfw.c;
				avg_sub[m] = (double) Haloes[i].n_satellites;
				lambda[m] = Haloes[i].lambda;
				triax[m] = Haloes[i].triax;
				shape[m] = Haloes[i].shape;

				//if(conc[m] == -1) 
				//	conc[m] = Haloes[i].c;
				m++;
			}
		}

			nHaloesCut=m;
		
			mMin = F_MIN * minimum(mass, nHaloesCut);
			mMax = F_MAX * maximum(mass, nHaloesCut);

			mass_bin = log_stepper(mMin, mMax, nBins);

			average_bin(mass, vel, mass_bin, vel_bin, vel_err, nBins, nHaloesCut);
			average_bin(mass, vir, mass_bin, vir_bin, vir_err, nBins, nHaloesCut);
			average_bin(mass, conc, mass_bin, conc_bin, conc_err, nBins, nHaloesCut);
			average_bin(mass, avg_sub, mass_bin, avg_sub_bin, avg_sub_err, nBins, nHaloesCut);
			average_bin(mass, triax, mass_bin, triax_bin, triax_err, nBins, nHaloesCut);
			average_bin(mass, shape, mass_bin, shape_bin, shape_err, nBins, nHaloesCut);
			average_bin(mass, lambda, mass_bin, lambda_bin, lambda_err, nBins, nHaloesCut);
			average_bin(mass, per, mass_bin, per_bin, per_err, nBins, nHaloesCut);
			average_bin(mass, gof, mass_bin, gof_bin, gof_err, nBins, nHaloesCut);
			average_bin(mass, chi, mass_bin, chi_bin, chi_err, nBins, nHaloesCut);
			lin_bin(mass, mass_bin, nBins, nHaloesCut, n_mass);	

			for(i=0; i<nBins-1; i++)	
			{
				HaloProperties[HALO_INDEX].mass[i]=0.5*(mass_bin[i]+mass_bin[i+1]);
				HaloProperties[HALO_INDEX].n_entry[i]=n_mass[i];
				HaloProperties[HALO_INDEX].vel[i]=vel_bin[i];
				HaloProperties[HALO_INDEX].halo.virial[i]=vir_bin[i];
				HaloProperties[HALO_INDEX].conc[i]=conc_bin[i];
				HaloProperties[HALO_INDEX].sub[i]=avg_sub_bin[i];
				HaloProperties[HALO_INDEX].halo.lambda[i]=lambda_bin[i];
				HaloProperties[HALO_INDEX].halo.shape[i]=shape_bin[i];
				HaloProperties[HALO_INDEX].halo.triax[i]=triax_bin[i];
				HaloProperties[HALO_INDEX].fit_nfw.chi[i]=chi_bin[i];
				HaloProperties[HALO_INDEX].fit_nfw.per[i]=per_bin[i];
				HaloProperties[HALO_INDEX].fit_nfw.gof[i]=gof_bin[i];
			}

			INFO_MSG("Fitting Mass-Concentration relation to a power law");
			param_conc[0] = 0.1; 
			param_conc[1] = pow(1.e-14, param_conc[0]);
			param_conc = best_fit_power_law(mass_bin, conc_bin, conc_err, nBins-1, param_conc);
			HaloProperties[HALO_INDEX].c_0 = param_conc[1] * pow(1.e14, param_conc[0]);
			HaloProperties[HALO_INDEX].c_beta = param_conc[0];
	
			INFO_MSG("Fitting Mass-Radial velocity relation to a power law");
			param_vel[0] = 0.5; 
			param_vel[1] = pow(1.e-14, param_vel[0]);
			param_vel = best_fit_power_law(mass_bin, vel_bin, vel_err, nBins-1, param_vel);
			HaloProperties[HALO_INDEX].vel_0 = param_vel[1] * pow(1.e14, param_vel[0]);
			HaloProperties[HALO_INDEX].vel_beta = param_vel[0];
			HaloProperties[HALO_INDEX].avgSub = average(avg_sub, nHaloesCut);

	free(mass_bin);
	free(mass);
	free(conc_err); 
	free(conc_bin); 
	free(conc); 
	free(avg_sub_err); 
	free(avg_sub_bin); 
	free(avg_sub); 
	free(triax_err); 
	free(triax_bin); 
	free(triax); 
	free(shape_err); 
	free(shape_bin); 
	free(shape); 
	free(lambda_err); 
	free(lambda_bin); 
	free(lambda); 
	free(vel_err); 
	free(vel_bin); 
	free(vel); 
	free(vir_err); 
	free(vir_bin); 
	free(vir); 
/*
*/
}


#ifdef GAS
void sort_T_mass_function()
{
	int nBins, nHaloesCut, nHaloes, i=0, m=0; 
	int *temp_bin_y, *cum_temp; 
	double *temp, *temp_bin_x; 
	double tMax, tMin, volume;

	INFO_MSG("Sorting X-ray temperature function"); 

	nBins = Settings.n_bins;
	nHaloesCut = n_haloes_per_criterion();

#ifdef WITH_MPI
		nHaloes = Settings.n_haloes; 
#else
		nHaloes = Settings.n_threshold; 
#endif

	temp = (double*) calloc(nHaloes, sizeof(double));	
	temp_bin_x = (double*) calloc(nBins, sizeof(double));	
	temp_bin_y = (int*) calloc(nBins-1, sizeof(int));	
	cum_temp = (int*) calloc(nBins-1, sizeof(int));	

		for(i=0; i<nHaloes; i++)
		{
			if(halo_condition(i) == 1)
			{
				temp[m] = Haloes[i].gas_only.T_mw;
				m++;
			}
		}

			tMax = F_MAX * maximum(temp, nHaloesCut);  
			tMin = F_MIN * minimum(temp, nHaloesCut);

			temp_bin_x = lin_stepper(tMin, tMax, nBins);
			lin_bin(temp, temp_bin_x, nBins, nHaloesCut, temp_bin_y);	
	
			cum_bin(temp_bin_y, cum_temp, nBins-1);
	
		volume = pow3(Settings.box_size);

		for(i=0; i<nBins-1; i++)
		{		
			HaloProperties[HALO_INDEX].T[i] = 0.5 * (temp_bin_x[i] + temp_bin_x[i+1]);
			HaloProperties[HALO_INDEX].n_T[i] = cum_temp[i] / volume;
		}

	free(temp_bin_y);
	free(temp_bin_x);
	free(temp);
}



void sort_gas_relations()
{
 	int nBins=0, nHaloes=0, nHaloesCut=0, i=0, m=0;
	double *params, *temperature, *temperature_bin, *temperature_err, *mass_bin, *mass;
	double *gas_fraction, *gas_fraction_bin, *gas_fraction_err;
	double *diff_cm, *diff_cm_bin, *diff_cm_err;
	double *lambda, *lambda_bin, *lambda_err;
	double *beta, *beta_bin, *beta_err;
	double *costh, *costh_bin, *costh_err;
	double *shape, *shape_bin, *shape_err;
	double *triax, *triax_bin, *triax_err;
	double *gas_ekin, *gas_ekin_bin, *gas_ekin_err;
	double *gas_virial, *gas_virial_bin, *gas_virial_err;
	double *dm_ekin, *dm_ekin_bin, *dm_ekin_err;
	double *dm_virial, *dm_virial_bin, *dm_virial_err;

	double M_0, mMax, mMin;
	
	INFO_MSG("Sorting mass temperature and baryon fraction");

	nBins = Settings.n_bins;
	nHaloesCut = n_haloes_per_criterion();

#ifdef WITH_MPI
		nHaloes = Settings.n_haloes; 
#else
		nHaloes = Settings.n_threshold; 
#endif
	
	Settings.tick=0;
	params = (double*) calloc(2,sizeof(double));

	mass = (double*) calloc(nHaloesCut, sizeof(double));	
	mass_bin = (double*) calloc(nBins, sizeof(double));	
	temperature = (double*) calloc(nHaloesCut, sizeof(double));	
	temperature_bin = (double*) calloc(nBins-1, sizeof(double));	
	temperature_err = (double*) calloc(nBins-1, sizeof(double));	
	gas_fraction = (double*) calloc(nHaloesCut, sizeof(double));	
	gas_fraction_bin = (double*) calloc(nBins-1, sizeof(double));	
	gas_fraction_err = (double*) calloc(nBins-1, sizeof(double));	
	diff_cm = (double*) calloc(nHaloesCut, sizeof(double));	
	diff_cm_bin = (double*) calloc(nBins-1, sizeof(double));	
	diff_cm_err = (double*) calloc(nBins-1, sizeof(double));	
	lambda = (double*) calloc(nHaloesCut, sizeof(double));	
	lambda_bin = (double*) calloc(nBins-1, sizeof(double));	
	lambda_err = (double*) calloc(nBins-1, sizeof(double));	
	beta = (double*) calloc(nHaloesCut, sizeof(double));	
	beta_bin = (double*) calloc(nBins-1, sizeof(double));	
	beta_err = (double*) calloc(nBins-1, sizeof(double));	
	triax = (double*) calloc(nHaloesCut, sizeof(double));	
	triax_bin = (double*) calloc(nBins-1, sizeof(double));	
	triax_err = (double*) calloc(nBins-1, sizeof(double));	
	shape = (double*) calloc(nHaloesCut, sizeof(double));	
	shape_bin = (double*) calloc(nBins-1, sizeof(double));	
	shape_err = (double*) calloc(nBins-1, sizeof(double));	
	costh = (double*) calloc(nHaloesCut, sizeof(double));	
	costh_bin = (double*) calloc(nBins-1, sizeof(double));	
	costh_err = (double*) calloc(nBins-1, sizeof(double));	
	gas_ekin = (double*) calloc(nHaloesCut, sizeof(double));	
	gas_ekin_bin = (double*) calloc(nBins-1, sizeof(double));	
	gas_ekin_err = (double*) calloc(nBins-1, sizeof(double));	
	gas_virial = (double*) calloc(nHaloesCut, sizeof(double));	
	gas_virial_bin = (double*) calloc(nBins-1, sizeof(double));	
	gas_virial_err = (double*) calloc(nBins-1, sizeof(double));	
	dm_ekin = (double*) calloc(nHaloesCut, sizeof(double));	
	dm_ekin_bin = (double*) calloc(nBins-1, sizeof(double));	
	dm_ekin_err = (double*) calloc(nBins-1, sizeof(double));	
	dm_virial = (double*) calloc(nHaloesCut, sizeof(double));	
	dm_virial_bin = (double*) calloc(nBins-1, sizeof(double));	
	dm_virial_err = (double*) calloc(nBins-1, sizeof(double));	

	INFO_MSG("Memory allocated");

		for(i=0; i<nHaloes; i++)
		{
			if(halo_condition(i) == 1 && Haloes[i].gas.N > 0)
			{
				mass[m] = Haloes[i].Mvir;
				gas_fraction[m] = Haloes[i].gas_only.b_fraction;
				lambda[m] = Haloes[i].gas_only.lambda;
				costh[m] = Haloes[i].gas_only.gas_dm_costh;
				beta[m] = Haloes[i].fit_beta.beta;
				temperature[m] = Haloes[i].gas_only.T_mw * Cosmo.h;
				gas_ekin[m] = Haloes[i].gas.Ekin;
				gas_virial[m] = Haloes[i].gas.vir;
				dm_ekin[m] = Haloes[i].dm.Ekin;
				dm_virial[m] = Haloes[i].dm.vir;
				shape[m] = Haloes[i].gas_only.shape;
				triax[m] = Haloes[i].gas_only.triax;
				diff_cm[m] = Haloes[i].gas_only.diff.cm;
				m++;
			}
		}
			nHaloesCut = m;

			mMin = F_MIN * minimum(mass, nHaloesCut);
			mMax = F_MAX * maximum(mass, nHaloesCut);

			//fprintf(stderr, "min=%e, max=%e\n", mMin, mMax);
			//fprintf(stderr, "n=%d, m=%d\n", n, m);
			mass_bin = log_stepper(mMin, mMax, nBins);

			average_bin(mass, temperature, mass_bin, temperature_bin, temperature_err, nBins, nHaloesCut);
			average_bin(mass, gas_fraction, mass_bin, gas_fraction_bin, gas_fraction_err, nBins, nHaloesCut);
			average_bin(mass, lambda, mass_bin, lambda_bin, lambda_err, nBins, nHaloesCut);
			average_bin(mass, beta, mass_bin, beta_bin, beta_err, nBins, nHaloesCut);
			average_bin(mass, triax, mass_bin, triax_bin, triax_err, nBins, nHaloesCut);
			average_bin(mass, shape, mass_bin, shape_bin, shape_err, nBins, nHaloesCut);
			average_bin(mass, gas_ekin, mass_bin, gas_ekin_bin, gas_ekin_err, nBins, nHaloesCut);
			average_bin(mass, gas_virial, mass_bin, gas_virial_bin, gas_virial_err, nBins, nHaloesCut);
			average_bin(mass, dm_ekin, mass_bin, dm_ekin_bin, dm_ekin_err, nBins, nHaloesCut);
			average_bin(mass, dm_virial, mass_bin, dm_virial_bin, dm_virial_err, nBins, nHaloesCut);
			average_bin(mass, diff_cm, mass_bin, diff_cm_bin, diff_cm_err, nBins, nHaloesCut);
			average_bin(mass, costh, mass_bin, costh_bin, costh_err, nBins, nHaloesCut);

			for(i=0; i<nBins-1; i++)
			{
				HaloProperties[HALO_INDEX].gas_T[i]=temperature_bin[i];
				HaloProperties[HALO_INDEX].gas_fraction[i]=gas_fraction_bin[i];
				HaloProperties[HALO_INDEX].gas_dm_costh[i]=costh_bin[i];
				HaloProperties[HALO_INDEX].gas_diff_cm[i]=diff_cm[i];
				HaloProperties[HALO_INDEX].gas.ekin[i]=gas_ekin_bin[i];
				HaloProperties[HALO_INDEX].gas.virial[i]=gas_virial_bin[i];
				HaloProperties[HALO_INDEX].dm.ekin[i]=dm_ekin_bin[i];
				HaloProperties[HALO_INDEX].dm.virial[i]=dm_virial_bin[i];
				HaloProperties[HALO_INDEX].gas.lambda[i]=lambda_bin[i];
				HaloProperties[HALO_INDEX].gas.beta[i]=beta_bin[i];
				HaloProperties[HALO_INDEX].gas.triax[i]=triax_bin[i];
				HaloProperties[HALO_INDEX].gas.shape[i]=shape_bin[i];
			}
	
			params[0] = 1.5; 
			params[1] = pow(1.e-14, params[0]);
			
			INFO_MSG("Fitting Mass-Temperature relation to a power law");
				
			params = best_fit_power_law(mass_bin, temperature_bin, temperature_err, nBins-1, params);

		M_0 = -log(params[1])/params[0]/14/log(10);

		HaloProperties[HALO_INDEX].T0=M_0;
		HaloProperties[HALO_INDEX].alpha=params[0];
		HaloProperties[HALO_INDEX].beta0=average(beta,nHaloesCut);;

//		fprintf(stdout, "M-Tx    a:%lf   M_0 10e+%f SM\n",a[0],M_0/log(10));
		
		// TODO free all stuff
	free(params);
	free(mass);
	free(mass_bin);
	free(lambda); 
	free(lambda_bin); 
	free(lambda_err); 
	free(beta); 
	free(beta_bin); 
	free(beta_err); 
	free(gas_ekin); 
	free(gas_ekin_bin); 
	free(gas_ekin_err); 
	free(gas_virial); 
	free(gas_virial_bin); 
	free(gas_virial_err); 
	free(dm_ekin); 
	free(dm_ekin_bin); 
	free(dm_ekin_err); 
	free(dm_virial); 
	free(dm_virial_bin); 
	free(dm_virial_err); 
	free(costh); 
	free(costh_bin); 
	free(costh_err); 
	free(triax); 
	free(triax_bin); 
	free(triax_err); 
	free(shape); 
	free(shape_bin); 
	free(shape_err); 
	free(temperature); 
	free(temperature_bin); 
	free(temperature_err); 
}



void sort_alignment_and_displacement()
{
	int i=0, m=0, nBins, nHaloesCut, nHaloes;
	int *array_costh_bin_y, *array_diff_cm_bin_y;
	double *array_costh, *array_diff_cm, *array_costh_bin, *array_diff_cm_bin; 
	double half_t, half_s, sMax, sMin, tMax, tMin, p_s, p_t;

	INFO_MSG("Sorting gas displacement and alignment");

	nBins = Settings.n_bins;
	nHaloesCut = n_haloes_per_criterion();

#ifdef WITH_MPI
	nHaloes = Settings.n_haloes; 
#else
	nHaloes = Settings.n_threshold; 
#endif

	Settings.tick=0;
	array_costh = (double*) calloc(nHaloesCut, sizeof(double));	
	array_diff_cm = (double*) calloc(nHaloesCut, sizeof(double));	
	array_costh_bin = (double*) calloc(nBins, sizeof(double));	
	array_diff_cm_bin = (double*) calloc(nBins, sizeof(double));	
	array_costh_bin_y = (int*) calloc(nBins-1, sizeof(int));	
	array_diff_cm_bin_y = (int*) calloc(nBins-1, sizeof(int));	

		for(i=0; i<nHaloes; i++)
		{
			if(halo_condition(i) == 1)
			{
				array_costh[m] = Haloes[i].gas_only.gas_dm_costh;
				array_diff_cm[m] = Haloes[i].gas_only.diff.cm;
				m++;
			}
		}

			sMax = F_MAX*maximum(array_costh, nHaloesCut); 
			sMin = F_MIN*minimum(array_costh, nHaloesCut);
			tMax = F_MAX*maximum(array_diff_cm, nHaloesCut); 
			tMin = F_MIN*minimum(array_diff_cm, nHaloesCut);

			array_costh_bin = lin_stepper(sMin, sMax, nBins);
			lin_bin(array_costh, array_costh_bin, nBins, nHaloesCut, array_costh_bin_y);	

			array_diff_cm_bin = log_stepper(tMin, tMax, nBins);
			lin_bin(array_diff_cm, array_diff_cm_bin, nBins, nHaloesCut, array_diff_cm_bin_y);	

			half_s = 0.5*(array_costh_bin[1]-array_costh_bin[0]);
			half_t = 0.5*(array_diff_cm_bin[1]-array_diff_cm_bin[0]);

			for(i=0; i<nBins-1; i++)	
			{
				p_s = (double) array_costh_bin_y[i]/nHaloesCut;
				p_t = (double) array_diff_cm_bin_y[i]/nHaloesCut;
				HaloProperties[HALO_INDEX].gas_dm_cth[i]=array_costh_bin[i]+half_s;
				HaloProperties[HALO_INDEX].cm[i]=array_diff_cm_bin[i]+half_t;
				HaloProperties[HALO_INDEX].p_gas_dm_cth[i]=p_s;
				HaloProperties[HALO_INDEX].p_cm[i]=p_t;
			}

	free(array_costh); 
	free(array_diff_cm);
	free(array_costh_bin); 
	free(array_diff_cm_bin);
	free(array_costh_bin_y); 
	free(array_diff_cm_bin_y);
}
#endif // Gas



void compute_halo_properties()
{
		sort_numerical_mass_function();
		sort_mass_relations();
		sort_nfw_parameters();
		sort_lambda_and_concentration();
		sort_shape_and_triaxiality();
#ifdef GAS
		sort_T_mass_function();
		sort_gas_relations();
		sort_alignment_and_displacement();
#endif
}
