#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "halo_properties.h"
#include "subhalo_general.h"
#include "../general_variables.h"
#include "../general_functions.h"
#include "../libmath/mathtools.h"
#include "../libmath/log_norm.h"
#include "../libmath/power_law.h"
#include "../libcosmo/cosmological_relations.h"


void sort_host_axis_alignment_and_spatial_anisotropy()
{
	int i=0, j=0, k=0, totSub=0, totSubNmin=0, nBins=0;
	int *costh_bin_y, *cosphi_bin_y;
	double R, cMax, cMin, halfstep, cpMax, cpMin, halfstep2, ct, anis;
	double *costh, *costh_bin, *cosphi, *cosphi_bin; 

	totSub = Settings.n_subhaloes;
	totSubNmin = Settings.n_subhaloes_nmin;
	nBins = Settings.n_bins; 

	fprintf(stdout, "\nSorting sub halo radial alignment for %d sub haloes\n", 
		totSubNmin);

	fprintf(stdout, "\nSorting sub halo spatial anisotropy for %d sub haloes\n", 
		totSub);

		cosphi = (double*) calloc(totSub, sizeof(double));
		cosphi_bin = (double*) calloc(nBins, sizeof(double));
		cosphi_bin_y = (int*) calloc(nBins-1, sizeof(int));

		costh = (double*) calloc(totSubNmin, sizeof(double));
		costh_bin = (double*) calloc(nBins, sizeof(double));
		costh_bin_y = (int*) calloc(nBins-1, sizeof(int));

		for(i=0; i< totSub; i++)
		{
			k = SubHaloes[i].host;

			R = sqrt (
				pow(SubHaloes[i].Xc - Haloes[k].Xc,2) +
				pow(SubHaloes[i].Yc - Haloes[k].Yc,2) +
				pow(SubHaloes[i].Zc - Haloes[k].Zc,2) );

			ct = (
				SubHaloes[i].Eax*(SubHaloes[i].Xc - Haloes[k].Xc)+
				SubHaloes[i].Eay*(SubHaloes[i].Yc - Haloes[k].Yc)+ 
				SubHaloes[i].Eaz*(SubHaloes[i].Zc - Haloes[k].Zc)
				) / R; 
	
			anis = (
				Haloes[k].Eax*(SubHaloes[i].Xc - Haloes[k].Xc)+
				Haloes[k].Eay*(SubHaloes[i].Yc - Haloes[k].Yc)+ 
				Haloes[k].Eaz*(SubHaloes[i].Zc - Haloes[k].Zc)
				) / R;
			
				if(SubHaloes[i].n_part > Settings.n_min-1)
				{
					costh[j] = sqrt(ct*ct);
					j++;
				}

					cosphi[i] = sqrt(anis*anis); // sqrt(ct*ct);
				}


			cMax = 1.01*maximum(costh, totSubNmin); 
			cMin = minumum(costh, totSubNmin); 
			cpMax = 1.01*maximum(cosphi, totSub); 
			cpMin = minimum(cosphi, totSub);

			costh_bin   = lin_stepper(cMin, cMax, nBins);
			lin_bin(costh, costh_bin, nBins, totSubNmin, costh_bin_y);

		        cosphi_bin   = lin_stepper(cpMin, cpMax, nBins);
			lin_bin(cosphi, cosphi_bin, nBins, totSub, cosphi_bin_y);

			halfstep = 0.5*(costh_bin[1]-costh_bin[0]);
			halfstep2 = 0.5*(cosphi_bin[1]-cosphi_bin[0]);
	
			SubHaloZ.costh0=average(costh, totSubNmin);
			SubHaloZ.cosphi0=average(cosphi, totSub);
	
		fprintf(stdout, "AvgCosphi: %lf\n", SubHaloZ.cosphi0); 

	for(j=0; j < nBins-1; j++)
	{
		SubHaloZ.costh[j] = costh_bin[j] + halfstep;
		SubHaloZ.costh_count[j] = (double) costh_bin_y[j] / totSub;
		SubHaloZ.cosphi[j] = cosphi_bin[j] + halfstep2;
		SubHaloZ.cosphi_count[j] = (double) cosphi_bin_y[j] / totSub;
	}

	free(cosphi);
	free(cosphi_bin);
	free(cosphi_bin_y);
	free(costh);
	free(costh_bin);
	free(costh_bin_y);

	fprintf(stdout, "\n");
}



void sort_sub_shape_and_triaxiality()
{
	int i=0, j=0, nBins=Settings.n_bins, nHaloes=Settings.n_subhaloes_nmin; 
	int *array_shape_bin_y, *array_triax_bin_y;
	double sMax, sMin, tMax, tMin, p_s, p_t, half_s, half_t;
	double *array_shape, *array_triax, *array_shape_bin, *array_triax_bin;

	fprintf(stdout, "\nSorting subhalo shape and triaxiality.\n"); 

	Settings.tick=0;
	array_shape = (double*) calloc(nHaloes, sizeof(double));	
	array_triax = (double*) calloc(nHaloes, sizeof(double));	
	array_shape_bin = (double*) calloc(nBins, sizeof(double));	
	array_triax_bin = (double*) calloc(nBins, sizeof(double));	
	array_shape_bin_y = (int*) calloc(nBins-1, sizeof(int));	
	array_triax_bin_y = (int*) calloc(nBins-1, sizeof(int));	

		for(i=0; i<nHaloes; i++)
		{
			if(SubHaloes[i].n_part > Settings.n_min -1)
			{
				array_shape[j] = SubHaloes[i].shape;
				array_triax[j] = SubHaloes[i].triax;
				j++;
				}
					}
	
				array_shape = shellsort(array_shape, nHaloes);
				array_triax = shellsort(array_triax, nHaloes);

				sMax = array_shape[nHaloes-1]; sMin = array_shape[0];
				tMax = array_triax[nHaloes-1]; tMin = array_triax[0];

				array_shape_bin = lin_stepper(sMin, sMax, nBins);
				lin_bin(array_shape, array_shape_bin, nBins, nHaloes, array_shape_bin_y);	

				array_triax_bin = lin_stepper(tMin, tMax, nBins);
				lin_bin(array_triax, array_triax_bin, nBins, nHaloes, array_triax_bin_y);	

			half_s = 0.5*(array_shape_bin[1] - array_shape_bin[0]);
			half_t = 0.5*(array_triax_bin[1] - array_triax_bin[0]);

			SubHaloZ.s0 = average(array_shape, nHaloes);
			SubHaloZ.t0 = average(array_triax, nHaloes);

		for(i=0; i<nBins-1; i++)	
		{
			p_s = (double) array_shape_bin_y[i]/nHaloes;
			p_t = (double) array_triax_bin_y[i]/nHaloes;
			SubHaloZ.shape[i]=array_shape_bin[i]+half_s;
			SubHaloZ.triax[i]=array_triax_bin[i]+half_t;
			SubHaloZ.p_shape[i]=p_s;
			SubHaloZ.p_triax[i]=p_t;
			SubHaloZ.n_shape[i]=array_shape_bin_y[i];
			SubHaloZ.n_triax[i]=array_triax_bin_y[i];
		}

	free(array_shape); 
	free(array_triax);
	free(array_shape_bin); 
	free(array_triax_bin);
	free(array_shape_bin_y); 
	free(array_triax_bin_y);

	fprintf(stdout, "\n");
}



void sort_sub_lambda()
{
	int nBins=Settings.n_bins, nHaloes, i=0, m=0, *lambda_int_y; 
	double *bin_x, *params, *lambda, *lambda_err_y, *lambda_double_y;
	double l_0, sig, halfstep, lMax, lMin, delta_l, norm, value;

	fprintf(stdout, "\nSorting subhalo spin parameter.\n"); 

	Settings.tick=0;
	nHaloes=Settings.n_subhaloes_nmin;
	lambda = (double*) calloc(nHaloes, sizeof(double));	
	params = (double*) calloc(2, sizeof(double));

	bin_x = (double*) calloc(nBins, sizeof(double));	
	lambda_int_y = (int*) calloc(nBins-1, sizeof(int));	
	lambda_double_y = (double*) calloc(nBins-1, sizeof(double));	
	lambda_err_y = (double*) calloc(nBins-1, sizeof(double));	

		for(i=0; i<nHaloes; i++)
		{
			if(SubHaloes[i].n_part > Settings.n_min-1)
			{
				lambda[m] = SubHaloes[i].lambda;
				m++;
				}
			}
		
					lambda = shellsort(lambda, nHaloes);

					lMax = lambda[nHaloes-1];  lMin = lambda[0];
					delta_l = (lMax-lMin)/nBins; norm = 1./(delta_l*nHaloes);

					bin_x = lin_stepper(lMin, lMax, nBins);

					lin_bin(lambda, bin_x, nBins, nHaloes, lambda_int_y);	

			for(i=0; i<nBins-1; i++)
			{	
				value = lambda_int_y[i];
				lambda_err_y[i]=sqrt(value*norm); 
				lambda_double_y[i]=norm*value; 
			}

				params = best_fit_lognorm(lambda, nHaloes, nBins, bin_x, lambda_double_y, lambda_err_y);

				l_0 = params[0]; 
				sig = params[1];
				SubHaloZ.l_0=l_0;
				SubHaloZ.l_sig=sig;

				halfstep=(bin_x[1]-bin_x[0])*0.5;

		for(i=0; i<nBins-1; i++)
		{		
			SubHaloZ.l[i]=bin_x[i]+halfstep;
			SubHaloZ.p_l[i]=lambda_double_y[i];
			SubHaloZ.err_p_l[i]=lambda_err_y[i];
		}	

	free(lambda);
	free(lambda_int_y);
	free(lambda_double_y);
	free(lambda_err_y);
	free(bin_x);
	free(params);

	fprintf(stdout, "\n");
}



double *generate_average_from_random_set(double* all_r)
{
	int n=0, j=0, i=0, k=0, host=0, TOT_ITER=10, subDim=0, totSub=0,*subset=NULL; 
	double r=0, *all_r_new=NULL; 

	subDim=Settings.min_subhaloes; 
	totSub=Settings.n_subhaloes;

	fprintf(stdout, "\nGenerating random subset from complete set of SubHaloes.\n");

	subset = (int*) calloc(subDim, sizeof(int)); 
	all_r = (double*) calloc(subDim, sizeof(double));
	all_r_new = (double*) calloc(subDim, sizeof(double));

		for(j=0; j<TOT_ITER; j++)
		{
			subset = generate_random_subset(totSub, subDim, subset);
			subset = int_shellsort(subset, subDim);
	
			for(i=0; i<totSub; i++)
			{
				host = SubHaloes[i].host;
				if(i==subset[k] && host != SubHaloes[i].id) 
				{
					r = sqrt(pow(Haloes[host].Xc - SubHaloes[i].Xc, 2) 
						+ pow(Haloes[host].Yc - SubHaloes[i].Yc, 2) 
						+ pow(Haloes[host].Zc - SubHaloes[i].Zc, 2));
					all_r[k] = r/Haloes[host].Rvir;
					k++;
				}
			}

			all_r = shellsort(all_r, subDim);

			for(n=0; n<subDim; n++)
			{
				all_r_new[n] += all_r[n];
			}
		}

		for(i=0; i<subDim; i++) 
			all_r_new[i] /= (double) TOT_ITER;

	return all_r_new;
	
	free(subset);
	free(all_r);
	free(all_r_new);

	fprintf(stdout, "\n");
}



void n_r_subhalo_subset()
{
	int i, cumul=0, subDim=0, nBins=0, *n_r=NULL, *cum_n_r=NULL; 
	double *all_r=NULL, *R=NULL, rMin, rMax;

	subDim=Settings.min_subhaloes;
	nBins=Settings.n_bins;

	fprintf(stdout, "\nSubhalo subset N(<R)\n");

	Settings.tick=0;
	R = (double*) calloc(nBins, sizeof(double));
	n_r = (int*) calloc(nBins-1, sizeof(int)); 
	cum_n_r = (int*) calloc(nBins-1, sizeof(int)); 

		all_r = generate_average_from_random_set(all_r);

		rMin = all_r[0]; 
		rMax = all_r[subDim-1];
		R = lin_stepper(rMin, rMax, nBins);
		lin_bin(all_r, R, nBins, subDim, n_r);

		for(i=0; i<nBins-1; i++) 
		{
			cumul += n_r[i];
			cum_n_r[i] = cumul;
		}

	for(i=0; i<nBins-1; i++) 
	{
		SubHaloZ.r_sub_subset[i]=R[i];
		SubHaloZ.n_r_sub_subset[i]=n_r[i];
		SubHaloZ.cum_n_r_sub_subset[i]=cum_n_r[i]; 
	}

	free(R); 
	free(n_r); 
	free(all_r); 
	free(cum_n_r);

	fprintf(stdout, "\n");
}



void n_r_subhalo()
{
	int i=0, host=0, cumul=0, totSub=0, nBins=0; 
	int *cum_n_r=NULL, *n_r=NULL; 
	double r, rMin, rMax;
	double *R=NULL, *all_r=NULL;

	totSub = Settings.n_subhaloes;
	nBins = Settings.n_bins;

	fprintf(stdout, "\nSubhalo N(<R)\n");
	Settings.tick=0;

	R = (double*) calloc(nBins, sizeof(double));
	all_r = (double*) calloc(totSub, sizeof(double));
	n_r = (int*) calloc(nBins-1, sizeof(int)); 
	cum_n_r = (int*) calloc(nBins-1, sizeof(int)); 

		for(i=0; i<totSub; i++)
		{
			host = SubHaloes[i].host;
			if(host != SubHaloes[i].id) 
			{
				r = sqrt(pow(Haloes[host].Xc - SubHaloes[i].Xc, 2) 
					+ pow(Haloes[host].Yc - SubHaloes[i].Yc, 2) 
					+ pow(Haloes[host].Zc - SubHaloes[i].Zc, 2));
				all_r[i] = r/Haloes[host].Rvir;
				}
			}

			all_r = shellsort(all_r, totSub);
			rMin = all_r[0]; rMax = all_r[totSub-1];
			R = lin_stepper(rMin, rMax, nBins);
			lin_bin(all_r, R, nBins, totSub, n_r);

		for(i=0; i<nBins-1; i++) 
		{
			cumul += n_r[i];
			cum_n_r[i] = cumul;
		}

	for(i=0; i<nBins-1; i++) 
	{
		SubHaloZ.r_sub[i]=R[i];
		SubHaloZ.n_r_sub[i]=n_r[i];
		SubHaloZ.cum_n_r_sub[i]=cum_n_r[i]; 
	}

	free(R); 	
	free(n_r); 
	free(all_r); 
	free(cum_n_r);

	fprintf(stdout, "\n");
}



void sort_velocity_distribution()
{
	int totSub=0, i=0, k=0, nBins=0;
	int *vel_y=NULL;	
	double vHost=0, vel_0=0, halfstep=0, velMax=0, velMin=0, vDiff=0; 
	double *vel=NULL, *vel_x=NULL;

	totSub = Settings.n_subhaloes;
	nBins = Settings.n_bins;
	
	fprintf(stdout, "\nSorting sub halo velocity distribution.\n");
	Settings.tick=0;
	
	vel = (double*) calloc(totSub, sizeof(double));
	vel_x = (double*) calloc(nBins, sizeof(double));
	vel_y = (int*) calloc(nBins-1, sizeof(int));

		for(i=0; i<totSub; i++) 
		{
			k = SubHaloes[i].host;
			if(k != SubHaloes[i].id) 
			{
				vHost = sqrt(pow(Haloes[k].VXc, 2) + pow(Haloes[k].VYc, 2) + pow(Haloes[k].VZc, 2));
				vDiff = sqrt( 
						pow(Haloes[k].VXc - SubHaloes[i].VXc, 2)+
						pow(Haloes[k].VYc - SubHaloes[i].VYc, 2)+
						pow(Haloes[k].VZc - SubHaloes[i].VZc, 2)
					);
				vel[i]=vDiff/vHost;
				}	
			}	

			vel = shellsort(vel, totSub);
			vel_0 = average(vel, totSub); 
			SubHaloZ.vel_0=vel_0;

			velMin = vel[0] ; 
			velMax = vel[totSub-1];
			vel_x = lin_stepper(velMin, velMax, nBins);
			lin_bin(vel, vel_x, nBins, totSub, vel_y);
			halfstep=(vel_x[1]-vel_x[0])*0.5;

		for(i=0; i<nBins-1; i++)
		{
			SubHaloZ.vel_sub[i]=vel_x[i];
			SubHaloZ.n_vel_sub[i]=vel_y[i];
			SubHaloZ.p_vel_sub[i]=(double) vel_y[i]/totSub;
		}
	
	free(vel);
	free(vel_x);
	free(vel_y);

	fprintf(stdout, "\n");
}



void sort_sub_mass_function()
{
	int totSub=Settings.n_subhaloes, nBins=Settings.n_bins, i=0, k=0, cumSat=0;
	int *cum_n_mass, *n_mass;
	double halfstep, mMax, mMin; 
	double *mass, *mass_bin; 

	fprintf(stdout, "\nSorting sub halo mass function");
	Settings.tick=0;
	
	mass = (double*) calloc(totSub, sizeof(double));
	mass_bin = (double*) calloc(nBins, sizeof(double));
	n_mass = (int*) calloc(nBins-1, sizeof(int));
	cum_n_mass = (int*) calloc(nBins-1, sizeof(int));
			
	for(i=0; i<totSub; i++) 
		mass[i] = SubHaloes[i].Mvir;

		mass = shellsort(mass, totSub);
		mMin = mass[0]; 
		mMax = mass[totSub-1];
		mass_bin = log_stepper(mMin, mMax, nBins);

		lin_bin(mass, mass_bin, nBins, totSub, n_mass);	

		SubHaloZ.avgMass=average(mass, totSub);

			for(i=0; i<nBins-1; i++)
			{
				k = nBins - i - 1;
				cumSat += n_mass[k];
				cum_n_mass[k] = cumSat;
			}


		for(i=0; i<nBins-1; i++)
		{
			halfstep=(mass_bin[i+1]-mass_bin[i])*0.5;
			SubHaloZ.mass_sub[i]=mass_bin[i]+halfstep;
			SubHaloZ.n_sub[i]=n_mass[i];
			SubHaloZ.cum_n_sub[i]=cum_n_mass[i];
		}

		free(mass); 
		free(mass_bin); 
		free(n_mass); 
		free(cum_n_mass);

	fprintf(stdout, "\n");
}


void sort_eccentricity()
{	// FIXME
	int totSub = Settings.n_subhaloes, nBins = Settings.n_bins, i=0; 
	int *cum_n_ecc=NULL, *n_ecc=NULL;
	double e=0, eMax=0, eMin=0; 
	double *ecc=NULL, *ecc_bin=NULL; 

	fprintf(stdout, "\nSorting eccentricity");
	Settings.tick=0;
	
	ecc = (double*) calloc(totSub, sizeof(double));
	ecc_bin = (double*) calloc(nBins, sizeof(double));
	n_ecc = (int*) calloc(nBins-1, sizeof(int));
	cum_n_ecc = (int*) calloc(nBins-1, sizeof(int));

		for(i=0; i<totSub; i++) 
		{
			ecc[i] = e;
		}	

				ecc = shellsort(ecc, totSub);

			eMin = ecc[0]; eMax = ecc[totSub-1];
			ecc_bin = log_stepper(eMin, eMax, nBins);
			lin_bin(ecc, ecc_bin, nBins, totSub, n_ecc);

	for(i=0; i<nBins-1; i++)
		fprintf(stdout, "m_bin:%e  n_bin:%d\n", ecc_bin[i], n_ecc[i]); 

	for(i=0; i<nBins-1; i++)
	{
		SubHaloZ.ecc[i]=ecc_bin[i];
		SubHaloZ.n_ecc[i]=n_ecc[i];
	}

	free(ecc);
	free(ecc_bin);
	free(n_ecc);
	free(cum_n_ecc);

	fprintf(stdout, "\n");
}


void compute_subhalo_properties()
{
	fprintf(stdout,"\nComputing subhalo properties.\n");

	initialize_subhalo_properties_structure(); // Allocates memory for subhalo properties distributions
	
	init_subhalo_struct(); // Allocates memory for the N subhalo structures

	load_subhalo_list(); // Selects within the halo catalogue the SubHaloes and copies their properties into subhalo struct
			
			/* Do the actual computation of the properties */
		sort_host_axis_alignment_and_spatial_anisotropy();
		sort_velocity_distribution();
		sort_sub_shape_and_triaxiality();
		sort_sub_lambda();
		sort_sub_mass_function();
		//sort_eccentricity(); // TODO fix the eccentricity calculation

	fprintf(stdout,"\n");
}

