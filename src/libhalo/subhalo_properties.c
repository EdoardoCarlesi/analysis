#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../libmath/math.h"
#include "../libcosmo/cosmo.h"
#include "../general_def.h"

#include "halo.h"


/*
 * Declare functions
 */
double* generate_average_from_random_set(double*);

void n_r_subhalo(void);
void n_r_subhalo_subset(void);
void sort_eccentricity(void);
void sort_velocity_distribution(void);
void sort_host_axis_alignment_and_spatial_anisotropy(void);


/*
 * Initialize functions
 */ 
void sort_host_axis_alignment_and_spatial_anisotropy()
{
	int i=0, j=0, k=0, m=0, totSub=0, totSubNmin=0, nBins=0;
	int *costh_bin_y, *cosphi_bin_y;
	double sum, R=0, cMax, cMin, halfstep, cpMax, cpMin, halfstep2, ct=0, anis=0;
	double *costh, *costh_bin, *cosphi, *cosphi_bin; 

	totSub = SubStructure.N_sub;
	totSubNmin = Settings.n_sub_threshold;
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
			k = Haloes[i].host;
	
			sum = 0; ct = 0; R = 0; anis = 0;

			for(m=0; m<3; m++)
				sum += pow2(Haloes[i].X[m] - Haloes[k].X[m]);

			R = sqrt(sum);

			for(m=0; m<3; m++)
				ct += Haloes[i].Ea[m]*(Haloes[i].X[m] - Haloes[k].X[m]);
	
			ct = ct / R; 
	
			for(m=0; m<3; m++)
				anis += Haloes[k].Ea[m]*(Haloes[i].X[m] - Haloes[k].X[m]);
	
			anis = anis / R;
			
				if(Haloes[i].n_part > Settings.n_min-1)
				{
					costh[j] = sqrt(ct*ct);
					j++;
				}

					cosphi[i] = sqrt(anis*anis); // sqrt(ct*ct);
				}


			cMax = 1.01*maximum(costh, totSubNmin); 
			cMin = minimum(costh, totSubNmin); 
			cpMax = 1.01*maximum(cosphi, totSub); 
			cpMin = minimum(cosphi, totSub);

			costh_bin   = lin_stepper(cMin, cMax, nBins);
			lin_bin(costh, costh_bin, nBins, totSubNmin, costh_bin_y);

		        cosphi_bin   = lin_stepper(cpMin, cpMax, nBins);
			lin_bin(cosphi, cosphi_bin, nBins, totSub, cosphi_bin_y);

			halfstep = 0.5*(costh_bin[1]-costh_bin[0]);
			halfstep2 = 0.5*(cosphi_bin[1]-cosphi_bin[0]);
	
			SubHaloProperties[HALO_INDEX].costh0=average(costh, totSubNmin);
			SubHaloProperties[HALO_INDEX].cosphi0=average(cosphi, totSub);
	
		fprintf(stdout, "AvgCosphi: %lf\n", SubHaloProperties[HALO_INDEX].cosphi0); 

	for(j=0; j < nBins-1; j++)
	{
		SubHaloProperties[HALO_INDEX].costh[j] = costh_bin[j] + halfstep;
		SubHaloProperties[HALO_INDEX].costh_count[j] = (double) costh_bin_y[j] / totSub;
		SubHaloProperties[HALO_INDEX].cosphi[j] = cosphi_bin[j] + halfstep2;
		SubHaloProperties[HALO_INDEX].cosphi_count[j] = (double) cosphi_bin_y[j] / totSub;
	}

	free(cosphi);
	free(cosphi_bin);
	free(cosphi_bin_y);
	free(costh);
	free(costh_bin);
	free(costh_bin_y);

	fprintf(stdout, "\n");
}



double* generate_average_from_random_set(double* all_r)
{
	int n=0, j=0, i=0, k=0, m=0, host=0, TOT_ITER=10, subDim=0, totSub=0,*subset=NULL; 
	double r=0, sum=0, *all_r_new=NULL; 

	subDim = Settings.n_sub_min; 
	totSub = SubStructure.N_sub;

	fprintf(stdout, "\nGenerating random subset from complete set of Haloes.\n");

	subset = (int*) calloc(subDim, sizeof(int)); 
	all_r = (double*) calloc(subDim, sizeof(double));
	all_r_new = (double*) calloc(subDim, sizeof(double));

		for(j=0; j<TOT_ITER; j++)
		{
			subset = generate_random_subset(totSub, subDim, subset);
			subset = int_shellsort(subset, subDim);
	
			for(i=0; i<totSub; i++)
			{
				host = Haloes[i].host;
				if(i==subset[k] && host != Haloes[i].id) 
				{
					sum = 0;

					for(m=0; m<3; m++)
						sum += pow2(Haloes[host].X[m] - Haloes[i].X[m]); 

					r = sqrt(sum);
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

	subDim=Settings.n_sub_min; 
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
		SubHaloProperties[HALO_INDEX].r_sub_subset[i]=R[i];
		SubHaloProperties[HALO_INDEX].n_r_sub_subset[i]=n_r[i];
		SubHaloProperties[HALO_INDEX].cum_n_r_sub_subset[i]=cum_n_r[i]; 
	}

	free(R); 
	free(n_r); 
	free(all_r); 
	free(cum_n_r);

	fprintf(stdout, "\n");
}



void n_r_subhalo()
{
	int i=0, m=0, host=0, cumul=0, totSub=0, nBins=0; 
	int *cum_n_r=NULL, *n_r=NULL; 
	double r, rMin, rMax, sum=0;
	double *R=NULL, *all_r=NULL;

	totSub = SubStructure.N_sub;
	nBins = Settings.n_bins;

	fprintf(stdout, "\nSubhalo N(<R)\n");
	Settings.tick=0;

	R = (double*) calloc(nBins, sizeof(double));
	all_r = (double*) calloc(totSub, sizeof(double));
	n_r = (int*) calloc(nBins-1, sizeof(int)); 
	cum_n_r = (int*) calloc(nBins-1, sizeof(int)); 

		for(i=0; i<totSub; i++)
		{
			host = Haloes[i].host;
			if(host != Haloes[i].id) 
			{
				sum = 0;

				for(m=0; m<3; m++)
					sum += pow2(Haloes[host].X[m] - Haloes[i].X[m]);

				r = sqrt(sum);

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
		SubHaloProperties[HALO_INDEX].r_sub[i]=R[i];
		SubHaloProperties[HALO_INDEX].n_r_sub[i]=n_r[i];
		SubHaloProperties[HALO_INDEX].cum_n_r_sub[i]=cum_n_r[i]; 
	}

	free(R); 	
	free(n_r); 
	free(all_r); 
	free(cum_n_r);

	fprintf(stdout, "\n");
}



void sort_velocity_distribution()
{
	int totSub=0, i=0, k=0, m=0, nBins=0;
	int *vel_y=NULL;	
	double vHost=0, vel_0=0, halfstep=0, velMax=0, velMin=0, vDiff=0, sum=0; 
	double *vel=NULL, *vel_x=NULL;

	totSub = SubStructure.N_sub;
	nBins = Settings.n_bins;
	
	fprintf(stdout, "\nSorting sub halo velocity distribution.\n");
	Settings.tick=0;
	
	vel = (double*) calloc(totSub, sizeof(double));
	vel_x = (double*) calloc(nBins, sizeof(double));
	vel_y = (int*) calloc(nBins-1, sizeof(int));

		for(i=0; i<totSub; i++) 
		{
			k = Haloes[i].host;
			if(k != Haloes[i].id) 
			{
				sum = 0;
				for(m=0; m<3; m++)
					sum += pow2(Haloes[k].V[m]);
				vHost = sqrt(sum);

				sum = 0;
				for(m=0; m<3; m++)
					sum += pow2(Haloes[k].V[m]- Haloes[i].V[m]);
				vDiff = sqrt( sum);

				vel[i]=vDiff/vHost;
				}	
			}	

			vel = shellsort(vel, totSub);
			vel_0 = average(vel, totSub); 
			SubHaloProperties[HALO_INDEX].vel_0=vel_0;

			velMin = vel[0] ; 
			velMax = vel[totSub-1];
			vel_x = lin_stepper(velMin, velMax, nBins);
			lin_bin(vel, vel_x, nBins, totSub, vel_y);
			halfstep=(vel_x[1]-vel_x[0])*0.5;

		for(i=0; i<nBins-1; i++)
		{
			SubHaloProperties[HALO_INDEX].vel_sub[i]=vel_x[i];
			SubHaloProperties[HALO_INDEX].n_vel_sub[i]=vel_y[i];
			SubHaloProperties[HALO_INDEX].p_vel_sub[i]=(double) vel_y[i]/totSub;
		}
	
	free(vel);
	free(vel_x);
	free(vel_y);

	fprintf(stdout, "\n");
}



void sort_eccentricity()
{	// FIXME
	int totSub=0, nBins=0, i=0; 
	int *cum_n_ecc=NULL, *n_ecc=NULL;
	double e=0, eMax=0, eMin=0; 
	double *ecc=NULL, *ecc_bin=NULL; 
	
	totSub = SubStructure.N_sub; 
	nBins = Settings.n_bins; 

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
		SubHaloProperties[HALO_INDEX].ecc[i]=ecc_bin[i];
		SubHaloProperties[HALO_INDEX].n_ecc[i]=n_ecc[i];
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


		Settings.use_sub = 1;
		compute_halo_properties();
		Settings.use_sub = 0;
	
/*		
		sort_host_axis_alignment_and_spatial_anisotropy();
		sort_velocity_distribution();
		n_r_subhalo();
		n_r_subhalo_subset();
		sort_eccentricity();
*/

	fprintf(stdout,"\n");
}

