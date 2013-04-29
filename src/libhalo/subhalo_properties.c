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
void n_r_subhalo(void);
void sort_velocity_distribution(void);
void sort_host_axis_alignment_and_spatial_anisotropy(void);

// FIXME
void sort_eccentricity(void);
void n_r_subhalo_subset(void);
double* generate_average_from_random_set(double*);


/*
 * Initialize functions
 */ 
void sort_host_axis_alignment_and_spatial_anisotropy()
{
	int h=0, i=0, j=0, k=0, l=0, m=0, n=0, totHost=0, totSub=0, totSubNmin=0, nBins=0;
	int *costh_bin_y, *cosphi_bin_y;
	double sum, R=0, cMax, cMin, halfstep, cpMax, cpMin, halfstep2, ct=0, anis=0;
	double *costh, *costh_bin, *cosphi, *cosphi_bin; 

	totSub = SubStructure.N_sub;
	totHost = SubStructure.N_host;
	totSubNmin = Settings.n_sub_threshold;
	nBins = (int) (F_SUB * Settings.n_bins); 

	fprintf(stdout, "\nSorting sub halo radial alignment for %d sub haloes in %d bins\n", totSubNmin, nBins);

	fprintf(stdout, "\nSorting sub halo spatial anisotropy for %d sub haloes in %d bins\n", totSub, nBins);

		cosphi = (double*) calloc(totSub, sizeof(double));
		cosphi_bin = (double*) calloc(nBins, sizeof(double));
		cosphi_bin_y = (int*) calloc(nBins-1, sizeof(int));

		costh = (double*) calloc(totSubNmin, sizeof(double));
		costh_bin = (double*) calloc(nBins, sizeof(double));
		costh_bin_y = (int*) calloc(nBins-1, sizeof(int));

		for(i=0; i<totHost; i++)
		{
			k = SubStructure.host[i].index;

			for(j=0; j<SubStructure.host[i].n_sub; j++)			 
			{
				sum = 0; ct = 0; R = 0; anis = 0;
				l = SubStructure.host[i].sub_index[j];

				// Center of mass distance host-satellite 
				for(m=0; m<3; m++)
					sum += pow2(Haloes[l].X[m] - Haloes[k].X[m]);
				R = sqrt(sum);

				// SubHalo axis alignment to the halo center
				for(m=0; m<3; m++)
					ct += Haloes[l].Ea[m]*(Haloes[l].X[m] - Haloes[k].X[m]);
				ct = ct / R; 

				// Alignment of subhaloes wrt the host major axis
				for(m=0; m<3; m++)
					anis += Haloes[k].Ea[m]*(Haloes[l].X[m] - Haloes[k].X[m]);
				anis = anis / R;

					if(halo_condition(l) == 1)
					{
				//		fprintf(stderr, "%d) ct=%f, an=%f R=%f\n", n, ct, anis, R);
						costh[n] = sqrt(ct*ct);
						n++;
					}

				cosphi[h] = sqrt(anis*anis); 
				h++;
			}
		}

			cMax = F_MAX * maximum(costh, totSubNmin); 
			cMin = F_MIN * minimum(costh, totSubNmin); 
			cpMax = F_MAX * maximum(cosphi, totSub); 
			cpMin = F_MIN * minimum(cosphi, totSub);

			costh_bin = lin_stepper(cMin, cMax, nBins);
			lin_bin(costh, costh_bin, nBins, totSubNmin, costh_bin_y);

		        cosphi_bin = lin_stepper(cpMin, cpMax, nBins);
			lin_bin(cosphi, cosphi_bin, nBins, totSub, cosphi_bin_y);

			halfstep  = 0.5 * (costh_bin[1]-costh_bin[0]);
			halfstep2 = 0.5 * (cosphi_bin[1]-cosphi_bin[0]);
	
			HaloProperties[HALO_INDEX].costh0 = average(costh, totSubNmin);
			HaloProperties[HALO_INDEX].cosphi0 = average(cosphi, totSub);
	
		fprintf(stdout, "AvgCosPhi: %lf\n", HaloProperties[HALO_INDEX].cosphi0); 
		fprintf(stdout, "AvgCosTh : %lf\n", HaloProperties[HALO_INDEX].costh0); 

	for(j=0; j<nBins-1; j++)
	{
		HaloProperties[HALO_INDEX].costh[j] = costh_bin[j] + halfstep;
		HaloProperties[HALO_INDEX].costh_count[j] = (double) costh_bin_y[j] / totSubNmin;
		HaloProperties[HALO_INDEX].cosphi[j] = cosphi_bin[j] + halfstep2;
		HaloProperties[HALO_INDEX].cosphi_count[j] = (double) cosphi_bin_y[j] / totSub;
	}

	free(cosphi);
	free(cosphi_bin);
	free(cosphi_bin_y);
	free(costh);
	free(costh_bin);
	free(costh_bin_y);
}


// Number of sub haloes per r-bin
void n_r_subhalo()
{
	int h=0, i=0, j=0, k=0, m=0, host=0, cum=0, totSub=0, totHost=0, nBins=0; 
	double r, sum=0, sub_frac=0;
	double *R=NULL, *all_r=NULL, *sub_r=NULL, *err_n_r=NULL, *all_n_r=NULL, *n_bin=NULL;
	int *bin_n_r;

	totSub = SubStructure.N_sub;
	totHost = SubStructure.N_host;
	nBins = (int) (F_SUB * Settings.n_bins); 

	fprintf(stdout, "\nSubhalo N(<R)\n");
	Settings.tick=0;

	all_r = (double*) calloc(totSub, sizeof(double));
	sub_r = (double*) calloc(totSub, sizeof(double));
	R = (double*) calloc(nBins, sizeof(double));
	n_bin = (double*) calloc(nBins-1, sizeof(double));
	bin_n_r = (int*) calloc(nBins-1, sizeof(int));
	all_n_r = (double*) calloc(nBins-1, sizeof(double));
	err_n_r = (double*) calloc(nBins-1, sizeof(double));

		R = log_stepper(RMIN, RMAX, nBins);

		for(i=0; i<totHost; i++)
		{
			host = SubStructure.host[i].index;

		if(SubStructure.host[i].n_sub > SUB_MIN)
		{
			for(j=0; j<SubStructure.host[i].n_sub; j++)
			{
				sum = 0;
				sub_frac = 1.; ///(double)SubStructure.host[i].n_sub;

				k = SubStructure.host[i].sub_index[j];

				for(m=0; m<3; m++)
					sum += pow2(Haloes[host].X[m] - Haloes[k].X[m]);

				r = sqrt(sum);

			//	fprintf(stderr, "host=%d sub=%d r=%f frac=%f\n", host, k, r, sub_frac);

				all_r[j] = r/Haloes[host].Rvir;
				sub_r[j] = sub_frac;
			}
	
			//average_bin(all_r, sub_r, R, bin_n_r, err_n_r, nBins, j);
			lin_bin(all_r, R, nBins, j, bin_n_r);

			for(k=0; k<nBins-1; k++)
			{
				if(bin_n_r[k] != bin_n_r[k])
				{
					bin_n_r[k] = 0;
				}

				if(bin_n_r[k] > 0.)
				{
					n_bin[k]++;
					all_n_r[k] += bin_n_r[k]; 
			//fprintf(stderr, "%d) bin_n_r=%lf, all_n_r=%lf n=%f\n", k, bin_n_r[k], all_n_r[k], n_bin[k]);
				}
				///SubStructure.host[i].n_sub;
			}
			}
		}

	//	RMIN = minimum(all_r, totSub); 
	//	RMAX = maximum(all_r, totSub);

	//	R = lin_stepper(RMIN, RMAX, nBins);

	//	average_bin(all_r, sub_r, R, bin_n_r, err_n_r, nBins, totSub);
			
	//	cum_bin(bin_n_r, all_n_r, nBins);
	
	sum = 0;
	double sum1 = 0;

	for(j=0; j<nBins-1; j++) 
	{
		//i = nBins -j -1;
		i = nBins -j -2;
		HaloProperties[HALO_INDEX].r_sub[i] = 0.5 * (R[i] + R[i+1]);
		HaloProperties[HALO_INDEX].n_r_sub[i] = all_n_r[i]/n_bin[i];
		sum1 += n_bin[i];
		sum += all_n_r[i]/sum1;
		//fprintf(stderr,"%d) sum=%f, sum1=%f\n", i, sum, sum1);
		HaloProperties[HALO_INDEX].cum_n_r_sub[i] = sum; 
	}

	free(R); 	
	free(all_r); 
	free(sub_r); 
	free(bin_n_r); 
	free(err_n_r); 
	free(all_n_r);

	fprintf(stdout, "\n");
}


// Sort subhalo velocity distribution and radial velocity distribution
void sort_velocity_distribution()
{
	int totSub=0, totHost=0, h=0, i=0, j=0, k=0, m=0, n=0, nBins=0;
	int *vel_y=NULL;	
	double vHost=0, vel_0=0, velMax=0, velMin=0, vDiff=0, sum=0, r=0; 
	double *vel=NULL, *vel_x=NULL, *all_r=NULL, *vel_r=NULL, *bin_r=NULL, *vel_err=NULL;

	totSub = SubStructure.N_sub;
	totHost = SubStructure.N_host;
	nBins = (int) (F_SUB * Settings.n_bins); 
	
	fprintf(stdout, "\nSorting sub halo velocity distribution.\n");
	Settings.tick=0;
	
	vel = (double*) calloc(totSub, sizeof(double));
	vel_x = (double*) calloc(nBins, sizeof(double));
	vel_y = (int*) calloc(nBins-1, sizeof(int));
	all_r = (double*) calloc(totSub, sizeof(double));
	bin_r = (double*) calloc(nBins, sizeof(double));
	vel_r = (double*) calloc(nBins-1, sizeof(double));
	vel_err = (double*) calloc(nBins-1, sizeof(double));

		for(i=0; i<totHost; i++) 
		{
			k = SubStructure.host[i].index;

			for(j=0; j<SubStructure.host[i].n_sub; j++) 
			{
				n = SubStructure.host[i].sub_index[j];
				r = 0;

				// Sort velocity difference
				sum = 0;
				for(m=0; m<3; m++)
					sum += pow2(Haloes[k].V[m]);
				vHost = sqrt(sum);

				sum = 0;
				for(m=0; m<3; m++)
					sum += pow2(Haloes[k].V[m] - Haloes[n].V[m]);
				vDiff = sqrt(sum);

				vel[h] = vDiff/vHost;

				// Sort velocity difference as a function of radius
				for(m=0; m<3; m++)
					sum += pow2(Haloes[k].X[m] - Haloes[n].X[m]);

				r = sqrt(sum);

				all_r[h] = r/Haloes[k].Rvir;
				h++;
			}	
		}	

			velMin = F_MIN * minimum(vel,totSub); 
			velMax = F_MAX * maximum(vel,totSub);
			vel_x = log_stepper(velMin, velMax, nBins);
			lin_bin(vel, vel_x, nBins, totSub, vel_y);

			//RMIN = minimum(all_r,totSub); 
			//RMAX = maximum(all_r,totSub);
			bin_r = log_stepper(RMIN, RMAX, nBins);

			average_bin(all_r, vel, bin_r, vel_r, vel_err, nBins, totSub);

			HaloProperties[HALO_INDEX].vel_0 = vel_0;

		for(i=0; i<nBins-1; i++)
		{
			HaloProperties[HALO_INDEX].vel_sub_r[i] = vel_r[i];
			HaloProperties[HALO_INDEX].vel_sub[i] = 0.5 * (vel_x[i]+vel_x[i+1]);
			HaloProperties[HALO_INDEX].n_vel_sub[i] = vel_y[i];
			HaloProperties[HALO_INDEX].p_vel_sub[i] = (double) vel_y[i]/totSub;
		}
	
	free(vel);
	free(vel_x);
	free(vel_y);

	fprintf(stdout, "\n");
}



double* generate_average_from_random_set(double* all_r)
{
		// FIXME
	int n=0, j=0, i=0, k=0, m=0, TOT_ITER=10, subDim=0, totSub=0,*subset=NULL; 
	uint64_t host=0;
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
}



void n_r_subhalo_subset()
{
	// FIXME
	int i, cumul=0, subDim=0, nBins=0, *n_r=NULL, *cum_n_r=NULL; 
	double *all_r=NULL, *R=NULL;

	subDim=Settings.n_sub_min; 
	nBins=Settings.n_bins;

	INFO_MSG("Subhalo subset N(<R)");

	Settings.tick=0;
	R = (double*) calloc(nBins, sizeof(double));
	n_r = (int*) calloc(nBins-1, sizeof(int)); 
	cum_n_r = (int*) calloc(nBins-1, sizeof(int)); 

		all_r = generate_average_from_random_set(all_r);

		R = lin_stepper(RMIN, RMAX, nBins);
		lin_bin(all_r, R, nBins, subDim, n_r);

		for(i=0; i<nBins-1; i++) 
		{
			cumul += n_r[i];
			cum_n_r[i] = cumul;
		}

	for(i=0; i<nBins-1; i++) 
	{
		HaloProperties[HALO_INDEX].r_sub_subset[i]=R[i];
		HaloProperties[HALO_INDEX].n_r_sub_subset[i]=n_r[i];
		HaloProperties[HALO_INDEX].cum_n_r_sub_subset[i]=cum_n_r[i]; 
	}

	free(R); 
	free(n_r); 
	free(all_r); 
	free(cum_n_r);

	fprintf(stdout, "\n");
}



void sort_eccentricity()
{
	// FIXME
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
		HaloProperties[HALO_INDEX].ecc[i]=ecc_bin[i];
		HaloProperties[HALO_INDEX].n_ecc[i]=n_ecc[i];
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

		sort_host_axis_alignment_and_spatial_anisotropy();
		sort_velocity_distribution();
		n_r_subhalo();

		//n_r_subhalo_subset();
		//sort_eccentricity();

	fprintf(stdout,"\n");
}

