#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "halo.h"
#include "../libio/io.h"
#include "../libmath/math.h"
#include "../general_def.h"

#include <omp.h>


int Ncross;


void read_cross_halo(char *cross_url)
{
	int i=0, N=0;
	double a;
	char dummy[512];
	FILE *cross_file = fopen(cross_url, "r");

		crossHaloes = (struct halo*) calloc(1, sizeof(struct halo));

		N = get_lines(cross_file, cross_url);
		fprintf(stderr, "Reading cross correlation file %s, %d lines\n", cross_url, N);
		N = 0;
		
		while(!feof(cross_file))
		{
			N++;
			
			if(N<2)
			fgets(dummy, 512, cross_file);

			if(N>1)
			{
				crossHaloes = (struct halo*) realloc(crossHaloes, N * sizeof(struct halo));
				fscanf(cross_file, "%lf\t", &crossHaloes[N-1].Mvir);
				fscanf(cross_file, "%f\t", &crossHaloes[N-1].Rvir);
				fscanf(cross_file, "%f\t", &crossHaloes[N-1].Msub);
				fscanf(cross_file, "%d\t", &crossHaloes[N-1].n_part);
				fscanf(cross_file, "%f\t", &crossHaloes[N-1].X[0]);
				fscanf(cross_file, "%f\t", &crossHaloes[N-1].X[1]);
				fscanf(cross_file, "%f\t", &crossHaloes[N-1].X[2]);
				fscanf(cross_file, "%f\t", &crossHaloes[N-1].c_nfw);
				fscanf(cross_file, "%f\t", &crossHaloes[N-1].abs_th_vir);
				fscanf(cross_file, "%f\t", &crossHaloes[N-1].lambda);
				fscanf(cross_file, "%f\t", &crossHaloes[N-1].shape);
				fscanf(cross_file, "%f\t", &crossHaloes[N-1].triax);
#ifdef GAS
				fscanf(cross_file, "%f\t", &crossHaloes[N-1].gas_only.T_mw);
				fscanf(cross_file, "%f\t", &crossHaloes[N-1].gas_only.b_fraction);
#endif
			}
		}

	Ncross = N;

	fprintf(stderr, "Found %d cross haloes\n", N);
}



int find_cross_correlated_halo(int index)
{
	int id=-1, j=0, k=0;
	double x[3], Rvir, sum=0;
		
		for(j=0; j<3; j++)
			x[j] = Haloes[index].X[j];
			Rvir = Haloes[index].Rvir;

		for(k=0; k<Ncross; k++)
		{
			sum = 0;

			for(j=0; j<3; j++)
				sum += pow2(x[j] - crossHaloes[k].X[j]);
			
					sum = sqrt(sum);

				if(sum < 2*Rvir) 
			id = k;

		//	fprintf(stderr, "x=%f, X=%f, sum=%f, id=%d\n", x[0], crossHaloes[k].X[0], sum, id);
		}
					
	return id;
}



int find_ratios()
{
	int i=0, j=0, index=-1, check=0, n=0, nHaloes, count=0, colTot=10, nBins=4;
	double thr, ratio[10], *all_ratios[10], avg[10], sig[10], *avg_bin[10], *mass_bin;
	char out_url[200];

 	sprintf(out_url, "%s%s", Urls.output_prefix, "cross_correlation.dat");
	FILE *out_file = fopen(out_url, "w");

		fprintf(out_file, 
	"#Mlcdm(1)\tM(2)\t\tMsub(3)\t\tR(4)\t\tc(5)\t\tlambda(6)\tshape(7)\ttriax(8)\tvirial(9)\tTemp(10)\tfrac(11)\n");

	nHaloes = Settings.n_haloes; 
	thr = 2.25;

	mass_bin = (double *) calloc(nBins, sizeof(double));

		for(j=0; j<colTot; j++)
		{
			avg_bin[j] = (double *) calloc(nBins, sizeof(double));
			all_ratios[j] = (double *) calloc(1, sizeof(double));
		}

		for(i=0; i<nHaloes; i++)
		{
			if(halo_condition(i) == 1)
			{
				index = find_cross_correlated_halo(i);
				count = 0;

				if(index >= 0)		        
				{

					ratio[count++] = Haloes[i].Mvir / crossHaloes[index].Mvir;
					ratio[count++] = 0.0; //Haloes[i].Msub / crossHaloes[index].Msub;
					ratio[count++] = Haloes[i].Rvir / crossHaloes[index].Rvir;
#ifndef NO_PROFILES
					ratio[count++] = Haloes[i].fit_nfw.c / crossHaloes[index].c_nfw;
#endif
					ratio[count++] = Haloes[i].lambda / crossHaloes[index].lambda;
					ratio[count++] = Haloes[i].shape / crossHaloes[index].shape;
					ratio[count++] = Haloes[i].triax / crossHaloes[index].triax;
#ifdef GAS
					ratio[count++] = (-2*Haloes[i].dm.Ekin/Haloes[i].dm.Epot)
								/ crossHaloes[index].abs_th_vir;
					ratio[count++] = Haloes[i].gas_only.T_mw 
								/ crossHaloes[index].gas_only.T_mw;
					ratio[count++] = Haloes[i].gas_only.b_fraction 
								/ crossHaloes[index].gas_only.b_fraction;
#endif					
					check = 0;

					for(j=0; j<colTot; j++)
					{
						if(ratio[j] > thr)
							check = 1;
					}

#ifdef DOCHECK
					if(check == 0)
					{
#endif
					n++;
					cross_correlated_index = realloc(cross_correlated_index, n * sizeof(int));
					cross_correlated_index[n-1] = i;

					for(j=0; j<colTot; j++)
						all_ratios[j] = realloc(all_ratios[j], n * sizeof(double));

	
					fprintf(out_file, "%e", crossHaloes[index].Mvir);

						for(j=0; j<colTot; j++)
						{
							all_ratios[j][n-1] = ratio[j];

							fprintf(out_file, "\t%lf", ratio[j]);
						}
					
					fprintf(out_file, "\n");
#ifdef DOCHECK
					}
#endif
				}
			}
		}
	
		double mMax = 1.2e15;
		double mMin = 1.e7;

		mass_bin = log_stepper(mMin,mMax,nBins);

		fprintf(out_file, "#Avg  =\t");
				for(j=0; j<colTot; j++)
				{	
					avg[j] = average(all_ratios[j], n);
//					avb_bin = 
					fprintf(out_file, "\t%lf", avg[j]);
				}
					fprintf(out_file, "\n");
			
		fprintf(out_file, "#sigma=\t");
				for(j=0; j<colTot; j++)
				{	
					sig[j] = sigma(all_ratios[j], avg[j], n);
					fprintf(out_file, "\t%lf", pow2(sig[j]));
				}
					fprintf(out_file, "\n");

	fprintf(stderr, "Correlated %d haloes\n", n);
		
	fclose(out_file);
}



void cross_correlation(char *url)
{
	cross_correlated_index = (int *) calloc(1, sizeof(int));

	read_cross_halo(url);

	find_ratios();		
}
