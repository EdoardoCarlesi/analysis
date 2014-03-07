#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_vector.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../general_def.h"

#include "math.h"



double get_interpolated_value(double x_array[], double y_array[], int npts, double A)
{
	double V; 

		gsl_interp_accel *gia = gsl_interp_accel_alloc();	
#ifdef INTERP_CSPLINE
		gsl_spline *gs = gsl_spline_alloc(gsl_interp_cspline, npts);
#else
		gsl_spline *gs = gsl_spline_alloc(gsl_interp_akima, npts);
#endif
		gsl_spline_init(gs, x_array, y_array, npts);
	
		V = gsl_spline_eval(gs,A,gia);
		//V = gsl_interp_eval(gsl_interp_akima, x_array, y_array, A, gia);
		gsl_spline_free (gs);
		gsl_interp_accel_free (gia);

	return V;	
}



double* lin_stepper(double min, double max, int bins)
{
#ifdef PRINT_INFO
	fprintf(stdout, "\nReturning %d linear steps between %e and %e...\n", bins, min, max);
#endif
	int i=0;
	double s, *steps;

		steps = (double*) calloc(bins,sizeof(double));
		s = (max-min)/(bins-1);
		steps[0] = min;

		for(i=1; i<bins; i++) 
			steps[i]=steps[i-1]+s;

	return steps;
}



double* log_stepper(double a_0, double a_1, int num)
{
#ifdef PRINT_INFO
	fprintf(stdout, "\nReturning %d logarithmic steps between %e and %e...\n", num, a_0,a_1);
#endif
	int i=0;
	double a=a_0, step, *steps;

	steps = (double *) calloc(num, sizeof(double));
	step = -log(a_0/a_1)/((double) num-1); 

		for(i=0; i<num; i++)
		{
        		a=a_0*exp(step*i);
			steps[i]=a;
		}

	return steps;
}



double remove_outliers(double* vec,int size,int index)
{
	//TODO
	return 1;
}



double maximum(double *array, int size)
{
	int i=0;
	int index=0;
	double max=array[0];

		for(i=0;i<size;i++) 
			if(array[i]>max)
			{ 
				max = array[i];
				index=i;
			}

//	fprintf(stderr, "size=%d, max=%e\n", size, max);

	return max;
}



void maxima(double *array, int size, double *val, int *pos, int NMax)
{
	INFO_MSG("Searching for func maxima...");
	size_t i=0, j=0, index=0;
	
	gsl_vector *v = gsl_vector_calloc(size);
	gsl_permutation *p = gsl_permutation_calloc(size);
	
		for(i=0; i<size; i++)
			gsl_vector_set(v, i, array[i]);

				gsl_sort_vector_index(p, v);

	for(j=0; j<NMax; j++)
	{
		index = gsl_permutation_get(p, size-j-1);
		val[j] = array[index];
		pos[j] = index;
	}

}



void minima(double *array, int size, double *val, int *pos, int NMax)
{
	size_t i=0, j=0, index=0;
	
	gsl_vector *v = gsl_vector_calloc(size);
	gsl_permutation *p = gsl_permutation_calloc(size);
	
		for(i=0; i<size; i++)
			gsl_vector_set(v, i, array[i]);

				gsl_sort_vector_index(p, v);

	for(j=0; j<NMax; j++)
	{
		index = gsl_permutation_get(p, j);
		val[j] = array[index];
		pos[j] = index;
	}

}



double minimum(double *array, int size)
{
	int i=0;
	double min=array[0];
 
			for(i=0;i<size;i++) 
				if(array[i]<min) 
					min = array[i];
	return min;
}



double nonzero_minimum(double *array, int size)
{
	int i=0;
	double min=array[0];
 
			for(i=0;i<size;i++) 
				if(array[i]<min && array[i]!=0) min = array[i];
	return min;
}



double median(double *array, int size)
{
	int i=0, half_size; 
	double *avg, med; 

	half_size = (int) size / 2;
	avg = (double *) calloc(size, sizeof(double));
	avg = shellsort(array, size);
	
	med = array[half_size];
	
	return med;
}
	



double average(double *array, int size)
{
	int i=0, true_size; 
	double avg=0; 
	true_size = size;

		for (i=0; i<size; i++) 
		{
				if(array[i] != array[i])
				{
					true_size--;
				}
				else
				{
					avg += array[i];
				}
		}
	
	return avg/(double) true_size;
}	



int check_array(int input, int *arr, int dim)
{
	int i=0, found=0, obj=0;

		for(i=0; i<dim; i++)
		{
			obj=arr[i];
		if (obj==input) 
			{
				found=1;
			}
		}

	return found;
}



int *generate_random_subset(int Nmax, int Nmin, int *subset)
{
	int i=0, num=0;

	if(Nmin > Nmax) 
	{
		fprintf(stdout, "\nTrying to generate a subset %d larger than the dataset %d.\n",Nmin,Nmax);

		} else {

			srand ((int)time(NULL));

			while (i <= Nmin)
			{
				num = rand()%Nmax;

				if(check_array(num,subset,i)==0)
				{
					subset[i]=num;
					i++;
				}
			}	
		}

	return subset;
}



int int_maximum(int *array, int size)
{
		int max=array[0], i=0;
	
		for(i=0;i<size;i++)

			if(array[i]>max) 
				max = array[i];

	return max;
}


		/* Number the entries per bin in a given array */
void lin_bin(double* array, double* bins, int bin_size, int array_size, int* binned_array)
{
	int i=0, j=0;
#ifdef PRINT_INFO
	fprintf(stdout, "\nBinning into %d bins an array of size:%d...", bin_size, array_size);
#endif
	gsl_histogram *h = gsl_histogram_alloc (bin_size-1);
	gsl_histogram_set_ranges (h, bins, bin_size);
		
		for(i=0; i<array_size; i++)
			gsl_histogram_increment(h, array[i]);
		
		for(j=0; j<bin_size-1; j++)
			binned_array[j] = (int) h->bin[j];

	gsl_histogram_free(h);
}



void median_bin (double* array_x, double* array_y, 
	double* bins, double* binned_array, double *error_array, int bin_size, int array_size)
{
#ifdef PRINT_INFO
	fprintf(stdout, "\nBinning and taking the median of %d bins an array of size:%d...", bin_size, array_size);
#endif
	int i=0, j=0, half, *n_bins; 
	double **y_bins;
	
	y_bins = (double**) calloc(bin_size, sizeof(double*)) ;
	n_bins = (int*) calloc(bin_size, sizeof(int)) ;

#	pragma omp parallel for \
	private(j,i) shared(bin_size, array_size, y_bins, n_bins, array_x, array_y)
	for(j=0; j<bin_size-1; j++)
	{
		y_bins[j] = (double*) calloc(1, sizeof(double));
		n_bins[j] = 0;

		for(i=0; i<array_size; i++)
		{
			if(array_x[i] >= bins[j] && array_x[i] < bins[j+1])
			{
				n_bins[j]++;
				y_bins[j] = (double *) realloc(y_bins[j], (n_bins[j]+1)*sizeof(double));
				y_bins[j][n_bins[j]-1] = array_y[i];
			}
		}
	}

	for(j=0; j<bin_size-1; j++)
	{
		binned_array[j] = median(y_bins[j], n_bins[j]);
	}

	for(j=0; j<bin_size; j++)
		free(y_bins[j]);

	free(y_bins);
	free(n_bins);
}



void average_bin(double* array_x, double* array_y, double* bins, 
	double* binned_array, double *error_array, int bin_size, int array_size)
{
#ifdef PRINT_INFO
	fprintf(stdout, "\nBinning and averaging into %d bins an array of size:%d...", bin_size, array_size);
#endif
	int i=0, j=0; 

	gsl_histogram *h = gsl_histogram_alloc (bin_size-1);
	gsl_histogram *g = gsl_histogram_alloc (bin_size-1);
	gsl_histogram_set_ranges (h, bins, bin_size);
	gsl_histogram_set_ranges (g, bins, bin_size);
		
		for(i=0; i<array_size; i++)
		{
			gsl_histogram_increment  (h, array_x[i]);
			gsl_histogram_accumulate (g, array_x[i], array_y[i]);
		}

			for(j=0; j<bin_size-1; j++)
			{
				binned_array[j] = g->bin[j]/h->bin[j];
					if(h->bin[j]>0)	// Assuming poissonian error
						error_array[j] = g->bin[j] / sqrt(h->bin[j]);
			}

	gsl_histogram_free(h);
	gsl_histogram_free(g);
}


void double_cum_bin(double *n, double *nCum, int size)
{
	int i=0;
	
	nCum[size-1] = n[size-1];		

	for(i=1; i<size; i++)
		nCum[size-i-1] = n[size-i-1] + nCum[size-i];
/*
	for(i=0; i<size; i++)
	{
		fprintf(stdout, "Cumulative:%d, bin:%d\n", nCum[i], n[i]);
	}
*/
}


void cum_bin(int *n, int *nCum, int size)
{
	int i=0;
	
	nCum[size-1] = n[size-1];		

	for(i=1; i<size; i++)
		nCum[size-i-1] = n[size-i-1] + nCum[size-i];
/*
	for(i=0; i<size; i++)
	{
		fprintf(stdout, "Cumulative:%d, bin:%d\n", nCum[i], n[i]);
	}
*/
}



int* int_shellsort(int *array, int n)
{
	int i=0,j=0,inc=3, m=0, temp=0;

     while(inc>0)
     {
          for(i=0;i<n;i++)
	{
		m++;

                   j=i;
                   temp=array[i];

                   while((j>=inc)&&array[j-inc]>temp)
                   {
                        array[j]=array[j-inc];
                        j=j-inc;
                   }

                  array[j]=temp;
          }
          if(inc/2!=0)
           inc=inc/2;
          else if(inc==1)
           inc=0;
          else
           inc=1;
     }

	return array;
}



double* shellsort(double *array, int n)
{
	int i,j,inc=3, m=0;
	double temp;

     while(inc>0)
     {
          for(i=0;i<n;i++)
          {
		m++;

                   j=i;
                   temp=array[i];

                   while((j>=inc)&&array[j-inc]>temp)
                   {
                        array[j]=array[j-inc];
                        j=j-inc;
                   }

                   array[j]=temp;
          }
          if(inc/2!=0)
           inc=inc/2;
          else if(inc==1)
           inc=0;
          else
           inc=1;
     }

	return array;
}



double *invert_array(double *array, int size)
{
	int i=0;
	double *inv_array;

		inv_array = (double *) calloc(size, sizeof(double));

		for(i=0; i<size; i++)
		{
			inv_array[i]=array[size-i-1];
		}	

	return inv_array;
}



double solid_angle(double theta, void *p)
{
	return sin(theta);
}



double integrate_solid_angle(double the1,double the2, double phi1, double phi2)
{
	double result, error;

	the1 *= PI/180;
	the2 *= PI/180;
	phi1 *= PI/180;
	phi2 *= PI/180;

		gsl_function F;
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
		F.function=&solid_angle;
		F.params=0;

			gsl_integration_qags(&F, the1, the2, 0, 1e-4, 1000, w, &result, &error);

	return result*(phi2-phi1);
}
