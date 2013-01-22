#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "statistics.h"
#include "log_norm.h"


struct data;

struct parameters;

struct chi2_reduced Chi2;


double mean(double *vec, int size)
{
	int i=0; 
	double v0=0;

		for(i=0; i<size; i++) 
				v0 +=vec[i];

			v0=v0/size;

	return v0;
}



double sigma(double* vec, double v0, int size)
{
	int i=0;
	double sig=0; 

		for(i=0; i<size; i++)
			sig += (vec[i]-v0)*(vec[i]-v0);
			

return sqrt(sig/size);
}



double log_mean(double *vec, int size)
{
	int i=0; 
	double v0=0;

	for(i=0; i<size; i++) 
			v0 +=log(vec[i]);

		v0=v0/size;

	return v0;
}



double log_sigma(double* vec, int size)
{
		double sig=0; double mu; int i=0;

			mu=mean(vec, size);

		for(i=0; i<size; i++)
		{
			if(vec[i]!=0)
				sig += (log(vec[i])-log(mu))*(log(vec[i])-log(mu));
		}

	return sqrt(sig/size);
}



gsl_vector* least_square_nl_fit(struct data dat, struct parameters par, gsl_multifit_function_fdf f)
{
	size_t n, p;
	int status, i, iter = 0; 
	double *x_init;
	const gsl_rng_type * type;
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	
	n = dat.n;
	p = par.n;
		
	gsl_matrix *covar = gsl_matrix_alloc(p,p);

	x_init = (double *) calloc(par.n, sizeof(double));

		for(i=0; i<par.n; i++) 
			x_init[i] = par.guess_p[i];
		

		gsl_vector_view x = gsl_vector_view_array (x_init, p);
		gsl_rng * r;
		gsl_rng_env_setup();
     
    	   	type = gsl_rng_default;
    	   	r = gsl_rng_alloc (type);

  	        T = gsl_multifit_fdfsolver_lmsder;
    		s = gsl_multifit_fdfsolver_alloc(T, n, p);
    		gsl_multifit_fdfsolver_set(s, &f, &x.vector);

#ifdef PRINT_INFO
	       print_state (iter, s);
#endif

       do
         {
           iter++;
           status = gsl_multifit_fdfsolver_iterate (s);

#ifdef PRINT_INFO
	       print_state (iter, s);
#endif

           if (status)
             break;
     
           status = gsl_multifit_test_delta (s->dx, s->x,
                                             1e-6, 1e-6);
         }
       while (status == GSL_CONTINUE && iter < 500);
     
#ifdef PRINT_INFO
	       print_state (iter, s);
#endif

       gsl_multifit_covar (s->J, 0.0, covar);
     
       { 
       gsl_matrix_free (covar);
       gsl_rng_free (r);
     }

	return s->x;
       gsl_multifit_fdfsolver_free (s);
}


	/* This should generate the poissonian error */
double *numerical_sigma(int size, double *x, double *y)
{
	int i=0;
	double *n_s;

	n_s = (double*) calloc(size, sizeof(double));

	for(i=0; i<size; i++)
	{
		if(y[i]!=0){
			n_s[i] = y[i];
		} else {
			n_s[i]=n_s[i-1];
		}
	}

	return n_s;
}



void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
	int i=0;
	int size = s->x->size;
	
		for(i=0; i<size; i++) 
			fprintf(stdout, "%zu. s->x[%d] %f ", iter, i,  gsl_vector_get(s->x, i));

	fprintf(stderr, "\n");
}



double chi_square(double *y_th, double *y_ex, double *err, int size, int skip)
{
		int i=0;
		double chi=0; 
		
			for(i=0; i<size; i++) 
			{
				if(i>skip-1){
					chi+=pow(y_th[i]-y_ex[i],2.0)*pow(err[i],-2.0);  
				}			
			}

	return chi; 
}



void calculate_and_sort_chi2()
{
// TODO Make it more general
/*
	int massCut = Settings.n_threshold;
	int nBins = NFW.bins;
	double vir = Cosmo.virial; 
	double *a_x; double *a_y; double v;
	double cc, mm, rr, rho0, rs, r_back;
	int n=0; int i=0; int k=0;
	int *vir_haloes; int n_vir=0;

	for(k=0; k<massCut; k++){
	v=haloes[k].th_vir;
	if(sqrt(v*v)<vir) n_vir++;
	}

fprintf(stderr, "Total of %d virialized haloes over %d. \n", n_vir, massCut);
	vir_haloes=(int*) calloc(n_vir, sizeof(int));

	for(k=0; k<massCut; k++){
	v=haloes[k].th_vir;
	if(sqrt(v*v)<vir) {
	vir_haloes[i]=k; i++;
	}
	}
	
//for(k=0; k<haloes[0].n_bins; k++)fprintf(stderr, "r:%lf, rho:%lf \n", haloes[0].radius[k], haloes[0].rho[k]);

Chi2.chi2s = (double*) calloc(n_vir, sizeof(double));
fit_and_store_nfw_parameters_from_list(n_vir, vir_haloes, massCut);

i=0;
for(k=0; k<n_vir; k++){
i=vir_haloes[k];
Chi2.chi2s[k] = haloes[i].chi_nfw;
//fprintf(stderr, "%d) chi2:%lf \n", k, Chi2.chi2s[k]);
}
i=0; 
k=0;

Chi2.chi2s = shellsort(Chi2.chi2s, n_vir);
Chi2.bins = nBins;
Chi2.binned_chi = (double*) calloc(Chi2.bins, sizeof(double));
Chi2.outcomes = (int*) calloc(Chi2.bins, sizeof(int));
Chi2.binned_chi = lin_stepper(Chi2.chi2s[0], Chi2.chi2s[massCut-1], nBins);

FILE *nfwfit = fopen(Urls_internal.output_prefix, "w");
//FILE *nfwfit = fopen("nfw_chi_2_distribution-vdez0_1.dat", "w");
//for(n=0; n<massCut; n++) fprintf(stderr, "%d, chi: %lf \n", n, Chi2.chi2s[n]);
fprintf(stderr, "\nBinning chi2...");
i=0; n=0; k=0;

Chi2.outcomes = lin_bin(Chi2.chi2s, Chi2.binned_chi, Chi2.bins, n_vir, Chi2.outcomes);
fprintf(stderr, "Binned.\n");

for(n=0; n<Chi2.bins; n++) {
fprintf(nfwfit, "%lf %d %lf\n", Chi2.binned_chi[n], Chi2.outcomes[Chi2.bins-n-1], (double) Chi2.outcomes[Chi2.bins-n-1]/(double)massCut);
//fprintf(stderr, "%lf %d %lf\n", 0, 0, 0); //Chi2.binned_chi[n], 0, 0); //Chi2.outcomes[Chi2.bins-n-1], (double) Chi2.outcomes[Chi2.bins-n-1]/(double)massCut);
}
fprintf(stderr, "Done.\n");
fclose(nfwfit);
*/
}
