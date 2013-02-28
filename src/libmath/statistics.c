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

#include "math.h"


/*
 * Declare functions
 */
struct data;
struct parameters;

void print_state(size_t, gsl_multifit_fdfsolver *);


/*
 * Initialize functions
 */ 
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



void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
	int i=0;
	int size = s->x->size;
	
		for(i=0; i<size; i++) 
			fprintf(stdout, "%zu. s->x[%d] %f \n", iter, i,  gsl_vector_get(s->x, i));
}



double chi_square(double *y_th, double *y_ex, double *err, int size, int skip)
{
	int i=0;
	double chi=0; 
		
		for(i=skip; i<size; i++) 
			chi+=pow2(y_th[i-skip]-y_ex[i])/pow2(err[i]);  

	return chi; 
}



double goodness_of_fit(double *y_th, double *y_ex, int size, int skip)
{
	int i=0;
	double gof=0; 
		
		for(i=skip; i<size; i++) 
			gof+=pow2(log(y_th[i-skip])-log(y_ex[i]));

	gof /= (double) (size-skip);
	
	return gof; 
}



double percentage_error(double *y_th, double *y_ex, int size, int skip)
{
	int i=0;
	double per=0; 
		
		for(i=skip; i<size; i++) 
			per+=sqrt(pow2(y_th[i-skip]-y_ex[i]))/y_ex[i];

	per /= (double) (size-skip);
	
	return per; 
}
