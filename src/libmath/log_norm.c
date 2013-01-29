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

#include "log_norm.h"
#include "statistics.h"


double lognorm(double l, double l_0, double sig)
{
	double p, norm, arg;

		norm = 1./(l*sig*sqrt(2*3.14)), 
		arg  = log(l/l_0)/sig;
		p = norm*exp(-arg*arg*0.5);

	return p;
}


	/* Derivative of the distribution wrt the l_0 parameter */
double dl_lognorm(double l, double l_0, double sig)
{
	return -1./(2*sig*sig)*(2/l_0 - 2*log(l)/l_0)*lognorm(l,l_0,sig);
}


	/* Derivative of the distribution wrt the sig parameter */
double ds_lognorm(double l, double l_0, double sig)
{	
	return (-1./sig + pow(log(l/l_0),2)*pow(sig,-3))*lognorm(l,l_0,sig);
}



int lognorm_f(const gsl_vector *x, void *data, gsl_vector *f)
{
	size_t n=0, i=0;
	double arg, Yi, fset, t, norm, l_0, sig, *vx, *vy, *err;

		n = ((struct data *)data)->n;
		vx = ((struct data *)data)->x;
		vy = ((struct data *)data)->y;
		err = ((struct data *)data)->err;

	l_0 = gsl_vector_get(x,0);
	sig = gsl_vector_get(x,1);

		for(i=1; i<n; i++)
	{
		t = vx[i];
 		norm = 1./(t*sig*sqrt(2*3.14));
		arg  = log(t/l_0)/sig;
		Yi = norm*exp(-arg*arg*0.5);
		fset =  (Yi - vy[i])/err[i];
	
		if(fset!=fset) fset=0;
	
			gsl_vector_set(f, i, fset);
	}

	return GSL_SUCCESS;
}



int d_lognorm_f(const gsl_vector *x, void *data, gsl_matrix *J)
{
	size_t n=0, i=0;
	double log_norm, e1, e2, s, t, l_0, sig, dl_log_norm, ds_log_norm;
	double *vx, *err;

	n = ((struct data*)data)->n;	
	vx = ((struct data *)data)->x;
	err = ((struct data*)data)->err;

	l_0 = gsl_vector_get(x,0);
	sig = gsl_vector_get(x,1);

		for(i = 1; i < n; i++)
	{
		s = err[i];
		t = vx[i]; 
		
		log_norm = lognorm(t, l_0, sig);
		dl_log_norm = dl_lognorm(t, l_0, sig);
		ds_log_norm = ds_lognorm(t, l_0, sig);

		e1 = dl_log_norm; 
		e2 = ds_log_norm; 

		gsl_matrix_set(J,i,0,e1/s);
		gsl_matrix_set(J,i,1,e2/s);
	}

	return GSL_SUCCESS;
}



int fd_lognorm_f(const gsl_vector *x, void *data, gsl_vector * f, gsl_matrix *J)
{
	lognorm_f(x,data,f);
	d_lognorm_f(x,data,J);

	return GSL_SUCCESS;
}



double* best_fit_lognorm (double *array_values, int nHaloes, int nBins, 
		 	double *array_bins, double *array_values_bin, double *error) 
{
	double *params, l_0, sig;
	struct parameters par;
	struct data dat;

	l_0 = mean(array_values, nHaloes);
	sig = log_sigma(array_values, nHaloes);

		fprintf(stdout, "\nbest_fit_lognorm(). First guess parameters: l_0=%lf sigma=%lf\n", l_0, sig);

		dat.n = (size_t) nBins;
		dat.x = array_bins;
		dat.y = array_values_bin;
		dat.err = error; 

		/* First guess initialization */
		par.n = 2;
		par.guess_p = (double *) calloc(par.n, sizeof(double));
		par.guess_p[0] = l_0;
		par.guess_p[1] = sig;
		par.fitted_p = gsl_vector_alloc(par.n);

		/* Initialize function */
		gsl_multifit_function_fdf f;
       		f.f = &lognorm_f;
     		f.df = &d_lognorm_f;
     		f.fdf = &fd_lognorm_f;
       		f.n = dat.n;
       		f.p = par.n;
       		f.params = &dat;

		/* Do the fit */
		par.fitted_p = least_square_nl_fit(dat, par, f);

		/* Set the correctly fitted parameters */
		l_0 = gsl_vector_get(par.fitted_p,0);
		sig = gsl_vector_get(par.fitted_p,1);

		fprintf(stdout, "The best fit parameters for this distribution are: \n");
		fprintf(stdout, "l_0: %lf sigma: %lf \n", l_0, sig);

		/* Return the best fit parameters*/
		params = (double*) calloc(2, sizeof(double));
		params[0]=l_0;
		params[1]=sig; 

	return params;
}
