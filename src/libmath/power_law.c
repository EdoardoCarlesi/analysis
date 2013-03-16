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

#include "math.h"

/*
 * Declare functions
 */
int power_law_f(const gsl_vector*, void *, gsl_vector*);
int d_power_law_f(const gsl_vector*, void *, gsl_matrix*);
int fd_power_law_f(const gsl_vector*, void *, gsl_vector*, gsl_matrix*);


/*
 * Initialize functions
 */ 
int power_law_f(const gsl_vector *x, void *data, gsl_vector *f)
{
	size_t n, i; 
	double t, a, b, Yi, fset, *vx, *vy, *err;

	n = ((struct data *)data)->n;
  	vx = ((struct data *)data)->x;
        vy = ((struct data *)data)->y;
   	err = ((struct data *)data)->err;

	      a = gsl_vector_get(x,0);
	      b = gsl_vector_get(x,1);

		for(i=1; i<n; i++)
		{
	     		t = vx[i];
	     		Yi = b*pow(t,a);
	      		fset =  (Yi - vy[i])/err[i];
	
			if(fset!=fset) 
				fset=0;

			gsl_vector_set(f, i, fset);
		}

	return GSL_SUCCESS;
}



int d_power_law_f(const gsl_vector *x, void *data, gsl_matrix *J)
{
	size_t n, i; 
	double s, t, e1, e2, a, b, *vx, *err;

	n = ((struct data*)data)->n;	
	vx = ((struct data *)data)->x;
	err = ((struct data*)data)->err;

	a = gsl_vector_get(x,0);
	b = gsl_vector_get(x,1);

	for(i = 1; i < n; i++)
	{
		s = err[i];
		t = vx[i]; 
	
		e1 = b*log(t)*pow(t,a);
		e2 = pow(t,a);

		gsl_matrix_set(J,i,0,e1/s);
		gsl_matrix_set(J,i,1,e2/s);
	}	

	return GSL_SUCCESS;
}



int fd_power_law_f(const gsl_vector *x, void *data, gsl_vector * f, gsl_matrix *J)
{

		power_law_f(x,data,f);
		d_power_law_f(x,data,J);

	return GSL_SUCCESS;
}



double* best_fit_power_law(double *x, double *y, double *err, int tot, double *guess)
{
	int i=0, n_eff=0;
	double a, b, *vec, *x_eff, *y_eff, *e_eff;
	struct data dat;
	struct parameters par;

	a = guess[0];
	b = guess[1]; 

	fprintf(stdout, "\nFinding best fit values using nlls fit, guess values a:%e, b:%e.\n", a, b);
	
	x_eff = (double*) calloc(1, sizeof(double));
	y_eff = (double*) calloc(1, sizeof(double));
	e_eff = (double*) calloc(1, sizeof(double));

	for(i=0; i<tot; i++)
	{
		if(x[i] != 0.)
		{
			n_eff++;
			x_eff = realloc(x_eff, n_eff*sieof(double));
			y_eff = realloc(y_eff, n_eff*sieof(double));
			e_eff = realloc(e_eff, n_eff*sieof(double));
			x_eff[n_eff] = x[i];
			y_eff[n_eff] = y[i];
			e_eff[n_eff] = err[i];
		}
	}


	dat.n = (size_t) n_eff;
	dat.x = x_eff;
	dat.y = y_eff;
	dat.err = e_eff; 

		par.n = 2;
		par.guess_p = (double *) calloc(par.n, sizeof(double));
		par.guess_p[0] = a;
		par.guess_p[1] = b;
		par.fitted_p = gsl_vector_alloc(par.n);

			/* Initialize function */
		gsl_multifit_function_fdf f;
       		f.f = &power_law_f;
     		f.df = &d_power_law_f;
     		f.fdf = &fd_power_law_f;
       		f.n = dat.n;
       		f.p = par.n;
       		f.params = &dat;

			/* Do the fit */
			par.fitted_p = least_square_nl_fit(dat, par, f);

			/* Set the correctly fitted parameters */
			a = gsl_vector_get(par.fitted_p,0);
			b = gsl_vector_get(par.fitted_p,1);
			fprintf(stdout, "The best fit parameters for this distribution are alpha: %e b:%e \n", a, b);

		vec = (double*) calloc(2,sizeof(double));
		vec[0] = a;
		vec[1] = b;

	return vec; 
}
