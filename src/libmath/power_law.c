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

#include "power_law.h"
#include "statistics.h"


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
	
				if(fset!=fset) fset=0;

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
	double a, b, *vec;
	struct data dat;
	struct parameters par;

	a = guess[0];
	b = guess[1]; 

	fprintf(stderr, "Finding best fit values using a non linear least square fit.\n");
	fprintf(stderr, "Initial values are, a:%lf, b:%lf.\n", a, b);

	dat.n = (size_t) tot;
	dat.x = x;
	dat.y = y;
	dat.err = err; 

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
			fprintf(stderr, "The best fit parameters for this distribution are: \n");
			fprintf(stderr, "alpha: %lf b:%lf \n", a, b);

		vec = (double*) calloc(2,sizeof(double));
		vec[0] = a;
		vec[1] = b;

	return vec; 
}
