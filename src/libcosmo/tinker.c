#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>

#include "../general_def.h"
#include "../libmath/math.h"

#include "cosmo.h"


/*
 * Declare functions
 */
struct tinker_mf T_mf;

double mf_tinker(double, double, double, double, double);
double dA_mf_tinker(double, double, double, double, double);
double da_mf_tinker(double, double, double, double, double);
double db_mf_tinker(double, double, double, double, double);
double de_mf_tinker(double, double, double, double, double);

int mf_tinker_f(const gsl_vector*, void *, gsl_vector*);
int d_mf_tinker_f(const gsl_vector*, void *, gsl_matrix*);
int fd_mf_tinker_f(const gsl_vector*, void *, gsl_vector*, gsl_matrix*);


/*
 * Initialize functions
 */ 
	/* Returns the differential number density of objects for the Tinker 2008 mf*/
double tinker(double M, void *p)
{
	double A, a, b, c, norm, sig;

		A = T_mf.A;
		a = T_mf.a;
		b = T_mf.b;
		c = T_mf.c;

		norm = mf_normalization(M);
		sig = sigmaM(M);

	return norm*A*(pow(sig/b,-a) + 1)*exp(-c/(sig*sig));
}



double integral_tinker(double Mmin, double Mmax)
{
	double result, error;
	gsl_function T;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

		T.function=&tinker;
		T.params=0;
		gsl_integration_qags(&T, Mmin, Mmax, 0, 1e-4, 1000, w, &result, &error);
		gsl_integration_workspace_free (w);

	return result;
}


	/* Same as before, but with user -given parameters */
double mf_tinker(double M, double A, double a, double b, double c)
{
	double norm, sig;

		norm = mf_normalization(M);
		sig = sigmaM(M);

	return norm*A*(pow(sig/b,-a) + 1)*exp(-c/(sig*sig));
}



double dA_mf_tinker(double M, double A, double a, double b, double c)
{
	double norm, sig;

		norm = mf_normalization(M);
		sig = sigmaM(M);

	return norm*(pow(sig/b,-a) + 1)*exp(-c/(sig*sig));
}



double da_mf_tinker(double M, double A, double a, double b, double c)
{
	double norm, sig;

		norm = mf_normalization(M);
		sig = sigmaM(M);

	return norm*A*( pow(b/sig,a)*log(b/sig)*ln_10)*exp(-c/(sig*sig));
}


double db_mf_tinker(double M, double A, double a, double b, double c)
{
	double norm, sig;
	
		norm = mf_normalization(M);
		sig = sigmaM(M);
	
	return norm*A*(a*pow(sig,-a)*pow(b,a-1))*exp(-c/(sig*sig));
}



double dc_mf_tinker(double M, double A, double a, double b, double c)
{
	double norm, sig;

		norm = mf_normalization(M);
		sig = sigmaM(M);
	
	return -norm*A*pow(sig,-2)*(pow(b/sig,a) + 1)*exp(-c/(sig*sig));
}



int mf_tinker_f(const gsl_vector *x, void *data, gsl_vector *f)
{
	size_t n, i;
	double A, a, b, c, t, s, Yi, fset, *vx, *vy, *err;

		n = ((struct data *)data)->n;
		vx = ((struct data *)data)->x;
		vy = ((struct data *)data)->y;
		err = ((struct data *)data)->err;

		A = gsl_vector_get(x,0);
		a = gsl_vector_get(x,1);
		b = gsl_vector_get(x,2);
		c = gsl_vector_get(x,3);

		for(i=0; i<n; i++)
		{
			if(vy[i]!=0 && vx[i]!=0)
			{
				t = vx[i];
				s = err[i];
				
					if(s==0) 
						s = 2*err[i-1];
					
				Yi = mf_tinker(t,A,a,b,c);
				fset = (Yi - vy[i])/s;
	
			if(fset!=fset)
				fset=0;

			gsl_vector_set(f, i, fset);
			}
		}

	return GSL_SUCCESS;
}



int d_mf_tinker_f(const gsl_vector *x, void *data, gsl_matrix *J)
{	
	size_t i, n;
	double A, a, b, c, s, t, e1, e2, e3, e4, *vx, *err;

	n = ((struct data*)data)->n;	
	vx = ((struct data *)data)->x;
	err = ((struct data*)data)->err;

	A = gsl_vector_get(x,0);
	a = gsl_vector_get(x,1);
	b = gsl_vector_get(x,2);
	c = gsl_vector_get(x,3);

	for(i = 1; i < n; i++)
	{
		s = err[i];

			if(s==0) 
				s = err[i-1];

		t = vx[i]; 

		e1 = dA_mf_tinker(t,A,a,b,c);
		e2 = da_mf_tinker(t,A,a,b,c);
		e3 = db_mf_tinker(t,A,a,b,c); 
		e4 = dc_mf_tinker(t,A,a,b,c); 

	gsl_matrix_set(J,i,0,e1/s);
	gsl_matrix_set(J,i,1,e2/s);
	gsl_matrix_set(J,i,2,e3/s);
	gsl_matrix_set(J,i,3,e4/s);
	}

	return GSL_SUCCESS;
}



int fd_mf_tinker_f(const gsl_vector *x, void *data, gsl_vector * f, gsl_matrix *J)
{

	mf_tinker_f(x,data,f);
	d_mf_tinker_f(x,data,J);

	return GSL_SUCCESS;
}



void best_fit_mf_tinker(double *x_array, double* y_array, double* y_err, int nBins)
{
	double A, a, b, c;
	struct data dat;
	struct parameters par;

	A = T_mf.A;
	a = T_mf.a;
	b = T_mf.b;
	c = T_mf.c;

	fprintf(stdout, "Finding best fit values for the Tinker mass function using a non linear least square fit.\n");

		dat.n = (size_t) nBins;
		dat.x = x_array;
		dat.y = y_array;
		dat.err = y_err;

		/* First guess initialization */
		par.n = 4;
		par.guess_p = (double *) calloc(par.n, sizeof(double));
		par.guess_p[0] = A;
		par.guess_p[1] = a;
		par.guess_p[2] = b;
		par.guess_p[3] = c;
		par.fitted_p = gsl_vector_alloc(par.n);

		/* Initialize function */
		gsl_multifit_function_fdf f;
       		f.f = &mf_tinker_f;
       		f.df = &d_mf_tinker_f;
       		f.fdf = &fd_mf_tinker_f;
       		f.n = dat.n;
       		f.p = par.n;
       		f.params = &dat;

			/* Do the fit */
			par.fitted_p = least_square_nl_fit(dat, par, f);

			/* Set the correctly fitted parameters */
			A = gsl_vector_get(par.fitted_p,0);
			a = gsl_vector_get(par.fitted_p,1);
			b = gsl_vector_get(par.fitted_p,2);
			c = gsl_vector_get(par.fitted_p,3);

		fprintf(stdout, "The best fit parameters for this distribution are: \n");
		fprintf(stdout, "A: %e a: %lf b: %lf c: %lf \n", A, a, b, c);

		T_mf.A = A;
		T_mf.a = a;
		T_mf.b = b;
		T_mf.c = c;

}
