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

#include "nfw.h"
#include "../general_variables.h"
#include "../libmath/statistics.h"


struct nfw NFW;


double nfw(double r, double r_s, double rho_0)
{
	double rrs; 

		rrs = r/r_s;

	return rho_0*(1./(rrs*(1.+rrs)*(1.+rrs)));
}



double nfw_drs(double r, double rs, double rho_0)
{
	return (3.*rs*rs*rho_0)/(r*(rs+r)*(rs+r)) - (2.*rs*rs*rs)/(r*(rs+r)*(rs+r)*(rs+r));
}



double nfw_drho0(double r, double rs, double rho_0)
{
	return 1./((r/rs)*(1.+r/rs)*(1.+r/rs));
}



int nfw_f(const gsl_vector *x, void *data, gsl_vector *f)
{
	size_t i, n;
	double rs, rho0, t, rrs, Yi, fset, *vx, *vy, *err;

	n = ((struct data *)data)->n;
	vx = ((struct data *)data)->x;
	vy = ((struct data *)data)->y;
	err = ((struct data *)data)->err;

	rs = gsl_vector_get(x,0);
	rho0 = gsl_vector_get(x,1);

	for(i=1; i<n; i++)
	{
		t = vx[i];
		rrs = t/rs;
		Yi = rho0*(1./(rrs*(1.+rrs)*(1.+rrs)));
		fset =  (Yi - vy[i])/err[i];
	
		if(fset!=fset) fset=0;

	gsl_vector_set(f, i, fset);
	}

	return GSL_SUCCESS;
}



int d_nfw_f(const gsl_vector *x, void *data, gsl_matrix *J)
{
	size_t i, n;
	double rs, rho0, s, r, e1, e2, *vx, *err;
		
	n = ((struct data*)data)->n;	
	vx = ((struct data *)data)->x;
	err = ((struct data*)data)->err;

	rs = gsl_vector_get(x,0);
	rho0 = gsl_vector_get(x,1);

	for(i = 1; i < n; i++)
	{
		s = err[i];
		r = vx[i]; 
	
		e1 = (3.*rs*rs*rho0)/(r*(rs+r)*(rs+r)) - (2.*rs*rs*rs*rho0)/(r*(rs+r)*(rs+r)*(rs+r));
		e2 = 1./((r/rs)*(1.+r/rs)*(1.+r/rs));

	gsl_matrix_set(J,i,0,e1/s);
	gsl_matrix_set(J,i,1,e2/s);
	}

	return GSL_SUCCESS;
}



int fd_nfw_f(const gsl_vector *x, void *data, gsl_vector * f, gsl_matrix *J)
{

	nfw_f(x,data,f);
	d_nfw_f(x,data,J);

	return GSL_SUCCESS;
}



double* best_fit_nfw(double rho0, double rs, int nBins, double *array_data_x, double *array_data_y, double *y_err)
{
	double *params;
	struct data dat;
	struct parameters par;

	dat.n = (size_t) nBins;
	dat.x = array_data_x;
	dat.y = array_data_y;
	dat.err = y_err; 

	NFW.err = dat.err;

		/* First guess initialization */
	par.n = 2;
	par.guess_p = (double *) calloc(par.n, sizeof(double));
	par.guess_p[0] = rs;
	par.guess_p[1] = rho0;
	par.fitted_p = gsl_vector_alloc(par.n);

		/* Initialize function */
	gsl_multifit_function_fdf f;
       	f.f = &nfw_f;
       	f.df = &d_nfw_f;
       	f.fdf = &fd_nfw_f;
       	f.n = dat.n;
       	f.p = par.n;
       	f.params = &dat;


		/* Do the fit */
	par.fitted_p = least_square_nl_fit(dat, par, f);

		/* Set the correctly fitted parameters */
	rs = gsl_vector_get(par.fitted_p,0);
	rho0 = gsl_vector_get(par.fitted_p,1);


			/* Return the best fit parameters*/
		params = (double*) calloc(2, sizeof(double));
		params[0]=rs;
		params[1]=rho0; 

		NFW.rs=rs;
		NFW.rho0=rho0;

	return params;
}



void fit_and_store_nfw_parameters(int massCut)
{
	fprintf(stderr, "fit_and_store_nfw_parameters(). Using the first %d haloes.\n", massCut);

	double cc, mm, rr, rho0, rs, r2, chisq, *y_th; 
	int bins, skip, k=0, j=0;

		for(k=0; k<massCut; k++)
	{
			rr = haloes[k].Rvir;
			cc = haloes[k].c; 
			mm = haloes[k].Mvir;
			bins = haloes[k].n_bins;
			skip = haloes[k].neg_r_bins;
			rho0 = cc*cc*cc*200/(3*(log(1+cc)-(cc/(1+cc))));
			rs = haloes[k].r2;
			r2 = haloes[k].r2;

				best_fit_nfw(rho0, rs, bins, haloes[k].radius, haloes[k].rho, haloes[k].err);

			y_th = (double*) calloc(haloes[k].n_bins,sizeof(double));
			haloes[k].rs_nfw = NFW.rs;
			haloes[k].rho0_nfw = NFW.rho0;

			for(j=0; j<haloes[k].n_bins; j++) 
				y_th[j] = nfw(haloes[k].radius[j], haloes[k].rs_nfw, haloes[k].rho0_nfw);
	
		chisq = chi_square(y_th,haloes[k].rho,haloes[k].err,bins,skip);
		haloes[k].chi_nfw = chisq/(haloes[k].n_bins-haloes[k].neg_r_bins-1);
	}
}



void fit_and_store_nfw_parameters_from_list(int massCut, int *list, int max)
{
	fprintf(stderr, "fit_and_store_nfw_parameters_from_list().\n");

	double cc, mm, rr, rho0, rs, chisq, *y_th;
	int bins, skip, k=0, j=0, i=0;

	for(j=0; j<massCut; j++)
		{
			k=list[j];
			rr = haloes[k].Rvir;
			cc = haloes[k].Rvir/haloes[k].r2;
			mm = haloes[k].Mvir;
			bins = haloes[k].n_bins;
			skip = haloes[k].neg_r_bins;
			rho0 = cc*cc*cc*200/(3*(log(1+cc)-(cc/(1+cc))));
			rs = haloes[k].r2;
		
		if(k<max-1) 
			{
				best_fit_nfw(rho0, rs, bins, haloes[k].radius, haloes[k].rho, haloes[k].err);
				haloes[k].rs_nfw = NFW.rs;
				haloes[k].rho0_nfw = NFW.rho0;

				y_th = (double*) calloc(haloes[k].n_bins,sizeof(double));

			for(i=0; i<bins; i++) 
				y_th[i] = nfw(haloes[k].radius[i], haloes[k].rs_nfw, haloes[k].rho0_nfw);

			chisq = chi_square(y_th,haloes[k].rho,haloes[k].err,bins,skip);
			haloes[k].chi_nfw = chisq /(haloes[k].n_bins-haloes[k].neg_r_bins-1);
			}
		}
}
