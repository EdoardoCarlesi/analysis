/* navarro frenk white density profile utilities */
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

#include "../libgen/vars.h"
#include "../libstat/statistics.h"

/* Navarro Frenk White profile */
double nfw(double r, double r_s, double rho_0){
double rrs = r/r_s;
double rho_r = rho_0*(1./(rrs*(1.+rrs)*(1.+rrs)));
return rho_r;
}

/* NFW derivatives */
double nfw_drs(double r, double rs, double rho_0){
double res = (3.*rs*rs*rho_0)/(r*(rs+r)*(rs+r)) - (2.*rs*rs*rs)/(r*(rs+r)*(rs+r)*(rs+r));
return res;
}

/* NFW derivatives */
double nfw_drho0(double r, double rs, double rho_0){
double res = 1./((r/rs)*(1.+r/rs)*(1.+r/rs));
return res;
}

int nfw_f(const gsl_vector *x, void *data, gsl_vector *f){
size_t n = ((struct data *)data)->n;
double *vx = ((struct data *)data)->x;
double *vy = ((struct data *)data)->y;
double *err = ((struct data *)data)->err;

	double rs = gsl_vector_get(x,0);
	double rho0 = gsl_vector_get(x,1);
	//fprintf(stderr, "nfw - rs:%lf, rho0:%lf \n", rs, rho0);
	size_t i;
	for(i=1; i<n; i++){
	double t = vx[i];
	double rrs = t/rs;
	double Yi = rho0*(1./(rrs*(1.+rrs)*(1.+rrs)));
	double fset =  (Yi - vy[i])/err[i];
//fprintf(stderr, "rs:%lf, rho0:%lf, vx: %lf, vy: %lf, rrs: %lf, Yi: %lf \n", rs,rho0, t, vy[i], rrs, Yi);
	//fprintf(stderr, "%s %lf %lf %lf \n", "fset, Yi, vy, err: ", fset, Yi, vy[i], err[i]);
	
	if(fset!=fset) fset=0;
	gsl_vector_set(f, i, fset);
	//fprintf(stderr, "%s %lf \n", "f_data(): ", f->data[i]);
	}
	return GSL_SUCCESS;
}

int d_nfw_f(const gsl_vector *x, void *data, gsl_matrix *J){
size_t n = ((struct data*)data)->n;	
double *vx = ((struct data *)data)->x;
double *err = ((struct data*)data)->err;

	double rs = gsl_vector_get(x,0);
	double rho0 = gsl_vector_get(x,1);

	//fprintf(stderr, "%s %lf %lf \n", "nfw - rs, rho0: ", rs, rho0);
	size_t i;

	for(i = 1; i < n; i++){
	double s = err[i];
	double r = vx[i]; 
		
	//double sig = 0;//sigma_M(r);

	//fprintf(stderr, "%s %lf %lf %lf \n", "norm, arg, lgn: ", norm, arg, nfw);
	//double da_nfw =  -1./(2*sig*sig)*(2/l_0 - 2*log(t)/l_0)*nfw;
	//double db_nfw = (-1./sig + pow(log(t/l_0),2)*pow(sig,-3))*nfw;
	//double de_nfw = (-1./sig + pow(log(t/l_0),2)*pow(sig,-3))*nfw;

	double drho0_nfw = 1./((r/rs)*(1.+r/rs)*(1.+r/rs));
	double drs_nfw = (3.*rs*rs*rho0)/(r*(rs+r)*(rs+r)) - (2.*rs*rs*rs*rho0)/(r*(rs+r)*(rs+r)*(rs+r));
	
	//double e1 = dl_nfw_distribution(t,l_0,sig);
	//double e2 = ds_nfw_distribution(t,l_0,sig);
	double e1 = drs_nfw; 
	double e2 = drho0_nfw; 
	//fprintf(stderr, "err: %lf, r: %lf, drs: %lf, drho0: %lf\n", s, r, e1, e2);
	gsl_matrix_set(J,i,0,e1/s);
	gsl_matrix_set(J,i,1,e2/s);
	}
	return GSL_SUCCESS;
}

int fd_nfw_f(const gsl_vector *x, void *data, gsl_vector * f, gsl_matrix *J){

	nfw_f(x,data,f);
	d_nfw_f(x,data,J);

	return GSL_SUCCESS;
}

double* best_fit_nfw(double rho0, double rs, int nBins, double *array_data_x, double *array_data_y, double *y_err){
//fprintf(stderr, "best_fit_nfw()\n");
	/* Best Fit Values */
 struct data dat;
dat.n = (size_t) nBins;
dat.x = array_data_x;
dat.y = array_data_y;
dat.err = y_err; 

NFW.err = dat.err;

/* First guess initialization */
struct parameters par;
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
//fprintf(stderr, "Finding best fit values using a non linear least square fit.\n");

/* Do the fit */
par.fitted_p = least_square_nl_fit(dat, par, f);
//	fprintf(stderr, "Finding best fit values using a non linear least square fit.\n");

/* Set the correctly fitted parameters */
rs = gsl_vector_get(par.fitted_p,0);
rho0 = gsl_vector_get(par.fitted_p,1);

//fprintf(stderr, "The best fit parameters for this distribution are: \n");
//fprintf(stderr, "A: %lf B: %lf \n", rs, rho0);

/* Return the best fit parameters*/
double *params;
params = (double*) calloc(2, sizeof(double));
params[0]=rs;
params[1]=rho0; 

NFW.rs=rs;
NFW.rho0=rho0;

return params;
}

void print_nfw(){
FILE * out_nfw = fopen("nfw_test.dat", "w");

int k=0;
for(k=0; k<NFW.bins; k++)
{
double rr=NFW.radius[k];
double od=NFW.overd[k];
double rs=NFW.rs;
double r0=NFW.rho0;

double nnffww = nfw(rr, rs, r0);
fprintf(out_nfw, "%lf %lf %lf %lf\n", rr, nnffww, od, NFW.err[k]);
}
fclose(out_nfw);
}

void fit_and_store_nfw_parameters(int massCut){
fprintf(stderr, "fit_and_store_nfw_parameters().\n");
	double cc, mm, rr, rho0, rs, r2;
	int bins; int skip; int k; int j;
	for(k=0; k<massCut; k++){
	rr = haloes[k].Rvir;
	cc = haloes[k].c; 
	mm = haloes[k].Mvir;
	bins = haloes[k].n_bins;
	skip = haloes[k].neg_r_bins;
	rho0 = cc*cc*cc*200/(3*(log(1+cc)-(cc/(1+cc))));
	rs = haloes[k].r2;
	r2 = haloes[k].r2;
//fprintf(stderr, "fit_and_store_nfw_parameters().\n");
best_fit_nfw(rho0, rs, bins, haloes[k].radius, haloes[k].rho, haloes[k].err);
//fprintf(stderr, "fit_and_store_nfw_parameters().\n");
	//best_fit_nfw(rho0, rs, bins, haloes[k].radius, haloes[k].over_rho, haloes[k].over_err);

	double *y_th; y_th = (double*) calloc(haloes[k].n_bins,sizeof(double));
	haloes[k].rs_nfw = NFW.rs;
	haloes[k].rho0_nfw = NFW.rho0;
//fprintf(stderr, "fit_and_store_nfw_parameters().\n");
	for(j=0; j<haloes[k].n_bins; j++) y_th[j] = nfw(haloes[k].radius[j], haloes[k].rs_nfw, haloes[k].rho0_nfw);
	//double chisq = chi_square(y_th,haloes[k].over_rho,haloes[k].over_err,bins,skip);
	double chisq = chi_square(y_th,haloes[k].rho,haloes[k].err,bins,skip);
//fprintf(stderr, "fit_and_store_nfw_parameters().\n");
	haloes[k].chi_nfw = chisq/(haloes[k].n_bins-haloes[k].neg_r_bins-1);
//fprintf(stderr, "fit_and_store_nfw(), chi:%lf, rs: %lf, r2: %lf,  rho0: %e,  NFW.rho0: %e,  c: %lf\n", chisq/(haloes[k].n_bins-haloes[k].neg_r_bins-1),NFW.rs, r2, rho0, NFW.rho0, haloes[k].Rvir/haloes[k].rs_nfw);
//fprintf(stderr, "fit_and_store_nfw(), chi:%lf, rs:%lf, rho0:%lf\n", chisq, NFW.rs, NFW.rho0);
}
}

void fit_and_store_nfw_parameters_from_list(int massCut,int *list, int max){
fprintf(stderr, "fit_and_store_nfw_parameters_from_list().\n");
	double cc, mm, rr, rho0, rs;
	int bins; int skip; int k; int j; int i;
	for(j=0; j<massCut; j++){
	k=list[j];
	//k=j; //list[j];
	rr = haloes[k].Rvir;
	cc = haloes[k].Rvir/haloes[k].r2;
	mm = haloes[k].Mvir;
	bins = haloes[k].n_bins;
	skip = haloes[k].neg_r_bins;
	rho0 = cc*cc*cc*200/(3*(log(1+cc)-(cc/(1+cc))));
	rs = haloes[k].r2;
//fprintf(stderr, "j:%d) k:%d) R: %lf \n", j,k,rr);
	if(k<max-1) {
	best_fit_nfw(rho0, rs, bins, haloes[k].radius, haloes[k].rho, haloes[k].err);
	haloes[k].rs_nfw = NFW.rs;
	haloes[k].rho0_nfw = NFW.rho0;
//fprintf(stderr, "fit_and_store_nfw(), %d, %d) rs:%lf, rho0:%lf\n", j, k, rs, rho0);
	double *y_th; y_th = (double*) calloc(haloes[k].n_bins,sizeof(double));
	for(i=0; i<bins; i++) y_th[i] = nfw(haloes[k].radius[i], haloes[k].rs_nfw, haloes[k].rho0_nfw);
	double chisq = chi_square(y_th,haloes[k].rho,haloes[k].err,bins,skip);
	haloes[k].chi_nfw = chisq /(haloes[k].n_bins-haloes[k].neg_r_bins-1);
}
//fprintf(stderr, "fit_and_store_nfw(), chi:%lf, rs:%lf, rho0:%lf\n", Chi2.chi, NFW.rs, NFW.rho0);
}
}

