#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#ifndef GENERAL_VARIABLES_H
#define GENERAL_VARIABLES_H
#endif

#ifndef GENERAL_VARIABLES_H
#include "../general_variables.h"
#endif

extern struct nfw
{
	int bins;
	double *radius;
	double *overd;
	double *err;
	double rs;
	double rho0;
	double *chi2r;
} NFW;

double nfw(double, double, double);
double nfw_drs(double, double, double);
double nfw_drho0(double, double, double);

int nfw_f(const gsl_vector*, void *, gsl_vector*);
int d_nfw_f(const gsl_vector*, void *, gsl_matrix*);
int fd_nfw_f(const gsl_vector*, void *, gsl_vector*, gsl_matrix*);

double* best_fit_nfw(double, double, int, double*, double*, double*);

void fit_and_store_nfw_parameters_from_list(int*, int);
void fit_and_store_nfw_parameters();
