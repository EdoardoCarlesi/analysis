#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

extern struct tinker_mf
{
	double *norm;
	double *sig;
	double A;
	double a;
	double b;
	double c;
} T_mf;

double tinker(double, void*);
double integral_tinker(double, double);
double mf_tinker(double, double, double, double, double);

double dA_mf_tinker(double, double, double, double, double);
double da_mf_tinker(double, double, double, double, double);
double db_mf_tinker(double, double, double, double, double);
double de_mf_tinker(double, double, double, double, double);

int mf_tinker_f(const gsl_vector*, void *, gsl_vector*);
int d_mf_tinker_f(const gsl_vector*, void *, gsl_matrix*);
int fd_mf_tinker_f(const gsl_vector*, void *, gsl_vector*, gsl_matrix*);

double* best_fit_mf_tinker(double *, double*, double*, int);
