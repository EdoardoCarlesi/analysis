#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>

struct data
{
	size_t n;
	double *err;
	double *x;
	double *y;
};

struct parameters
{
	size_t n;
	double *guess_p;
	gsl_vector *fitted_p;
};

struct chi2_reduced
{
	double chi;
	int ntot;
	int bins;
	int *outcomes;
	double *binned_chi;
	double *chi2s;
} Chi2;

double mean(double*, int);
double log_mean(double*, int);

double sigma(double*, double, int);
double log_sigma(double*, int);
double* numerical_sigma(int, double*, double*);

gsl_vector* least_square_nl_fit(struct data dat, struct parameters, gsl_multifit_function_fdf);

void print_state (size_t, gsl_multifit_fdfsolver*);

double chi_square(double*, double*, double*, int, int);
void calculate_and_sort_chi2();
