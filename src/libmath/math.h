#include <gsl/gsl_vector.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>

#define PI 3.141593
#define inv_ln_10 1./log(10)
#define ln_10 log(10)
#define D_H 3000 
#define pow2(a) ((a)*(a))
#define pow3(a) ((a)*(a)*(a))


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


double average(double*, int);
double maximum(double*, int);
double minimum(double*, int);
double nonzero_minimum(double*, int);
void maxima(double*, int, double*, int*, int);
void minima(double*, int, double*, int*, int);

int int_maximum(int*, int);

double mean(double*, int);
double sigma(double*, double, int);
double log_mean(double*, int);
double log_sigma(double*, int);

double* shellsort(double* , int);
int* int_shellsort(int* , int);

int* generate_random_subset(int, int, int*);

double* invert_array(double*, int);

double* log_stepper(double, double, int);
double* lin_stepper(double, double, int);

void cum_bin(int*, int*, int);
void double_cum_bin(double*, double*, int);
void lin_bin(double*, double*, int, int, int*);
void average_bin(double*, double*, double *, double*, double*, int, int);
void median_bin(double*, double*, double *, double*, double*, int, int);

double get_interpolated_value(double*, double*, int, double);

double solid_angle(double, void*);
double integrate_solid_angle(double,double,double,double);

double chi_square(double*, double*, double*, int);
double goodness_of_fit(double*, double*, int);
double percentage_error(double*, double*, int);

gsl_vector* least_square_nl_fit(struct data, struct parameters, gsl_multifit_function_fdf);

double lognorm(double, double, double);
double* best_fit_lognorm(double*, int, int, double*, double*, double*);
double* best_fit_power_law(double *, double*, double *, int, double*);
