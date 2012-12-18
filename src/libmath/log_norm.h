#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

double lognorm_distribution(double, double, double);
double dl_lognorm_distribution(double, double, double);
double ds_lognorm_distribution(double, double, double);

int lognorm_f(const gsl_vector*, void *, gsl_vector*);
int d_lognorm_f(const gsl_vector*, void *, gsl_matrix*);
int fd_lognorm_f(const gsl_vector*, void *, gsl_vector*, gsl_matrix*);

double* best_fit_lognorm(double*, int, int, double*, double*, double*);
