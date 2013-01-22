/* nfw.h */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

double nfw(double, double, double);

double nfw_drs(double, double, double);

double nfw_drho0(double, double, double);

int nfw_f(const gsl_vector*, void *, gsl_vector*);

int d_nfw_f(const gsl_vector*, void *, gsl_matrix*);

int fd_nfw_f(const gsl_vector*, void *, gsl_vector*, gsl_matrix*);

void print_nfw();

double* best_fit_nfw(double, double, int, double*, double*, double*);
