#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

int power_law_f(const gsl_vector*, void *, gsl_vector*);
int d_power_law_f(const gsl_vector*, void *, gsl_matrix*);
int fd_power_law_f(const gsl_vector*, void *, gsl_vector*, gsl_matrix*);

double* best_fit_power_law(double *, double*, double *, int, double*);
