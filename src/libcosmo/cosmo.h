extern struct tinker_mf
{
	double *norm;
	double *sig;
	double A;
	double a;
	double b;
	double c;
} T_mf;

extern int PK_INDEX;
extern int MF_INDEX;

double tinker(double, void*);
double integral_tinker(double, double);
void best_fit_mf_tinker(double *, double*, double*, int);

double default_rho0(void);
double mass_temperature(double);
double convert_u_to_T(double);

double* integrate_number_density(double, double);

void compute_growth_factor(void);

void fit_correlation_function(void);

double power_k(double, int);
double correlation_r(double, int);

void get_n_M_z(double);
void compute_correlation_function(int);
void compute_theoretical_mass_function(void);
