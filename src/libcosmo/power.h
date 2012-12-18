double power_k(double, int);

void compute_growth_factor(void);

void normalize_all_power_spectra_to_sigma8(void);

void compute_correlation_function(int);
void fit_correlation_function(void);
double correlation_r(double, int);
double correlation_integral(double, void*);
