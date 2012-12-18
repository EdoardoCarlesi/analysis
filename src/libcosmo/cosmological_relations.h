double default_rho0(void);

double H_z(double, void*);
double inv_H_z(double, void*);

double comoving_distance(double, double);
double integrate_comoving_volume(double, double);
double comoving_vol(double, void*);

double mass_temperature(double);
void mass_variance(void);

double convert_u_to_T(double);
