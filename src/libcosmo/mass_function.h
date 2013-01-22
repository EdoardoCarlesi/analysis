extern int PK_INDEX;

double M_r(double);
double W_r(double, double);
double R_m(double);

double sigma_integ(double, void*);
double integrate_sigma2(double);
void normalize_sigma8(void);

void store_mf(int);
double mf_normalization(double);

void init_MassFunc(void);
void init_ThMassFunc(void);
void init_sigmas(void);
void init_dlnsdlnm(void);

double sigma_M(double);
double ln_sigma_M(double, void*);
double d_ln_sigma_d_ln_M(double);

void compute_numerical_mass_function(void);
void compute_theoretical_mass_function(int);

void initialize_mass_function_datastore(void);
