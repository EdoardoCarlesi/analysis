extern int HALO_INDEX;

void fit_and_store_nfw_parameters(void);
void average_nfw_profile(void);
double nfw(double, double, double);

#ifdef GAS
void fit_and_store_gas_parameters(void);
void average_gas_profiles(void);

	// Polytropic (Adiabatic) Temperature profile
double polytropic_T(double, double, double);

	// X-Ray surface brightness based on isothermal Beta model for the gas density profile
double rhoBeta(double, double, double, double);
double I_X(double, double, double, double);

	// King profile
double king(double, double, double);
#endif

void sort_axis_alignment(void);
void compute_subhalo_properties(void);
void compute_halo_properties(void);
void compute_halo_and_subhalo_properties(void);

void initialize_halo_storage(void);
void initialize_halo_properties_structure(void);

void find_substructure(void);

void free_halo_profiles(void);

void free_subhalo_properties(void);
void free_halo_properties(void);
void free_halo_profiles(void);
