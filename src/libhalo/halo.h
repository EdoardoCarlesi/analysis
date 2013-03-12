extern int HALO_INDEX;

void fit_and_store_nfw_parameters(void);
double nfw(double, double, double);

#ifdef GAS
void average_gas_fraction_profile(void);
	// Polytropic (Adiabatic) Temperature profile
void fit_and_store_polytropic_T_parameters(void);
double polytropic_T(double, double, double);

	// Beta model for the gas distribution
void fit_and_store_rhoBeta_parameters(void);
double rhoBeta(double, double, double);

	// X-Ray surface brightness
void fit_and_store_I_X_parameters(void);
double I_X(double, double, double);

	// King profile
void fit_and_store_king_parameters(void);
double king(double, double, double);
#endif

void compute_subhalo_properties(void);
void compute_halo_properties(void);
void compute_halo_and_subhalo_properties(void);

void initialize_halo_storage(void);
void initialize_halo_properties_structure(void);

void find_substructure(void);

void free_subhalo_properties(void);
void free_halo_properties(void);
void free_halo_profiles(void);
