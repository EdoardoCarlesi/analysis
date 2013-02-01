void initialize_halo_properties_structure(void);

void free_halo_properties(void);
void free_halo_profiles(int);

void sort_axis_alignement(void);
void sort_lambda(void);
void sort_concentration(void);
void sort_radial_velocity(void);
void sort_shape_and_triaxiality(void);

#ifdef GAS
void sort_gas_fraction(void);
void sort_and_fit_mass_temperature_relation(void);
#endif

void compute_halo_properties(void);
