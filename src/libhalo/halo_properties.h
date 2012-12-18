void initialize_halo_properties_structure();

void free_halo_properties();

void fit_concentration_relations();

void sort_axis_alignement();
void sort_lambda();
void sort_concentration();
void sort_radial_velocity();
void sort_shape_and_triaxiality();

#ifdef GAS
void sort_gas_fraction();
void sort_and_fit_mass_temperature_relation();
#endif

void compute_halo_properties();
