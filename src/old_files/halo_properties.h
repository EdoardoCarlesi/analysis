void pinitialize_halo_properties_structure();

void pfree_halo_properties();

void pfit_concentration_relations();

void psort_axis_alignement();
void psort_lambda();
void psort_concentration();
void psort_radial_velocity();
void psort_shape_and_triaxiality();

#ifdef GAS
void sort_gas_fraction();
void sort_and_fit_mass_temperature_relation();
#endif

void pcompute_halo_properties(void);
void pcompute_numerical_mass_function(void);
