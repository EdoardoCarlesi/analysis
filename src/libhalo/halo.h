extern int HALO_INDEX;

void initialize_halo_properties_structure(void);

void fit_and_store_nfw_parameters(void);

void compute_subhalo_properties(void);
void compute_halo_properties(void);
void compute_halo_and_subhalo_statistics(void);

void initialize_halo_evolution_structures(void);
void initialize_halo_storage(void);

void initialize_subhalo_properties_structure(void);
void find_substructure(void);

void free_subhalo_properties(void);
void free_halo_properties(void);
void free_halo_profiles(void);
