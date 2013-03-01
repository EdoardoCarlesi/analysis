extern int HALO_INDEX;

void fit_and_store_nfw_parameters(void);

void compute_subhalo_properties(void);
void compute_halo_properties(void);
void compute_halo_and_subhalo_properties(void);

void initialize_halo_storage(void);
void initialize_halo_properties_structure(void);

void find_substructure(void);

void free_subhalo_properties(void);
void free_halo_properties(void);
void free_halo_profiles(void);
