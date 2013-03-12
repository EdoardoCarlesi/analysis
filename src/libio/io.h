void set_halo_url(void);
void get_halo_files_urls(void);

void read_and_store_halo_urls(void);

void read_redshift_file(void);
void read_halo_file(void);
void read_profiles_file(void);
int get_lines(FILE *, char *url);
void read_pk_snapshots(void);
void init_pks(void);

void print_axis_alignment(void);

void print_all_halo_properties_to_one_file(void);
void print_all_subhalo_properties_to_one_file(void);
void print_theoretical_mass_function(void);
void print_numerical_mass_function(void);

void print_correlation_function(void);
void print_average_profiles(void);

void print_evolution_to_file(void);
void print_number_densities(void);
void print_growth_factor(void);

