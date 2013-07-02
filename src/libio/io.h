void set_halo_url(void);
void get_halo_files_urls(void);

int get_lines(FILE *, char *url);

void read_and_store_halo_urls(void);

void read_redshift_file(void);
void read_halo_file(void);
void read_profiles_file(void);
void read_pk_snapshots(void);
void init_pks(void);

void print_sub_per_host(int,int,int,double*,int*,int*,double*,double*,double*,double*,int*,double*,int*,double*,double*,double*,double*);

void print_axis_alignment(void);

void print_all_halo_properties_to_one_file(void);
void print_all_haloes(void);
void print_halo_profile(int);
void print_subhalo_only_properties(void);
void print_theoretical_mass_function(void);
void print_numerical_mass_function(void);
void print_halo_best_fit_results(void);

void print_correlation_function(void);
void print_average_profiles(void);

void print_evolution_to_file(void);
void print_number_densities(void);
void print_growth_factor(void);

