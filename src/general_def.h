#include <stdint.h>

#define WARNING(str1, str2) fprintf(stderr, "\n\t\tWARNING! %s: %s.\n", str1, str2)
#define ERROR(str1, str2) fprintf(stderr, "\n\t\tERROR! %s: %s.\n", str1, str2)
#define INFO_MSG(str) fprintf(stdout, "\n\t%s.\n", str)

#ifdef WITH_MPI
#define TASK_INFO_MSG(task, str) fprintf(stdout, "\n\tTask=%d, %s.\n", task, str)
#else 
#define ThisTask 0
#define NTask 1
#endif

void initialize_internal_variables(char **);
void default_init(void);

void set_halo_selection_criterion(void);
void check_condition_consistency(void);

int n_haloes_per_criterion(void);
int halo_condition(int);


extern struct power_spectrum 
{
	int npts;

	double z;
	double a;

	double *k; 
	double *pk;
	double *sigma;

} *Pks;


extern struct correlation_function
{
	int npts;

	double z;
	double a;

	double *r;
	double *xi_r;
	double *xi_fit;

} Xi, *Xis;


extern struct growth_factor 
{
	int npts;
	double scale_k;

	double *a;
	double *z;
	double *gf;

} GrowthFac;


extern struct general_settings
{
	int fit;	
	int tick;
	int inverse_read;
	int n_part_1D;
	int halo_skip;
	int pk_skip;
	int use_one_pk;
	int cat_number;
	int use_cat;

	double box_size;
	double zStart;
	
	int n_threshold; // Threshold can be mass or particle number per halo
	int n_virialized;
	int n_concentration;
	int n_spin;
	int n_all;
	int n_haloes;
	int n_haloes_size; // Used for MPI communication

	int min_subhaloes;
	int n_subhaloes;
	int n_subhaloes_nmin;
	int n_bins;
	int n_bins_th;
	int r_bins;

	int use_n_min; // If =1 using n_min as threshold criterion, otherwise using mass_min
	int use_criterion; // 0 use all haloes, 1 use mass/number, 2 use spin, 3 use concentration, 
			   // 4 use virialization, 5 use all combined criteria
	int n_min;
	double mass_min;
	int n_haloes_to_use; // Use the first N haloes
	
	int use_none;
	int use_mass;
	int use_spin;
	int use_vir;
	int use_conc;
	int use_all;
	int use_sub;

	double rho_c;
	double rho_0;
	double Rmin;
	double Rmax;
	double Mmin;
	double Mmax;

		struct
		{
			double X[3][2];	// X[0] min, X[1] max
		} box;

#ifndef GAS
	double pMass;
#else
	double dmMass; 
	double gasMass;
	double rho_dm;
	double rho_b;
#endif
} Settings, *pSettings;


extern struct cosmology
{
	double OmegaM;
	double OmegaL;
	double OmegaDE;
#ifdef GAS
	double OmegaB;
	double OmegaDM;
#endif
	double G; 
	double H_0;
	double h;
	double delta_c;
	double sigma8;
	double norm_sigma8;

	double err;
	double virial;
	double spin;

	int npts;
	double *z_hub;
	double *a_hub;
	double *w;
	double *Hubble;	
	
} Cosmo;


extern struct mass_function
{
	int bins;
	int *n_tot;	
	int *n_bin;

	double z;
	double Mmin;
	double Mmax;

	double *mass;
	double *mass_halfstep;
	double *sigma;

	double *ln_mass;
	double *ln_mass_halfstep;
	double *ln_sigma;
	double *dln_sigma_dln_mass;

	double *dn;
	double *n;

	double *err;
	double *err_dn;

} *MassFunc, *ThMassFunc; 


extern struct num_density
{
	int npts;
	double zMax;
	double zMin;

	double *z;
	double *n_th;
	double *n_num;

} NumDen; 


extern struct halo
{
	uint64_t id;
	uint64_t host;
	int n_part;
	int n_satellites;
        int n_bins;
	int neg_r_bins;

	// These variables are set to 1 if the halo satisfies the condition
	int vir;
	int spin;
	int conc;
	int mass; 
	int all;

	double X[3];
	double V[3];
	double L[3];
	double Ea[3]; 
	double a[3];

	double Mvir;
	double Rvir;
	double Rmax;
	double Vmax;
	double Ekin;
	double Epot;
	double th_vir;
	double abs_th_vir;
	double AngMom;
	double Jcirc;
	double ecc;
	double lambda;
	double lambdaE;

	double shape;
	double triax;
	double z_form;

	// NFW parameters
	double c;
	double rho0;
	double r2;
	double c_nfw;
	double rho0_nfw;
	double rs_nfw;
	
	struct
	{
		double c_nfw;
		double rho0_nfw;
		double rs_nfw;
		double chi_nfw;
		double gof_nfw;
		double per_nfw;
	} fit;

	double *radius;
	double *rho;
	double *err;

#ifdef GAS
	struct
	{
		int N;
		int M;
#ifdef EXTRA_GAS
		double X[3];
		double V[3];
#endif
	} gas, dm;

	struct
	{
		double a[3];
		double Ea[3];

		double lambda;
		double lambdaE;
		double Ekin;
		double Epot;
		double b_fraction;
		double Cum_u;
		double T;

		double *m;
		double *u;
	} gas_only;

#endif  // GAS
} *Haloes, **pHaloes;


extern struct sub_structure
{
	int N_host;
	int N_sub;

	struct sub_halo
	{
		uint64_t id;
		int index;
		int host_id;
		int host_index;
	} *sub;

	struct host_halo
	{
		uint64_t id;
		int index;
		int n_sub;
		int *sub_index;	
	} *host;
	
} SubStructure;


extern struct internal_urls
{
	int nCatalogueFiles;
	int nPkFiles;

	char *a_outputs;
	char *output_prefix;
	char *hubble_file;
	
	// Using a single halo/profile/P(k) file for the analysis
	char *halo_file;
	char *profiles_file;
	char *pk_file;

	// Files that store a list of URLS to the actual files
	char *pk_list;
	char *halo_list;
	char *profile_list;
	char *subhalo_list;
	
	// Arrays containing all the available urls for the analysis
	char **urls;
	char **urls_pks;
	char **urls_profiles;
	char **urls_satellites;
	char **urls_particles;

} Urls, *pUrls;


extern struct halo_properties
{
	int z_bins;
	int r_bins;
	int n_bins;
	double avgSub;
	double z;

	double c_0;
	double c_02;
	double c_sig;
	double c_sig2;

	double *mass; // The masses stored here are in different in principle, since include the mass/npart threshold

	double *c;
	double *c_c0;
	double *c_avg;
	double *p_c;
	double *err_p_c;

	double l_0;
	double l_sig;
	double *l;
	double *p_l;
	double *err_p_l;

	int *N_pairs;
	double *R;
	double *Th_c;
	double *Th_p;

	int *n_shape;
	double s0;
	double *shape;
	double *p_shape;
	
	double *radVel;
	double *err_radVel;

	int *n_triax;
	double t0;
	double *triax;
	double *p_triax;

	int *n_r_sub;
	int *cum_n_r_sub;
	double *r_sub;

	int *n_r_sub_subset;
	int *cum_n_r_sub_subset;
	double *r_sub_subset;

	double costh0;
	double *costh;
	double *costh_count;

	double cosphi0;
	double *cosphi;
	double *cosphi_count;
	double *cum_cosphi;

	int *n_vel_sub;
	double *vel_sub;
	double *p_vel_sub;
	double vel_0;

	double avgMass;
	double *mass_sub;
	
	double *ecc;

	int *n_sub;
	int *cum_n_sub;
	int *n_ecc;

#ifdef GAS
	double ln_M0;
	double T_alpha;

	double *gas_T;
	double *gas_u;
	double *gas_fraction;
#endif
} *HaloProperties, *SubHaloProperties;
