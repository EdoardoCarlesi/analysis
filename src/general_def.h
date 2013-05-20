#include <stdint.h>

#define WARNING(str1, str2) fprintf(stderr, "\n\t\tWARNING! %s: %s.\n", str1, str2)
#define ERROR(str1, str2) fprintf(stderr, "\n\t\tERROR! %s: %s.\n", str1, str2)
#define INFO_MSG(str) fprintf(stdout, "\n%s.\n", str)
#define F_PRINT(str, num) fprintf(stdout, "%s %e.\n", str, num)
#define D_PRINT(str, num) fprintf(stdout, "%s %d.\n", str, num)

	// How to bin the cosmic web
#define BIN_SIZE 15
	// Factors that modify the binning
#define F_MAX 1.0001
#define F_MIN 0.9999
	// Number of bins for subhaloes is multiplied by this
#define F_SUB 0.8
	// Subhalo distribution parameters
#define SUB_MIN 5
#define RMIN 0.2
#define RMAX 1.1
	// Density profile parameters
#define BIN_PROFILE 15
	// Halo density (dm, gas, Ix) profiles will start from  2 * this fraction of Rvir, and gas fraction from 1
#define Rvir_frac_min 0.005

#ifdef WITH_MPI
#define TASK_INFO_MSG(task, str) fprintf(stdout, "\nTask=%d, %s.\n", task, str)
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
	int c_web_size;

	double totGasMassInHalo;
	double totHaloMass;
	double totSubMass;
	double totDarkMass;

	double box_size;
	double zStart;
	
	int n_threshold; // Threshold can be mass or particle number per halo
	int n_virialized;
	int n_concentration;
	int n_spin;
	int n_all;
	int n_haloes;
	int n_haloes_size; // Used for MPI communication
	int n_cweb_type[4];

	int n_sub_min; // Minimum number of subhaloes per halo
	int n_sub_threshold;
	
	int n_bins;
	int n_bins_th;
	int r_bins;

	int use_n_min; // If =1 using n_min as threshold criterion, otherwise using mass_min
	int use_criterion; // 0 use all haloes, 1 use mass/number, 2 use spin, 3 use concentration, 
			   // 4 use virialization, 5 use all combined criteria
	int n_min;
	double mass_min;
	int n_haloes_to_use; // Use the first N haloes
	
	// Conditions in halo subset to be used
	int use_none;
	int use_mass;
	int use_spin;
	int use_vir;
	int use_conc;
	int use_all;
	int use_sub;
	int use_web;
	int use_web_type; // 0, 1, 2, 3

	double l_web;
	double rho_c;
	double rho_0;
	double Rmin;
	double Rmax;
	double Mmin;
	double Mmax;
	double Mprint;

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

} *MassFunc, *VelFunc, *ThMassFunc, *GasFunc, *NoGasFunc; 


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

	// 0 void, 1 sheet, 2 filament, 3 node
	int web_type[4];

	// These variables are set to 1 if the halo satisfies the condition
	int vir;
	int spin;
	int conc;
	int mass; 
	int all;

	float X[3];
	float V[3];
	float L[3];
	float Ea[3]; 
	float a[3];

	double Mvir;
	float Rvir;
	float Rmax;
	float Vmax;
	double Ekin;
	double Epot;
	float th_vir;
	float abs_th_vir;
	float AngMom;
	float Jcirc;
	float ecc;
	float lambda;
	float lambdaE;

	float shape;
	float triax;
	float z_form;

	// NFW parameters
	float c;
	float rho0;
	float r2;
	float c_nfw;
	float rho0_nfw;
	float rs_nfw;
	
	struct
	{
		float c; // NFW only
		float rs;
		float rho0;
		float ix0; // X-ray surface brightness parameter
		float beta; // X-ray surface brightness parameter / beta model
		float chi;
		float gof;
		float per;
	} fit_nfw, fit_poly, fit_beta, fit_IX;

	struct
	{
		float x[BIN_PROFILE];
		float y[BIN_PROFILE];
	} f_gas, nfw, rho_gas, i_x, hydro_m, T;

	int *npart;
	float *radius;
	double *mass_r;
	float *rho;
	float *err;

#ifdef GAS
	struct
	{
		int N;
		float M;
		float vir;
		float Ekin;
		float Epot;
#ifdef EXTRA_GAS
		float a[3];
		float Ea[3];
		float X[3];
		float V[3];
#endif
	} gas, dm;

	struct
	{
		float lambda;
		float lambdaE;
		float b_fraction;
		double Cum_u;
		float T_mw; // Mass-weighted temperature
		float T_0; // Central cluster temperature
		float I_X0; // X ray thermal emission
		float shape;
		float triax;
		
		// Alignment angle of the cluster/gas major axis
		float gas_dm_costh;

		struct
		{
			float cm;
			float shape;
			float triax;
			float lambda;
		} diff;

		double *m;
		double *u;
	
		// Extra info
		float *rho;
		float *frac;
		float *i_x;
		float *hydro_m;
		float *T;
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
	char *c_web_file;
	char *c_web_gas_file;
	
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
	double z;

	// Axis alignment
	double *R;
	double *Th_c;
	double *Th_p;
	int *N_pairs;

	// The masses stored here are in different than massfunc, since include the mass/npart threshold
	double *mass;
	
	// Number of entries per mass bin 
	int *n_entry; 	

	// Parameters that correlate with the mass
	double *vel;
	double *conc;

	struct
	{
		double s0;
		double t0;

		double *shape;
		double *triax;
		double *lambda;
		double *beta;
		double *virial;
		double *ekin;

		double *s;
		double *p_s;
		double *t;
		double *p_t;
		double *l;
		double *p_l;
		double *b;
		double *p_b;

	} halo, dm, gas, diff;

	struct
	{
		double *chi;	
		double *per;	
		double *gof;
		double *p_chi;	
		double *p_gof;	
		double *p_per;	

	} fit_nfw, fit_king, p_fit_nfw;

	struct
	{
		int n[BIN_PROFILE];
		double x[BIN_PROFILE];
		double y[BIN_PROFILE];

	} f_gas, nfw, rho_gas, i_x;

	// Distributions
	double *c;
	double *p_c;
	
	// Subhaloes
	double avgSub;
	double *sub;
	double *avg_sub;
	double *p_avg_sub;

	// Gas displacement and alignment
	double *cm;
	double *p_cm;
	double *gas_dm_cth;
	double *p_gas_dm_cth;

	// Best fit parameters
	double beta0; // Gas density
	double l_0;
	double l_sig;
	double c_0;
	double c_beta;
	double vel_0;
	double vel_beta;

	// Best fit parameters for the M-T relation
	double T0;
	double alpha;

	// n(T)
	double *T;
	double *n_T;

	// Subhalo stuff
	double *n_r_sub;
	double *cum_n_r_sub;

	double *r_sub;
	double avgMass;
	double *mass_sub;

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

	double *n_vel_sub;
	double *cum_n_vel_sub;
	double *vel_sub;
	double *vel_sub_r;
	double *p_vel_sub;
	double vel_0_sub;
	
	double *ecc;

	double *n_sub;
	double *cum_n_sub;
	int *n_ecc;

	double *gas_T;
	double *gas_u;
	double *gas_fraction;
	double *gas_dm_costh;
	double *gas_diff_cm;

} *HaloProperties;


void allocate_mass_function(struct mass_function *, int);
