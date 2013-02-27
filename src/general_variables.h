#ifndef WITH_MPI
#define ThisTask 0
#define NTask 1
#endif

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

	double rho_c;
	double rho_0;
	double Rmin;
	double Rmax;
	double Mmin;
	double Mmax;

		struct
		{
			double X[2];	// X[0] min, X[1] max
			double Y[2];
			double Z[2];
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

	// When loading cosmology parameter from file, store them here
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

} MassFunc, ThMassFunc, *MassFuncZ, *ThMassFuncZ;


extern struct num_density
{
	int npts;
	double zMax;
	double zMin;

	double *z;
	double *n_th;
	double *n_num;

} NumDen; 


extern struct nfw
{
	int bins;
	double rs;
	double rho0;
	double *radius;
	double *overd;
	double *profile;
	double *err;
	double *chi2r;

} NFW;


extern struct merger_tree
{
	int n_haloes;
	int n_z;
	int *id;
	double *z;
	double *M;
	double *c;

} *MergerTree;


extern struct halo
{
	long int id;
	int n_part;
	int n_satellites;
        int n_bins;
	int neg_r_bins;
	int host;

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

	double c;
	double rho0;
	double r2;
	double c_nfw;
	double rho0_nfw;
	double rs_nfw;
	double chi_nfw;

	int *id_satellites;

	double *radius;
	double *rho;
	double *over_rho;
	double *err;
	double *over_err;
	double *err_dn;
	double *bin;

#ifdef GAS
	int N_dm;
	int N_gas;
	double M_dm;
	double M_gas;

	double a_gas[3];
	double Ea_gas[3];

	double lambda_gas;
	double lambdaE_gas;
	double b_gas;
	double c_gas;
	double Ekin_gas;
	double Epot_gas;
	double b_fraction;
	double Cum_u_gas;
	double T_gas;

	double *m_gas;
	double *u_gas;

#ifdef EXTRA_GAS
	double X_dm[3];
	double V_dm[3];
	double X_gas[3];
	double V_gas[3];
#endif
#endif  // GAS
} *Haloes, *SubHaloes, **pHaloes, **pSubHaloes;


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
} SubHaloZ, HaloZ, *HaloProperties, *SubHaloProperties;
