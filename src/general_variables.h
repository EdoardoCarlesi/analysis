extern struct power_spectrum 
{
	int n_pk_entries;
	double z;
	double a;
	double lna;
	double *sigma;
	double *k; 
	double *pk;
	char *url;

} *Pks;


extern struct correlation_function
{
	int n_xi_entries;
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
	double *gf_z;
	double *gf_over_a_z;
	//double *ln_gf_a;
	//double *ln_gf_z;
} GrowthFac;


extern struct general_settings
{
	int fit;	
	int tick;
	int inverse_read;
	int nP_1D;
	int halo_skip;
	int pk_skip;
	int use_one_pk;
	int n_pk_files;
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
	int r_bins;

	int use_n_min; // If =1 using n_min, otherwise using mass_min
	int n_min;
	double mass_min;
	int n_haloes_to_use; // Use the first N haloes

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
	int npts;
#ifdef GAS
	double OmegaB;
	double OmegaDM;
#endif
	double OmegaM;
	double OmegaL;
	double OmegaVDE;
	double H_0;
	double h;
	double delta_c;
	double s8;
	double n_s8;
	double Gn; 

	double err;
	double virial;
	double spin;

	double *a_hub;
	double *z_hub;
	double *Hubble;
	double *w;

} Cosmo;


extern struct mass_function
{
	int bins;
	int *n_tot;	
	int *n_bin;

	double Mmin;
	double Mmax;

	double *err;
	double *err_dn;
	double *ln_masses;
	double *sigmas;
	double *ln_sigmas;
	double *dlnSdlnM;
	double *num_masses;
	double *dn_num_masses; // The derivative is taken in the middle of the mass bin

	double *dn;
	double *n;
	double *th_masses;
	double *diff_tin;
	double *tin;

} MassFunc, ThMassFunc, *MassFuncZ;


extern struct num_density
{
	int npts;
	double *z;
	double zMax;
	double zMin;
	double *n_tin;
	double *n_num;

} NumDen; 


extern struct nfw
{
	int bins;
	double rs;
	double rho0;
	double *radius;
	double *overd;
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
	int virial;
	int spin;
	int conc;
	int all;

	double Xc;
	double Yc;
	double Zc;
	double VXc;
	double VYc;
	double VZc;
	double Lx;
	double Ly;
	double Lz;
	double Eax; 
	double Eay; 
	double Eaz; 
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
	double rho0;
	double c;
	double r2;
	double aa;
	double bb;
	double cc;
	double shape;
	double triax;
	double z_form;
	double c_nfw;
	double chi_nfw;
	double rs_nfw;
	double rho0_nfw;

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
	double Eax_gas;
	double Eay_gas;
	double Eaz_gas;
	double Ebx_gas;
	double Eby_gas;
	double Ebz_gas;
	double Ecx_gas;
	double Ecy_gas;
	double Ecz_gas;
	double lambda_gas;
	double lambdaE_gas;
	double b_gas;
	double c_gas;
	double Ekin_gas;
	double Epot_gas;
	double b_fraction;
	double Cum_u_gas;
	double T_gas;
#ifdef EXTRA_GAS
	double X_dm;
	double Y_dm;
	double Z_dm;
	double VX_dm;
	double VY_dm;
	double VZ_dm;
	double X_gas;
	double Y_gas;
	double Z_gas;
	double VX_gas;
	double VY_gas;
	double VZ_gas;
#endif
#endif  // GAS
} *Haloes, *SubHaloes, **pHaloes, **pSubHaloes;


extern struct full_catalogue
{
	int numFiles;
	int *numHaloes;
	double *z;
	char **urls;
	char **urls_profiles;
	char **urls_satellites;
	char **urls_particles;

} FullCat, *pFullCat;


extern struct internal_urls
{
	char *halo_file;
	char *profiles_file;
	char *pk_file;
	char *a_outputs;
	char *output_prefix;
	char *hubble_file;
	char *pk_list;
	char *halo_list;
	char *profile_list;
	char *subhalo_list;

} Urls, *pUrls;


extern struct halo_properties
{
	int ID;
	int z_bins;
	int r_bins;
	int n_bins;
	double avgSub;
	double z;

	double c_0;
	double c_02;
	double c_sig;
	double c_sig2;
	double *c;
	double *c_c0;
	double *p_c;
	double *err_p_c;

	double *costh;
	double *costh_count;
	double costh0;

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
	
	double *mass;
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

	double *cosphi;
	double *cosphi_count;
	double *cum_cosphi;
	double cosphi0;

	int *n_vel_sub;
	double *vel_sub;
	double *p_vel_sub;
	double vel_0;

	double *mass_sub;
	double avgMass;
	double *ecc;
	int *n_sub;
	int *cum_n_sub;
	int *n_ecc;

#ifdef GAS
	double *gas_T;
	double *gas_u;
	double *gas_fraction;
	double ln_M0;
	double T_alpha;
#endif
} SubHaloZ, HaloZ, *HaloProperties, *SubHaloProperties;
