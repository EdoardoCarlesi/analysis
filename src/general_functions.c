#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "libio/io.h"
#include "libhalo/halo.h"
#include "libcosmo/cosmo.h"

#include "general_def.h"

#ifdef WITH_MPI
#include <mpi.h>
#include "libparallel/general.h"
#endif


/*
 * Define functions
 */
int subhalo_condition(int);


/*
 * Initialize functions
 */
void initialize_internal_variables(char **argv){

	INFO_MSG("Initializing internal variables");
	
#ifdef PRINT_INFO
	int kk=0;
 	INFO_MSG("Printing info");
	for(kk=1; kk<39; kk++) 
		fprintf(stdout, "argv[%d]: %s \n", kk, argv[kk]);
#endif

	int count=1;

		// Basic working file urls
		Urls.a_outputs = argv[count++];
		Urls.halo_file = argv[count++];
		Urls.profiles_file = argv[count++];
		Urls.pk_file = argv[count++];
		Urls.c_web_file = argv[count++];
		Urls.c_web_gas_file = argv[count++];

		// Basic simulation and analysis settings
		Settings.box_size = atof(argv[count++]); 
		Settings.n_part_1D = atoi(argv[count++]); 
		Settings.c_web_size = atoi(argv[count++]); 
		Settings.n_bins = atoi(argv[count++]);
		Settings.n_bins_th = atoi(argv[count++]);
		Settings.r_bins = atoi(argv[count++]);
		Settings.pk_skip = atoi(argv[count++]);
		Settings.halo_skip = atoi(argv[count++]);
		Settings.cat_number = atoi(argv[count++]);
		Settings.fit = atoi(argv[count++]);
		Settings.zStart = atof(argv[count++]);
		Settings.mass_min = atof(argv[count++]);
		Settings.Mmin = atof(argv[count++]);
		Settings.Mmax = atof(argv[count++]);
		Settings.Rmin = atof(argv[count++]);
		Settings.Rmax = atof(argv[count++]);
		Settings.l_web = atof(argv[count++]);
		Settings.n_min = atof(argv[count++]);
		Settings.use_n_min = atof(argv[count++]);
		Settings.n_haloes_to_use = atoi(argv[count++]); 
		Settings.use_criterion = atoi(argv[count++]); 
	
		// Cosmological parameters
		Cosmo.h = atof(argv[count++]);
		Cosmo.sigma8 = atof(argv[count++]);
		Cosmo.OmegaM = atof(argv[count++]);
		Cosmo.OmegaL = atof(argv[count++]);
		Cosmo.delta_c = atof(argv[count++]);
		Cosmo.spin = atof(argv[count++]);
		Cosmo.virial = atof(argv[count++]);
		GrowthFac.scale_k = atof(argv[count++]);
		NumDen.zMax = atof(argv[count++]);

		// Other urls related parameters
		Urls.output_prefix = argv[count++];
		Urls.halo_list = argv[count++];
		Urls.profile_list = argv[count++];
		Urls.subhalo_list = argv[count++];
		Urls.pk_list = argv[count++];
		Urls.nCatalogueFiles = atoi(argv[count++]);

	set_halo_selection_criterion();
	
	default_init();
}



void default_init()
{
	// Setting extra useful variables
	
	MF_INDEX = 0;
	PK_INDEX = 0;

	Settings.use_cat = Urls.nCatalogueFiles-Settings.cat_number;
	Settings.use_web = 0;

	MassFunc = malloc(30*sizeof(struct mass_function));
	SubMassFunc = malloc(30*sizeof(struct mass_function));
	ThMassFunc = malloc(30*sizeof(struct mass_function));
	HaloProperties = malloc(30*sizeof(struct halo_properties));
	SubHaloProperties = malloc(30*sizeof(struct halo_properties));

	Cosmo.H_0=Cosmo.h*100;
	Cosmo.G=6.672e-8;
	ThMassFunc[MF_INDEX].Mmin=Settings.Mmin;
	ThMassFunc[MF_INDEX].Mmax=Settings.Mmax;
	ThMassFunc[MF_INDEX].bins=Settings.n_bins_th;
	MassFunc[MF_INDEX].bins=Settings.n_bins; 

	// Init box - set min 75, max 0 as initial values that will be 
	// surely overwritten when the read_halo routine is called
	Settings.box.X[0][0] = 100;
	Settings.box.X[0][1] = 0;
	Settings.box.X[1][0] = 100;
	Settings.box.X[1][1] = 0;
	Settings.box.X[2][0] = 100;
	Settings.box.X[2][1] = 0;
	
	// Init some commonly used structures to default values
	GrowthFac.z = (double *) calloc(1, sizeof(double));
	GrowthFac.a = (double *) calloc(1, sizeof(double));
	GrowthFac.gf= (double *) calloc(1, sizeof(double));
	Pks = (struct power_spectrum *) calloc(1, sizeof(struct power_spectrum));
}


	// Can be useful to check wether the criteria are changed somewhere else
void check_condition_consistency()
{
	int count=0;
		
		if(Settings.use_all == 1)
			count++;
			else 
				Settings.use_all = 0;			

			if(Settings.use_vir == 1)
				count++;
				else 
					Settings.use_vir = 0;			

				if(Settings.use_conc == 1)
					count++;
					else 
					Settings.use_conc = 0;			

			if(Settings.use_spin == 1)
				count++;
				else 
			Settings.use_spin = 0;			

		if(Settings.use_mass == 1)
			count++;
			else 
		Settings.use_mass=0;			

	if(count > 1)
		WARNING("More than one halo condition has been set", "check_condition_consistency()");
}



int subhalo_condition(int i)
{

	if(Settings.use_sub == 1 && Haloes[i].host > 0)
	{	
		//fprintf(stderr, "i=%d", i);
		//fprintf(stderr, " UseSub=%d", Settings.use_sub);
		//fprintf(stderr, " Host=%llu\n", Haloes[i].host);
		return 1;
	} else {
		return 0;
	}
}



int halo_condition(int i)
{
	int j=0, condition=0;
	
	if(subhalo_condition(i) || Settings.use_sub == 0) 	
	{	

	if(Settings.use_spin == 1)
	{
		if(Haloes[i].spin == 1)
				condition = 1;
		else 
				condition = 0;
	}
		
		else if(Settings.use_conc == 1)
		{
			if(Haloes[i].conc == 1)
					condition = 1;
			else 
					condition = 0;
		}

			else if(Settings.use_mass == 1)
			{
				if(Haloes[i].mass == 1)
						condition = 1;
				else 
						condition = 0;
			}

				else if(Settings.use_vir == 1)
				{
					if(Haloes[i].vir == 1)
							condition = 1;
					else 
							condition = 0;
				}

			else if(Settings.use_all == 1)
			{
				if(Haloes[i].all == 1)
						condition = 1;
				else 
						condition = 0;
			}
	
		else 
			condition = 0;
	}

	if(condition == 1 && Settings.use_web == 1)
	{
		if(Haloes[i].web_type[Settings.use_web_type] == 1)
			condition = 1;
	}

	return condition;
}



int n_haloes_per_criterion()
{
	int nTot=0;

		if(Settings.use_spin == 1)
			nTot = Settings.n_spin;

			else if(Settings.use_all == 1)
				nTot = Settings.n_all;

				else if(Settings.use_vir == 1)
					nTot = Settings.n_virialized;

					else if(Settings.use_conc == 1)
						nTot = Settings.n_concentration;

				else if(Settings.use_mass == 1)
					nTot = Settings.n_threshold;

			else
				nTot = Settings.n_haloes;
		// When using subhaloes neglect the former choices
		if(Settings.use_sub == 1)
			nTot = Settings.n_sub_threshold;

	return nTot;
}



void set_halo_selection_criterion()
{
		// Init all criteria to zero
	Settings.use_mass = 0;
	Settings.use_spin = 0;
	Settings.use_conc = 0;
	Settings.use_vir = 0;
	Settings.use_all = 0;
	Settings.use_sub = 0;
	
	switch(Settings.use_criterion)
	{	
		case 1:
			Settings.use_mass = 1;
			INFO_MSG("Using halo mass/particle number selection criterion");
				break;	
		case 2:
			Settings.use_vir = 1;
			INFO_MSG("Using halo virialization selection criterion");
				break;	
	
		case 3:
			Settings.use_spin = 1;
			INFO_MSG("Using halo spin selection criterion");
				break;	
	
		case 4:
			Settings.use_conc = 1;
			INFO_MSG("Using halo concentration selection criterion");
				break;	
	
		case 5:
			Settings.use_all = 1;
			INFO_MSG("Using all halo combined selection criteria");
				break;	

		default:
			WARNING("No halo selection criterion specified",
				"set_halo_selection_criterion()");
	}
}
