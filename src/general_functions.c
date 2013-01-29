#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "general_variables.h"
#include "general_functions.h"

#ifdef WITH_MPI
#include <mpi.h>
#include "libparallel/general.h"
#endif



void initialize_internal_variables(char **argv){

	INFO_MSG("Initializing internal variables");
	
#ifdef PRINT_INFO
	int kk=0;
 	INFO_MSG("Printing info");
	for(kk=1; kk<40; kk++) 
		fprintf(stdout, "argv[%d]: %s \n", kk, argv[kk]);
#endif

	int count=1;

	// Basic working file urls
	Urls.a_outputs = argv[count++];
	Urls.halo_file = argv[count++];
	Urls.profiles_file = argv[count++];
	Urls.pk_file = argv[count++];

	// Basic simulation and analysis settings
	Settings.box_size = atof(argv[count++]); 
	Settings.n_part_1D = atoi(argv[count++]); 
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

//	fprintf(stdout, "Prefix=%s\n", Urls.output_prefix);

	set_halo_selection_criterion();
	
	default_init();
}



void default_init()
{
	// Setting extra useful variables
	Cosmo.H_0=Cosmo.h*100;
	Cosmo.G=6.672e-8;
	ThMassFunc.Mmin=Settings.Mmin;
	ThMassFunc.Mmax=Settings.Mmax;
	MassFunc.bins  = Settings.n_bins; 
	ThMassFunc.bins = Settings.n_bins_th;
	Settings.use_cat=Urls.nCatalogueFiles-Settings.cat_number;
	
	// Init some commonly used structures to default values
	GrowthFac.z = (double *) calloc(1, sizeof(double));
	GrowthFac.a = (double *) calloc(1, sizeof(double));
	GrowthFac.gf= (double *) calloc(1, sizeof(double));
	Pks = (struct power_spectrum *) calloc(1, sizeof(struct power_spectrum));
}



void print_counter(int freq){

	if(Settings.tick==freq) 
	{
		fprintf(stdout,"."); 
		Settings.tick=0;
	} else {
		Settings.tick++;
	}
}



void normalize_to_one(char * url_in, char * url_out)
{
		int inv=0, k=0, size=0;
		double norm=1.; 
		double *x=NULL, *y=NULL;
		char line[256];

		FILE *file_in = fopen(url_in, "r");
		FILE *file_out = fopen(url_out, "w");

		size = get_lines(file_in, url_in);

		x   = (double*) calloc(size,sizeof(double));
		y   = (double*) calloc(size,sizeof(double));

		for(k=0; k<size; k++)
		{
			fgets(line,256,file_in);
			sscanf(line,"%lf %lf", &x[k], &y[k]);
		}

		if (inv==1) 
		{
			norm = y[size-1];
			fprintf(stdout, "norm: %lf \n", y[size-1]);
		} else {
			norm = y[0];
			fprintf(stdout, "norm: %lf \n", y[0]);
		}

		for(k=0; k<size; k++)
			fprintf(file_out, "%e  %e\n", x[k], y[k]/norm);

	fclose(file_in);
	fclose(file_out);
}


	// Can be useful to check wether the criteria are changed somewhere else
void check_condition_consistency()
{
	int count=0;
		
		if(Settings.use_all == 1)
			Settings.use_all = 0;			
			else 
				count++;

			if(Settings.use_vir == 1)
					Settings.use_vir = 0;			
				else 
					count++;

				if(Settings.use_conc == 1)
					Settings.use_conc = 0;			
					else 
						count++;

			if(Settings.use_spin == 1)
				Settings.use_spin = 0;			
			else 
				count++;

		if(Settings.use_mass ==1)
			Settings.use_mass=0;			
			else 
				count++;

	if(count > 1)
		WARNING("More than one halo condition has been set", "check_condition_consistency()");
}



int halo_condition(int i)
{
	int condition=0;

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

	return nTot;
}



void set_halo_selection_criterion()
{
		// Init all criteria to zero
	Settings.use_mass = 0;
	Settings.use_vir = 0;
	Settings.use_spin = 0;
	Settings.use_conc = 0;
	Settings.use_all = 0;
	
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
