#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "read_io.h"
#include "halo_io.h"
#include "../general_functions.h"
#include "../general_variables.h"
#include "../libcosmo/cosmological_relations.h"

#ifdef WITH_MPI
#include "../libparallel/general.h"
#endif

#define LINE_SIZE 4096



void determine_simulation_settings()
{
	struct halo *HALO;
	struct general_settings *SETTINGS;

#ifdef WITH_MPI
	HALO = pHaloes[ThisTask];
	SETTINGS = &pSettings[ThisTask];
#else
	HALO = Haloes;
	SETTINGS = &Settings;
#endif

#ifdef WITH_MPI
	if(ThisTask==0)
#endif
	{
#ifndef GAS
	Settings.pMass = HALO[0].Mvir/HALO[0].n_part;
	Settings.rho_0 = Settings.pMass*pow(Settings.n_part_1D,3)*pow(Settings.box_size, -3);
	fprintf(stdout, "\nSetting particle mass: %e, rho_0: %e\n", Settings.pMass, Settings.rho_0);
#else
	Settings.dmMass  = HALO[0].M_dm/HALO[0].N_dm;
	Settings.gasMass = HALO[0].M_gas/HALO[0].N_gas; 
	Settings.rho_dm = Settings.dmMass*pow(Settings.nP_1D,3)*pow(Settings.box_size, -3);
	Settings.rho_b = Settings.gasMass*pow(Settings.nP_1D,3)*pow(Settings.box_size, -3);
	Settings.rho_0 = Settings.rho_b + Settings.rho_dm;
	fprintf(stdout, "\nSetting gasMass: %e, rho_b: %e\n", Settings.dmMass, Settings.rho_dm);
	fprintf(stdout, "\nSetting dmMass: %e, rho_dm: %e\n", Settings.gasMass, Settings.rho_b);
#endif

		Settings.rho_c = (3*HALO[0].Mvir) /
		(4*3.14*HALO[0].Rvir*HALO[0].Rvir*HALO[0].Rvir*200);
	}

		// Now set the mass_min to the value determined by the n_part
		// Check how to do this with DM and GAS particles
#ifndef GAS
		if(Settings.use_n_min == 1)
			Settings.mass_min = Settings.n_min * Settings.pMass; 
#else
		if(Settings.use_n_min == 1)
		{
			// Error message
		}
#endif
	
		fprintf(stdout, 
#ifdef WITH_MPI
		"\nTask=%d has %d haloes over the %e mass threshold of which:\n\
		- %d virialized\n\
		- %d with the right concentration\n\
		- %d satisfying the spin criterion\n\
		and %d haloes complying with all criteria.\n", 
		ThisTask, 
#else
		"\nThere are %d haloes over the %e mass threshold of which:\n\
		- %d virialized\n\
		- %d with the right concentration\n\
		- %d satisfying the spin criterion\n\
		and %d haloes complying with all criteria.\n", 
#endif
		SETTINGS->n_threshold, Settings.mass_min,
		SETTINGS->n_virialized,
		SETTINGS->n_concentration, 
		SETTINGS->n_spin,
		SETTINGS->n_all
		);

		fprintf(stdout, 
#ifdef WITH_MPI
			"\tTask=%d has box edges\n \tX=%lf  |  %lf\n \tY=%lf  |  %lf\n \tZ=%lf  |  %lf\n", 
				ThisTask,
#else
			"\tFound box edges\n \tX=%lf  |  %lf\n \tY=%lf  |  %lf\n \tZ=%lf  |  %lf\n",
#endif
				SETTINGS->box.X[0], SETTINGS->box.X[1], 
				SETTINGS->box.Y[0], SETTINGS->box.Y[1], 
				SETTINGS->box.Z[0], SETTINGS->box.Z[1] 
			);

}



void set_box(int i)
{
	struct halo *HALO;
	struct general_settings *SETTINGS;

	double x, y, z;

#ifdef WITH_MPI
	HALO = pHaloes[ThisTask];
	SETTINGS = &pSettings[ThisTask];
#else
	HALO = Haloes;
	SETTINGS = &Settings;
#endif

	x = HALO[i].Xc;
	y = HALO[i].Yc;
	z = HALO[i].Zc;

//		if(x != 0 && y != 0 && z != 0)
		{		// Check X
			if( x <	SETTINGS->box.X[0])
				SETTINGS->box.X[0] = x;

			if( x >	SETTINGS->box.X[1])
				SETTINGS->box.X[1] = x;

				// Check Y
			if( y <	SETTINGS->box.Y[0])
				SETTINGS->box.Y[0] = y;

			if( y >	SETTINGS->box.Y[1])
				SETTINGS->box.Y[1] = y;

				// Check Z
			if( z <	SETTINGS->box.Z[0])
				SETTINGS->box.Z[0] = z;

			if( z >	SETTINGS->box.Z[1])
				SETTINGS->box.Z[1] = z;
		}
}


	 /* Set additional halo properties derived from the basic ones */
void set_additional_halo_properties(int n)
{
	double c = 0.0; 
	struct halo *HALO;

#ifdef WITH_MPI
	HALO = pHaloes[ThisTask];
#else
	HALO = Haloes;
#endif

	c = HALO[n].Rvir/HALO[n].r2;

		// NFW initial guess parameters
	HALO[n].c = c;
	HALO[n].rho0 = (200/3)*c*c*c*(1./(log(1+c) - c/(1+c)));
	HALO[n].AngMom = HALO[n].lambda * sqrt(2)*HALO[n].Mvir *
			HALO[n].Rvir * sqrt(HALO[n].Mvir/HALO[n].Rvir);
	HALO[n].shape = HALO[n].cc/HALO[n].aa;
	HALO[n].triax = (pow(HALO[n].aa, 2.0) - pow(HALO[n].bb, 2.0))/
				(pow(HALO[n].aa, 2.0) - pow(HALO[n].cc, 2.0));
	HALO[n].ecc = sqrt(1 - 2*pow(HALO[n].lambdaE,2));
	HALO[n].th_vir=2*HALO[n].Ekin/HALO[n].Epot;	
	HALO[n].abs_th_vir=sqrt(HALO[n].th_vir*HALO[n].th_vir);	

#ifdef GAS
		HALO[n].N_dm = HALO[n].n_part - HALO[n].N_gas;
		HALO[n].M_dm = HALO[n].Mvir - HALO[n].M_gas;
		HALO[n].b_fraction = HALO[n].M_gas/HALO[n].Mvir;
#ifndef EXTRA_GAS
		HALO[n].T_gas = 0.0; 
#else
		HALO[n].T_gas = convert_u_to_T(HALO[n].Cum_u_gas);

		HALO[n].X_dm = 
	(HALO[n].Mvir*HALO[n].Xc - HALO[n].X_gas)/HALO[n].M_dm;

		HALO[n].Y_dm =
	(HALO[n].Mvir*HALO[n].Yc - HALO[n].Y_gas)/HALO[n].M_dm;

		HALO[n].Z_dm = 
	(HALO[n].Mvir*HALO[n].Zc - HALO[n].Z_gas)/HALO[n].M_dm;

		HALO[n].VX_dm = 
	(HALO[n].Mvir*HALO[n].VXc - HALO[n].VX_gas)/HALO[n].M_dm;

		HALO[n].VY_dm = 
	(HALO[n].Mvir*HALO[n].VYc - HALO[n].VY_gas)/HALO[n].M_dm;

		HALO[n].VZ_dm = 
	(HALO[n].Mvir*HALO[n].VZc - HALO[n].VZ_gas)/HALO[n].M_dm;
#endif // EXTRA_GAS
#endif // GAS
}



void get_halo_files_urls()
{
	int n=0, lin=0, lin_pro=0, lin_sub=0;
	char dummyline[200], url_fc[200], url_fc_pro[200], url_fc_sub[200];
	FILE *fc=NULL, *fc_pro=NULL, *fc_sub=NULL;

	struct internal_urls *URLS;

	fprintf(stdout, "\nget_halo_files_urls() for %d halo files.\n", Urls.nCatalogueFiles);

#ifdef WITH_MPI
	URLS = &pUrls[ThisTask];
#else
	URLS = &Urls;
#endif

	strcpy(url_fc, URLS->halo_list);
	strcpy(url_fc_pro, URLS->profile_list);
	strcpy(url_fc_sub, URLS->subhalo_list); 

	fc 	= fopen(url_fc,"r");
	fc_pro  = fopen(url_fc_pro,"r");
	fc_sub  = fopen(url_fc_sub,"r");
	
	if(fc==NULL) 
		ERROR("Could not find file:", url_fc);

	if(fc_pro==NULL) 
		ERROR("Could not find file:", url_fc_pro);

	if(fc_sub==NULL) 
		ERROR("Could not find file:", url_fc_sub);

		lin = get_lines(fc, url_fc);
		lin_pro = get_lines(fc_pro, url_fc_pro);
		lin_sub = get_lines(fc_sub, url_fc_sub);

		URLS->nCatalogueFiles = lin;

		URLS->urls = 
			(char **) calloc(URLS->nCatalogueFiles, sizeof(char *));
		URLS->urls_profiles = 
			(char **) calloc(URLS->nCatalogueFiles, sizeof(char *));
		URLS->urls_satellites = 
			(char **) calloc(URLS->nCatalogueFiles, sizeof(char *));

	if(lin > 0)
	{
		for(n=0; n<URLS->nCatalogueFiles; n++)
		{
			fgets(dummyline,200,fc);
			URLS->urls[n] = 
				(char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(URLS->urls[n], dummyline);
			URLS->urls[n][strlen(dummyline)-1]='\0';
		}
	fclose(fc);

		} else {
#ifdef WITH_MPI
			if(ThisTask==0)
#endif				
				WARNING("Trying to read empty file", url_fc);
		}

	if(lin_sub > 0)
	{
		for(n=0; n<URLS->nCatalogueFiles; n++)
		{
			fgets(dummyline,200,fc_sub);
			URLS->urls_satellites[n] = 
				(char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(URLS->urls_satellites[n], dummyline);
			URLS->urls_satellites[n][strlen(dummyline)-1]='\0';
		}
	fclose(fc_sub);

		} else {
#ifdef WITH_MPI
			if(ThisTask==0)
#endif
				WARNING("Trying to read empty file", url_fc_sub);
		}

	if(lin_pro > 0)
	{
		for(n=0; n<URLS->nCatalogueFiles; n++)
		{
			fgets(dummyline,200,fc_pro);
			URLS->urls_profiles[n] = 
				(char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(URLS->urls_profiles[n], dummyline);
			URLS->urls_profiles[n][strlen(dummyline)-1]='\0';
		}
	fclose(fc_pro);

		} else {
#ifdef WITH_MPI
			if(ThisTask==0)
#endif
				WARNING("Trying to read empty file", url_fc_pro);
		}

}



void set_halo_url()
{
	int n=0;
	struct internal_urls *URLS;
	
	n = Settings.use_cat;

#ifdef WITH_MPI
	URLS = &pUrls[ThisTask];
#else
	URLS = &Urls;
#endif

	URLS->halo_file = (char*) calloc(strlen(URLS->urls[n])-1, sizeof(char));
	strcpy(URLS->halo_file, URLS->urls[n]);
}



void read_halo_file()
{
	int n=0, j=0, thr=0, vir=0, conc=0, spin=0, skip=0, all=0, condition=0;
	double a=0; // Dummy variable to read useless columns
	char dummyline[LINE_SIZE]; 
	FILE *h_file;

	struct halo *HALO;
	struct general_settings *SETTINGS;
	struct internal_urls *URLS;
	
	Settings.tick=0;

#ifdef WITH_MPI
	if(ThisTask==0)
		skip = 1;
		else
			skip = 0;

	URLS = &pUrls[ThisTask];
	SETTINGS = &pSettings[ThisTask];
#else
	skip = 1;
	URLS = &Urls;
	SETTINGS = &Settings;
#endif

		h_file = fopen(URLS->halo_file, "r");
		SETTINGS->n_haloes = get_lines(h_file, URLS->halo_file) - skip;

#ifdef WITH_MPI
		pHaloes[ThisTask] = (struct halo*) calloc(SETTINGS->n_haloes, 
				sizeof(struct halo));
		HALO = pHaloes[ThisTask];
#else
		Haloes = (struct halo*) calloc(SETTINGS->n_haloes, sizeof(struct halo));
		HALO = Haloes;
#endif

	if(h_file==NULL) 
	{	
		ERROR("File not found", URLS->halo_file);

		} else {

#ifdef WITH_MPI
			fprintf(stdout, "Task=%d is reading from halo file:%s\n", 
				ThisTask, URLS->halo_file);
#else
			fprintf(stdout, "Reading from halo file:%s\n", URLS->halo_file);
#endif
		}


#ifdef WITH_MPI
	fprintf(stdout, "\nTask=%d is allocating memory for %d haloes...\n",
			ThisTask, SETTINGS->n_haloes);
#else
	fprintf(stdout, "\nAllocating memory for %d haloes...\n", SETTINGS->n_haloes);
#endif

	while(!feof(h_file) && n < SETTINGS->n_haloes)
	{
		fgets(dummyline, LINE_SIZE, h_file);
			
			if(j>=skip) 
			{
				sscanf(dummyline,
#ifndef GAS
	"%ld %d  %d  %lf %d  %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %d  %lf %lf %lf \
	 %lf %lf %lf \
	",
#else // There is a gas component
#ifndef EXTRA_GAS
	"%ld %d  %d  %lf %d  %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %d  %lf %lf %lf \
	 %lf %lf %lf \
	 %d  %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	", 
#else	// Use extra gas columns
	"%ld %d  %d  %lf %d  %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %d  %lf %lf %lf \
	 %lf %lf %lf \
	 %d  %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf \
	", 
#endif	// Extra Gas
#endif // GAS
	&HALO[n].id,		&HALO[n].host,		&HALO[n].n_satellites,	
	&HALO[n].Mvir, 		&HALO[n].n_part,	&HALO[n].Xc, 	
	&HALO[n].Yc, 		&HALO[n].Zc, 		&HALO[n].VXc, 	
	&HALO[n].VYc, 		&HALO[n].VZc, 		&HALO[n].Rvir, 
	&HALO[n].Rmax, 		&HALO[n].r2, 		&a, 				
	&a,			&HALO[n].Vmax, 		&a, 			
	&a, 			&HALO[n].lambda, 	// First 20 entries
	&HALO[n].lambdaE,	&HALO[n].Lx, 		&HALO[n].Ly, 			
	&HALO[n].Lz, 		&HALO[n].bb,		&HALO[n].cc, 		
	&HALO[n].Eax, 		&HALO[n].Eay, 		&HALO[n].Eaz,		
	&a, 			&a, 			&a, 
	&a, 			&a,			&a, 				
	&a,			&HALO[n].n_bins,	&a, 			
	&HALO[n].Ekin,	 	&HALO[n].Epot, 		// First 40 entries
	&a, 			&a,			&HALO[n].c_nfw
#ifdef GAS
	, &HALO[n].N_gas, 	&HALO[n].M_gas, 	&HALO[n].lambda_gas, 		
	&HALO[n].lambdaE_gas, 	&a,			&a, 			
	&a, 			&HALO[n].b_gas, 
	&HALO[n].c_gas, 	&HALO[n].Eax_gas, 	&HALO[n].Eay_gas, 		
	&HALO[n].Eaz_gas, 	&a,			&a, 			
	&a,			&a,			&a, 			
	&a, 			&HALO[n].Ekin_gas, 	&HALO[n].Epot_gas
#ifdef EXTRA_GAS
	, &HALO[n].X_gas,	&HALO[n].Y_gas, 	&HALO[n].Z_gas, 		
	&HALO[n].VX_gas, 	&HALO[n].VY_gas, 	&HALO[n].VZ_gas,	
	&HALO[n].Cum_u_gas
#endif
#endif
	); // */
#ifdef GAS
#ifndef EXTRA_GAS
	// The cumulative u has to be read from the profile catalogues
	HALO[n].Cum_u_gas = 0.0;
#endif
#endif
	// In the new catalogues haloes' major axis is normalized to one
	HALO[n].aa = 1.0; 

	set_additional_halo_properties(n);

#ifdef USE_UNIT_MPC
	HALO[n].Rvir *= 1.e-3;
	HALO[n].Xc *= 1.e-3;
	HALO[n].Yc *= 1.e-3;
	HALO[n].Zc *= 1.e-3;
	HALO[n].VXc *= 1.e-3;
	HALO[n].VYc *= 1.e-3;
	HALO[n].VZc *= 1.e-3;
#ifdef GAS
#ifdef EXTRA_GAS
	HALO[n].X_gas *= 1.e-3;
	HALO[n].Y_gas *= 1.e-3;
	HALO[n].Z_gas *= 1.e-3;
	HALO[n].VX_gas *= 1.e-3;
	HALO[n].VY_gas *= 1.e-3;
	HALO[n].VZ_gas *= 1.e-3;
	HALO[n].X_dm *= 1.e-3;
	HALO[n].Y_dm *= 1.e-3;
	HALO[n].Z_dm *= 1.e-3;
	HALO[n].VX_dm *= 1.e-3;
	HALO[n].VY_dm *= 1.e-3;
	HALO[n].VZ_dm *= 1.e-3;
#endif // Extra Gas
#endif // Gas
#endif // Use Mpc
		
		// Check the size of the box corresponding to the current task
	set_box(n);

		// Checking the various threshold conditions
	if(Settings.use_n_min == 1)
		{
			if(HALO[n].n_part>Settings.n_min) 
				condition = 1;
					else 
				condition = 0;
		} else {
			if(HALO[n].Mvir>Settings.mass_min) 
				condition = 1;
					else 
				condition = 0;
		}

	if(condition==1)
	{
		HALO[n].mass=1;
		thr++;

		if(HALO[n].abs_th_vir < Cosmo.virial) 
		{
			if(HALO[n].c_nfw == -1) 
			{
				HALO[n].conc=0;		
			} else {
				HALO[n].conc=1;
				conc++;
			}

				if(HALO[n].lambda > Cosmo.spin)
				{
					HALO[n].spin=0;		
				} else {
					HALO[n].spin=1;
					spin++;	
			}

			vir++;
		HALO[n].vir=1;
		} else {
			HALO[n].conc=0;
			HALO[n].vir=0;
	 	}
			if(	HALO[n].conc==1 &&
				HALO[n].vir==1 &&
				HALO[n].spin==1 ) 
	
				{		
					HALO[n].all = 1; 
						all++;
					} else {
				HALO[n].all = 0; 
				}
	}	// Condition == 1 */
			n++;
			j++;
		} else { // Skip this line

		HALO[n].mass = 0;
		j++;
	}
} // Finished counting haloes over threshold conditions

		SETTINGS->n_threshold=thr;
		SETTINGS->n_virialized=vir;
		SETTINGS->n_concentration=conc;
		SETTINGS->n_spin=spin;
		SETTINGS->n_all=all;
	
	determine_simulation_settings();	

		
	fclose(h_file);
}



void read_profiles_file()
{
	int i=0, j=0, k=0; 
	int npart=0, halo_bins=0, counter=0;
	int neg_r=0, npart_old=0;

	double radius, overd, dens, v_circ, a;
	double overd_old, radius_old;

#ifdef GAS
	double m_gas, u_gas;
#endif

	char dummyline[LINE_SIZE]; 
	FILE *p_file=NULL;

	struct halo *HALO;
	struct general_settings *SETTINGS;
	struct internal_urls *URLS;

#ifdef WITH_MPI
	URLS = &pUrls[ThisTask];
	SETTINGS = &pSettings[ThisTask];
#else
	URLS = &Urls;
	SETTINGS = &Settings;
#endif

		p_file = fopen(URLS->profiles_file, "r");
		SETTINGS->n_haloes = get_lines(p_file, URLS->profiles_file);

		if(p_file==NULL)
			ERROR("File not found", URLS->profiles_file);
		else
			fprintf(stdout, "Found profiles file:%s\n", 
				URLS->profiles_file);

#ifdef WITH_MPI
		HALO = pHaloes[ThisTask];
#else
		HALO = Haloes;
#endif

		while(counter < SETTINGS->n_threshold)
		{
			fgets(dummyline, LINE_SIZE, p_file);

			if(j>= SETTINGS->halo_skip) 
			{
				sscanf(dummyline, 
#ifndef GAS
		"%lf  %d   %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf \
		%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf \
		%lf  %lf  %lf  %lf", 
		&radius, &npart, &a, &overd, &dens, &v_circ, &a, &a, &a, &a, 
		&a,      &a, 	 &a, &a,     &a,    &a,      &a, &a, &a, &a, 
		&a, 	 &a,     &a, &a	// 24
#else
		"%lf  %d   %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf \
		%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf \
		%lf  %lf  %lf  %lf  %lf  %lf  %lf", 
		&radius, &npart, &a, &overd, &dens, &v_circ, &a, &a, &a, &a, 
		&a,      &a, 	 &a, &a,     &a,    &a,      &a, &a, &a, &a, 
		&a, 	 &a,     &a, &a,     &m_gas,&a	     &u_gas  // 27
#endif
	);
				halo_bins = HALO[counter].n_bins;

				if(i==0)
				{	
					HALO[counter].radius = 
						(double *) calloc(halo_bins,sizeof(double));
					HALO[counter].rho = 
						(double *) calloc(halo_bins,sizeof(double));
					HALO[counter].over_rho = 
						(double *) calloc(halo_bins,sizeof(double));
					HALO[counter].bin = 
						(double *) calloc(halo_bins,sizeof(double));
					HALO[counter].err = 
						(double *) calloc(halo_bins,sizeof(double));
					HALO[counter].over_err = 
						(double *) calloc(halo_bins,sizeof(double));
#ifdef GAS
					HALO[counter].u_gas = 
						(double *) calloc(halo_bins,sizeof(double));
					HALO[counter].m_gas = 
						(double *) calloc(halo_bins,sizeof(double));
#endif
				}

					HALO[counter].radius[i] = sqrt(radius*radius);	
					HALO[counter].rho[i] = dens;
					HALO[counter].over_rho[i] = Cosmo.OmegaM*overd;
					HALO[counter].bin[i] = npart-npart_old;
					HALO[counter].err[i] = overd/sqrt(npart-npart_old); 
					HALO[counter].over_err[i] = overd/sqrt(npart-npart_old); 

			if(HALO[counter].err[i] == 0.) 
				HALO[counter].err[i] = HALO[counter].err[i-1]; 
	
			if(HALO[counter].over_err[i] == 0.) 
				HALO[counter].over_err[i] = HALO[counter].over_err[i-1]; 

			if(radius<0.) 
			{
				neg_r++;
			}

			npart_old=npart;
			i++;
	
		if(i == halo_bins) 
		{
			HALO[counter].neg_r_bins=neg_r; 
	
			i=0;
			npart=0;
			neg_r=0;
			counter++;
		}

	overd_old=overd;
	radius_old=radius;
	npart_old=npart;

	j++;

	} 
		else 
	{ // Skip this line
		j++;
	}

	k++;

	}

	fclose(p_file);
}



int get_files_redshift_position(char *fileName)
{ 
	return (int) (atoi(&fileName[0])*10 + atoi(&fileName[1]))/10;
}



int* read_cross_correlated_haloes(char* halo_cc_list, int n_halos_comp, int* cc_ids){
	int hh=0, aa=0, bb=0; 
	char dummyline[128]; 
	FILE *cc_list = fopen(halo_cc_list, "r");
	
	fprintf(stdout, "Reading cross correlation file:%s\n", halo_cc_list);

	do {
		fgets(dummyline, 128, cc_list);
		sscanf(dummyline,"%d %d", &aa, &bb); 

		if(Haloes[bb].abs_th_vir<Cosmo.virial)
		{
			cc_ids[hh] = bb;
			hh++;
		}
	} while(hh<n_halos_comp);

	return cc_ids;
}
