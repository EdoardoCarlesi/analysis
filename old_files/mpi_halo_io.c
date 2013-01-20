#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>
#include "halo_io.h"
#include "general.h"
#include "../general_variables.h"
#include "../general_functions.h"
#include "../libio/read_io.h"

#define LINE_SIZE 2048


void mpi_determine_simulation_settings()
{
	struct halo *HALO;
	struct general_settings *SETTINGS;

#ifdef WITH_MPI
	HALO = pHaloes[ThisTask];
	SETTINGS = &pSettings[ThisTask];
#else
	HALO = haloes;
	SETTINGS = &Settings;
#endif

	if(ThisTask==0)
	{
#ifndef GAS
	Settings.pMass = HALO[0].Mvir/HALO[0].n_part;
	Settings.rho_0 = Settings.pMass*pow(Settings.nP_1D,3)*pow(Settings.box_size, -3);
	fprintf(stderr, "\nSetting particle mass: %e, rho_0: %e\n", Settings.pMass, Settings.rho_0);
#else
	Settings.dmMass  = HALO[0].M_dm/HALO[0].N_dm;
	Settings.gasMass = HALO[0].M_gas/HALO[0].N_gas; 
	Settings.rho_dm = Settings.dmMass*pow(Settings.nP_1D,3)*pow(Settings.box_size, -3);
	Settings.rho_b = Settings.gasMass*pow(Settings.nP_1D,3)*pow(Settings.box_size, -3);
	Settings.rho_0 = Settings.rho_b + Settings.rho_dm;
	fprintf(stderr, "\nSetting gasMass: %e, rho_b: %e\n", Settings.dmMass, Settings.rho_dm);
	fprintf(stderr, "\nSetting dmMass: %e, rho_dm: %e\n", Settings.gasMass, Settings.rho_b);
#endif

		Settings.rho_c = 
			(3*HALO[0].Mvir) /
		(4*3.14*HALO[0].Rvir*HALO[0].Rvir*HALO[0].Rvir*200);
	}

	//MPI_Bcast(&Settings, sizeof(struct general_settings), MPI_BYTE, 0, MPI_COMM_WORLD);

	if(Settings.use_n_min == 1)
	{
		fprintf(stderr, "\nTask=%d has %d haloes over the %d particles threshold of which:\n\
		- %d virialized\n\
		- %d with the right concentration\n\
		- %d satisfying the spin criterion\n\
		and %d haloes complying with all criteria.\n", 
		ThisTask, 
		SETTINGS->n_threshold, Settings.n_min,
		SETTINGS->n_virialized,
		SETTINGS->n_concentration, 
		SETTINGS->n_spin,
		SETTINGS->n_all
		);		

	} else {
		fprintf(stderr, "\nTask=%d has %d haloes over the %e mass threshold of which:\n\
		- %d virialized\n\
		- %d with the right concentration\n\
		- %d satisfying the spin criterion\n\
		and %d haloes complying with all criteria.\n", 
		ThisTask, 
		SETTINGS->n_threshold, Settings.mass_min,
		SETTINGS->n_virialized,
		SETTINGS->n_concentration, 
		SETTINGS->n_spin,
		SETTINGS->n_all
		);
	}
}



void mpi_read_halo_file()
{
	int b=0, n=0, j=0, thr=0, vir=0, conc=0, spin=0, skip=0, all=0, condition=0;
	double a=0; // Dummy variable to read columns
	char dummyline[4096]; 
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
	URLS = &Urls_internal;
	SETTINGS = &Settings;
#endif

		h_file = fopen(URLS->halo_file, "r");
		SETTINGS->n_haloes = get_lines(h_file, URLS->halo_file) - skip;

#ifdef WITH_MPI
		pHaloes[ThisTask] = (struct halo*) calloc(SETTINGS->n_haloes, 
				sizeof(struct halo));
		HALO = pHaloes[ThisTask];
#else
		HALO = haloes;
#endif


	if(h_file==NULL) 
	{
		fprintf(stderr, "Task=%d, halo file not found:%s\n", 
			ThisTask, URLS->halo_file);

		} else {
			fprintf(stderr, "Task=%d is reading from halo file:%s\n", 
				ThisTask, URLS->halo_file);
		}


	fprintf(stdout, "\nTask=%d is allocating memory for %d haloes...\n",
			ThisTask, SETTINGS->n_haloes);

	while(!feof(h_file) && n < SETTINGS->n_haloes)
	{
		fgets(dummyline, 4096, h_file);
			
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
	"%ld  %d  %d  %lf %d  %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %d  %lf %lf %lf \
	 %lf %lf %lf \
	 %d  %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	", 
#else	// Use extra gas columns
	"%ld  %d  %d  %lf %d  %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %d  %lf %lf %lf \
	 %lf %lf %lf \
	 %d  %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf \
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
	&HALO[n].Ekin,	 	&HALO[n].Epot, 	// First 40 entries
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

	mpi_set_additional_halo_properties(n);

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
		HALO[n].virial=1;
		} else {
			HALO[n].conc=0;
			HALO[n].virial=0;
	 	}
			if(	HALO[n].conc==1 &&
				HALO[n].virial==1 &&
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
		j++;
	}
} // Finished counting haloes over threshold conditions

		SETTINGS->n_threshold=thr;
		SETTINGS->n_virialized=vir;
		SETTINGS->n_concentration=conc;
		SETTINGS->n_spin=spin;
		SETTINGS->n_all=all;
	
	mpi_determine_simulation_settings();	

fclose(h_file);
}



void mpi_set_additional_halo_properties(int n)
{
	struct halo *HALO;

#ifdef WITH_MPI
	HALO = pHaloes[ThisTask];
#else
	HALO = haloes;
#endif

	double c = HALO[n].Rvir/HALO[n].r2;

	HALO[n].c = c;
	HALO[n].AngMom = HALO[n].lambda * sqrt(2)*HALO[n].Mvir *
			HALO[n].Rvir * sqrt(HALO[n].Mvir/HALO[n].Rvir);
	HALO[n].shape = HALO[n].cc/HALO[n].aa;
	HALO[n].triax = (pow(HALO[n].aa, 2.0) - pow(HALO[n].bb, 2.0))/
				(pow(HALO[n].aa, 2.0) - pow(HALO[n].cc, 2.0));
	HALO[n].ecc = sqrt(1 - 2*pow(HALO[n].lambdaE,2));
	HALO[n].th_vir=2*HALO[n].Ekin/HALO[n].Epot;	
	HALO[n].abs_th_vir=sqrt(HALO[n].th_vir*HALO[n].th_vir);	
	HALO[n].delta_c = (200/3)*c*c*c*(1./(log(1+c) - c/(1+c)));
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
#endif
#endif
}



void mpi_get_halo_files_urls()
{
	int n=0, lin=0, lin_pro=0, lin_sub=0;
	char dummyline[200], url_fc[200], url_fc_pro[200], url_fc_sub[200];
	FILE *fc=NULL, *fc_pro=NULL, *fc_sub=NULL;

	fprintf(stdout, "\nget_halo_files_urls() for %d halo files.\n", FC.numFiles);

	strcpy(url_fc, pUrls[ThisTask].halo_list);
	strcpy(url_fc_pro, pUrls[ThisTask].profile_list);
	strcpy(url_fc_sub, pUrls[ThisTask].subhalo_list); 

	fc 	= fopen(url_fc,"r");
	fc_pro  = fopen(url_fc_pro,"r");
	fc_sub  = fopen(url_fc_sub,"r");
	
	if(fc==NULL) 
		fprintf(stderr, "%s %s\n","Could not find halo list file: ", url_fc);

	if(fc_pro==NULL) 
		fprintf(stderr, "%s %s\n","Could not find profiles list file: ", url_fc_pro);

	if(fc_sub==NULL) 
		fprintf(stderr, "%s %s\n","Could not find subhalo list file: ", url_fc_sub);

		lin = get_lines(fc, url_fc);
		lin_pro = get_lines(fc_pro, url_fc_pro);
		lin_sub = get_lines(fc_sub, url_fc_sub);

		pFC[ThisTask].numFiles = lin;

		pFC[ThisTask].urls = (char **) calloc(pFC[ThisTask].numFiles, sizeof(char *));
		pFC[ThisTask].urls_profiles = (char **) calloc(pFC[ThisTask].numFiles, sizeof(char *));
		pFC[ThisTask].urls_satellites = (char **) calloc(pFC[ThisTask].numFiles, sizeof(char *));

	if(lin > 0)
	{
		for(n=0; n<pFC[ThisTask].numFiles; n++)
		{
			fgets(dummyline,200,fc);
			pFC[ThisTask].urls[n] = (char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(pFC[ThisTask].urls[n], dummyline);
			pFC[ThisTask].urls[n][strlen(dummyline)-1]='\0';
		}
	fclose(fc);

		} else {
			if(ThisTask==0)
				fprintf(stderr, "\nFile %s contains no URLS.\n", url_fc);
		}

	if(lin_sub > 0)
	{
		for(n=0; n<pFC[ThisTask].numFiles; n++)
		{
			fgets(dummyline,200,fc_sub);
			pFC[ThisTask].urls_satellites[n] = (char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(pFC[ThisTask].urls_satellites[n], dummyline);
			pFC[ThisTask].urls_satellites[n][strlen(dummyline)-1]='\0';
		}
	fclose(fc_sub);

		} else {
			if(ThisTask==0)
				fprintf(stderr, "\nFile %s contains no URLS.\n", url_fc_sub);
		}

	if(lin_pro > 0)
	{
		for(n=0; n<pFC[ThisTask].numFiles; n++)
		{
			fgets(dummyline,200,fc_pro);
			pFC[ThisTask].urls_profiles[n] = (char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(pFC[ThisTask].urls_profiles[n], dummyline);
			pFC[ThisTask].urls_profiles[n][strlen(dummyline)-1]='\0';
		}
	fclose(fc_pro);

		} else {
			if(ThisTask==0)
				fprintf(stderr, "\nFile %s contains no URLS.\n", url_fc_pro);
		}

}



void mpi_use_halo_url(int n)
{
	pUrls[ThisTask].halo_file = (char*) calloc(strlen(pFC[ThisTask].urls[n])-1, sizeof(char));
	strcpy(pUrls[ThisTask].halo_file, pFC[ThisTask].urls[n]);
}
