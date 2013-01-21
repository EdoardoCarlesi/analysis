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
	int ThisTask = 0;
	HALO = Haloes;
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

		// NFW approximate parameters
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

	struct full_catalogue *FULLCATALOGUE;
	struct internal_urls *URLS;

	fprintf(stdout, "\nget_halo_files_urls() for %d halo files.\n", FullCat.numFiles);

#ifdef WITH_MPI
	FULLCATALOGUE =	&pFullCat[ThisTask];
	URLS = &pUrls[ThisTask];
#else
	FULLCATALOGUE = &FullCat;
	URLS = &Urls;
#endif

	strcpy(url_fc, URLS->halo_list);
	strcpy(url_fc_pro, URLS->profile_list);
	strcpy(url_fc_sub, URLS->subhalo_list); 

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

		FULLCATALOGUE->numFiles = lin;

		FULLCATALOGUE->urls = 
			(char **) calloc(FULLCATALOGUE->numFiles, sizeof(char *));
		FULLCATALOGUE->urls_profiles = 
			(char **) calloc(FULLCATALOGUE->numFiles, sizeof(char *));
		FULLCATALOGUE->urls_satellites = 
			(char **) calloc(FULLCATALOGUE->numFiles, sizeof(char *));

	if(lin > 0)
	{
		for(n=0; n<FULLCATALOGUE->numFiles; n++)
		{
			fgets(dummyline,200,fc);
			FULLCATALOGUE->urls[n] = 
				(char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(FULLCATALOGUE->urls[n], dummyline);
			FULLCATALOGUE->urls[n][strlen(dummyline)-1]='\0';
		}
	fclose(fc);

		} else {
			if(ThisTask==0)
				fprintf(stderr, "\nFile %s contains no URLS.\n", url_fc);
		}

	if(lin_sub > 0)
	{
		for(n=0; n<FULLCATALOGUE->numFiles; n++)
		{
			fgets(dummyline,200,fc_sub);
			FULLCATALOGUE->urls_satellites[n] = 
				(char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(FULLCATALOGUE->urls_satellites[n], dummyline);
			FULLCATALOGUE->urls_satellites[n][strlen(dummyline)-1]='\0';
		}
	fclose(fc_sub);

		} else {
			if(ThisTask==0)
				fprintf(stderr, "\nFile %s contains no URLS.\n", url_fc_sub);
		}

	if(lin_pro > 0)
	{
		for(n=0; n<FULLCATALOGUE->numFiles; n++)
		{
			fgets(dummyline,200,fc_pro);
			FULLCATALOGUE->urls_profiles[n] = 
				(char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(FULLCATALOGUE->urls_profiles[n], dummyline);
			FULLCATALOGUE->urls_profiles[n][strlen(dummyline)-1]='\0';
		}
	fclose(fc_pro);

		} else {
			if(ThisTask==0)
				fprintf(stderr, "\nFile %s contains no URLS.\n", url_fc_pro);
		}

}



void use_halo_url(int n)
{
	struct internal_urls *URLS;
	struct full_catalogue *FULLCATALOGUE;

#ifdef WITH_MPI
	FULLCATALOGUE =	&pFullCat[ThisTask];
	URLS = &pUrls[ThisTask];
#else
	FULLCATALOGUE = &FullCat;
	URLS = &Urls;
#endif

	URLS->halo_file = (char*) calloc(strlen(FULLCATALOGUE->urls[n])-1, sizeof(char));
	strcpy(URLS->halo_file, FULLCATALOGUE->urls[n]);
}



void read_halo_file()
{
	int b=0, n=0, j=0, thr=0, vir=0, conc=0, spin=0, skip=0, all=0, condition=0;
	double a=0; // Dummy variable to read columns
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
	int ThisTask = 0;
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
	
	determine_simulation_settings();	

	fclose(h_file);
}



void read_profiles_file()
{
	int nr=0, k=0, j=0, np1=0, np2=0, halo_bins=0, counter=0, kk=0;
	double a, over, over1, over2, over3, radius1, radius2, err_p;
	char dummyline[LINE_SIZE]; 
	FILE *p_file=NULL;

	struct halo *HALO;
	struct general_settings *SETTINGS;
	struct internal_urls *URLS;
	
	Settings.tick=0;

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
			fprintf(stderr, "Profiles file not found:%s\n", 
				URLS->profiles_file);
		else
			fprintf(stderr, "Found profiles file:%s\n", 
				URLS->profiles_file);

#ifdef WITH_MPI
		HALO = pHaloes[ThisTask];
#else
		HALO = Haloes;
#endif

		Settings.tick=0;
		err_p=1; over3=500;

		if(Cosmo.err >0. && Cosmo.err <5.) 
			err_p = Cosmo.err;

		while(counter < SETTINGS->n_threshold)
		{
			fgets(dummyline, LINE_SIZE, p_file);
			if(j>= SETTINGS->halo_skip) 
			{
				sscanf(dummyline, "%lf %d %lf %lf %lf", 
					&radius1, &np1, &a, &over1, &over);	

				halo_bins = HALO[counter].n_bins;

				if(kk==0)
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
				}

				HALO[counter].radius[kk] = sqrt(radius1*radius1);	
				HALO[counter].rho[kk] = over1;
				HALO[counter].over_rho[kk] = Cosmo.OmegaM*over;
				HALO[counter].bin[kk] = np1-np2;
				HALO[counter].err[kk] = err_p*over1/sqrt(np1-np2); 
				HALO[counter].over_err[kk] = err_p*over/sqrt(np1-np2); 

			if(HALO[counter].err[kk] == 0.) 
				HALO[counter].err[kk] = 2*err_p*HALO[counter].err[kk-1]; 
	
			if(HALO[counter].over_err[kk] == 0.) 
				HALO[counter].over_err[kk] = 2*err_p*HALO[counter].over_err[kk-1]; 

			if(radius1<0.) 
			{
				nr++;
			}

			np2=np1;
			kk++;
	
		if(kk == halo_bins) 
		{
			HALO[counter].neg_r_bins=nr; 
			kk=0;
			np1=0;
			nr=0;
			counter++;
		}

	over2=over1;
	radius2=radius1;
	np2=np1;
	print_counter(10000);

	j++;
	} else { // Skip this line
	j++;
	}
	k++;
	}

	fprintf(stdout,"\n");

	close(p_file);
}



int get_files_redshift_position(char *fileName)
{ 
	return (int) (atoi(&fileName[0])*10 + atoi(&fileName[1]))/10;
}



int* read_cross_correlated_haloes(char* halo_cc_list, int n_halos_comp, int* cc_ids){
	int hh=0, aa=0, bb=0; 
	char dummyline[128]; 
	FILE *cc_list = fopen(halo_cc_list, "r");
	
	fprintf(stderr, "Reading cross correlation file:%s\n", halo_cc_list);

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
