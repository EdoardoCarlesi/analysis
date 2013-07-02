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
	if(ThisTask==0)
	{
#ifndef GAS
	Settings.pMass = pHaloes[ThisTask][0].Mvir/pHaloes[ThisTask][0].n_part;
	Settings.rho_0 = Settings.pMass*pow(Settings.nP_1D,3)*pow(Settings.box_size, -3);
	fprintf(stderr, "\nSetting particle mass: %e, rho_0: %e\n", Settings.pMass, Settings.rho_0);
#else
	Settings.dmMass  = pHaloes[ThisTask][0].M_dm/pHaloes[ThisTask][0].N_dm;
	Settings.gasMass = pHaloes[ThisTask][0].M_gas/pHaloes[ThisTask][0].N_gas; 
	Settings.rho_dm = Settings.dmMass*pow(Settings.nP_1D,3)*pow(Settings.box_size, -3);
	Settings.rho_b = Settings.gasMass*pow(Settings.nP_1D,3)*pow(Settings.box_size, -3);
	Settings.rho_0 = Settings.rho_b + Settings.rho_dm;
	fprintf(stderr, "\nSetting gasMass: %e, rho_b: %e\n", Settings.dmMass, Settings.rho_dm);
	fprintf(stderr, "\nSetting dmMass: %e, rho_dm: %e\n", Settings.gasMass, Settings.rho_b);
#endif

		Settings.rho_c = 
			(3*pHaloes[ThisTask][0].Mvir) /
		(4*3.14*pHaloes[ThisTask][0].Rvir*pHaloes[ThisTask][0].Rvir*pHaloes[ThisTask][0].Rvir*200);
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
		pSettings[ThisTask].n_threshold, Settings.n_min,
		pSettings[ThisTask].n_virialized,
		pSettings[ThisTask].n_concentration, 
		pSettings[ThisTask].n_spin,
		pSettings[ThisTask].n_all
		);		

	} else {
		fprintf(stderr, "\nTask=%d has %d haloes over the %e mass threshold of which:\n\
		- %d virialized\n\
		- %d with the right concentration\n\
		- %d satisfying the spin criterion\n\
		and %d haloes complying with all criteria.\n", 
		ThisTask, 
		pSettings[ThisTask].n_threshold, Settings.mass_min,
		pSettings[ThisTask].n_virialized,
		pSettings[ThisTask].n_concentration, 
		pSettings[ThisTask].n_spin,
		pSettings[ThisTask].n_all
		);
	}
}



void mpi_read_halo_file()
{
	int b=0, n=0, j=0, thr=0, vir=0, conc=0, spin=0, skip=0, all=0, condition=0;
	double a=0; // Dummy variable to read columns
	char dummyline[4096]; 
	FILE *h_file;

	if(ThisTask==0)
		skip = 1;
	else
		skip = 0;
	
	Settings.tick=0;
	h_file = fopen(pUrls[ThisTask].halo_file, "r");

	if(h_file==NULL) 
	{
		fprintf(stderr, "Task=%d, halo file not found:%s\n", ThisTask, pUrls[ThisTask].halo_file);

		} else {
			fprintf(stderr, "Task=%d is reading from halo file:%s\n", ThisTask, pUrls[ThisTask].halo_file);
		}
		
		pSettings[ThisTask].n_haloes = get_lines(h_file, pUrls[ThisTask].halo_file) - skip;

		pHaloes[ThisTask] = (struct halo*) calloc(pSettings[ThisTask].n_haloes, sizeof(struct halo));

	fprintf(stdout, "\nTask=%d is allocating memory for %d haloes...\n",
			ThisTask, pSettings[ThisTask].n_haloes);

	while(!feof(h_file) && n < pSettings[ThisTask].n_haloes)
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
	&pHaloes[ThisTask][n].id,	&pHaloes[ThisTask][n].host,	&pHaloes[ThisTask][n].n_satellites,	
	&pHaloes[ThisTask][n].Mvir, 	&pHaloes[ThisTask][n].n_part,	&pHaloes[ThisTask][n].Xc, 	
	&pHaloes[ThisTask][n].Yc, 	&pHaloes[ThisTask][n].Zc, 	&pHaloes[ThisTask][n].VXc, 	
	&pHaloes[ThisTask][n].VYc, 	&pHaloes[ThisTask][n].VZc, 	&pHaloes[ThisTask][n].Rvir, 
	&pHaloes[ThisTask][n].Rmax, 	&pHaloes[ThisTask][n].r2, 	&a, 				
	&a,				&pHaloes[ThisTask][n].Vmax, 	&a, 			
	&a, 				&pHaloes[ThisTask][n].lambda, 	// First 20 entries
	&pHaloes[ThisTask][n].lambdaE, 	&pHaloes[ThisTask][n].Lx, 	&pHaloes[ThisTask][n].Ly, 			
	&pHaloes[ThisTask][n].Lz, 	&pHaloes[ThisTask][n].bb,	&pHaloes[ThisTask][n].cc, 		
	&pHaloes[ThisTask][n].Eax, 	&pHaloes[ThisTask][n].Eay, 	&pHaloes[ThisTask][n].Eaz,		
	&a, 				&a, 				&a, 
	&a, 				&a,				&a, 				
	&a,				&pHaloes[ThisTask][n].n_bins,	&a, 			
	&pHaloes[ThisTask][n].Ekin, 	&pHaloes[ThisTask][n].Epot, 	// First 40 entries
	&a, 				&a, 				&pHaloes[ThisTask][n].c_nfw
#ifdef GAS
	, &pHaloes[ThisTask][n].N_gas, 	&pHaloes[ThisTask][n].M_gas, 	&pHaloes[ThisTask][n].lambda_gas, 		
	&pHaloes[ThisTask][n].lambdaE_gas, &a, 				&a, 			
	&a, 				&pHaloes[ThisTask][n].b_gas, 
	&pHaloes[ThisTask][n].c_gas, 	&pHaloes[ThisTask][n].Eax_gas, 	&pHaloes[ThisTask][n].Eay_gas, 		
	&pHaloes[ThisTask][n].Eaz_gas, 	&a, 				&a, 			
	&a,				&a, 				&a, 			
	&a, 				&pHaloes[ThisTask][n].Ekin_gas, &pHaloes[ThisTask][n].Epot_gas
#ifdef EXTRA_GAS
	, &pHaloes[ThisTask][n].X_gas,	&pHaloes[ThisTask][n].Y_gas, 	&pHaloes[ThisTask][n].Z_gas, 		
	&pHaloes[ThisTask][n].VX_gas, 	&pHaloes[ThisTask][n].VY_gas, 	&pHaloes[ThisTask][n].VZ_gas,	
	&pHaloes[ThisTask][n].Cum_u_gas
#endif
#endif
	); // */
#ifdef GAS
#ifndef EXTRA_GAS
	// The cumulative u has to be read from the profile catalogues
	pHaloes[ThisTask][n].Cum_u_gas = 0.0;
#endif
#endif
	// In the new catalogues haloes' major axis is normalized to one
	pHaloes[ThisTask][n].aa = 1.0; 

	mpi_set_additional_halo_properties(n);

#ifdef USE_UNIT_MPC
	pHaloes[ThisTask][n].Rvir *= 1.e-3;
	pHaloes[ThisTask][n].Xc *= 1.e-3;
	pHaloes[ThisTask][n].Yc *= 1.e-3;
	pHaloes[ThisTask][n].Zc *= 1.e-3;
	pHaloes[ThisTask][n].VXc *= 1.e-3;
	pHaloes[ThisTask][n].VYc *= 1.e-3;
	pHaloes[ThisTask][n].VZc *= 1.e-3;
#ifdef GAS
#ifdef EXTRA_GAS
	pHaloes[ThisTask][n].X_gas *= 1.e-3;
	pHaloes[ThisTask][n].Y_gas *= 1.e-3;
	pHaloes[ThisTask][n].Z_gas *= 1.e-3;
	pHaloes[ThisTask][n].VX_gas *= 1.e-3;
	pHaloes[ThisTask][n].VY_gas *= 1.e-3;
	pHaloes[ThisTask][n].VZ_gas *= 1.e-3;
	pHaloes[ThisTask][n].X_dm *= 1.e-3;
	pHaloes[ThisTask][n].Y_dm *= 1.e-3;
	pHaloes[ThisTask][n].Z_dm *= 1.e-3;
	pHaloes[ThisTask][n].VX_dm *= 1.e-3;
	pHaloes[ThisTask][n].VY_dm *= 1.e-3;
	pHaloes[ThisTask][n].VZ_dm *= 1.e-3;
#endif // Extra Gas
#endif // Gas
#endif // Use Mpc

		// Checking the various threshold conditions
	if(Settings.use_n_min == 1)
		{
			if(pHaloes[ThisTask][n].n_part>Settings.n_min) 
				condition = 1;
					else 
				condition = 0;
		} else {
			if(pHaloes[ThisTask][n].Mvir>Settings.mass_min) 
				condition = 1;
					else 
				condition = 0;
		}

	if(condition==1)
	{
		thr++;

		if(pHaloes[ThisTask][n].abs_th_vir < Cosmo.virial) 
		{
			if(pHaloes[ThisTask][n].c_nfw == -1) 
			{
				pHaloes[ThisTask][n].conc=0;		
			} else {
				pHaloes[ThisTask][n].conc=1;
				conc++;
			}

				if(pHaloes[ThisTask][n].lambda > Cosmo.spin)
				{
					pHaloes[ThisTask][n].spin=0;		
				} else {
					pHaloes[ThisTask][n].spin=1;
					spin++;	
			}

			vir++;
		pHaloes[ThisTask][n].virial=1;
		} else {
			pHaloes[ThisTask][n].conc=0;
			pHaloes[ThisTask][n].virial=0;
	 	}
			if(	pHaloes[ThisTask][n].conc==1 &&
				pHaloes[ThisTask][n].virial==1 &&
				pHaloes[ThisTask][n].spin==1 ) 
	
				{		
					pHaloes[ThisTask][n].all = 1; 
						all++;
					} else {
				pHaloes[ThisTask][n].all = 0; 
				}
	}	// Condition == 1 */
			n++;
			j++;
		} else { // Skip this line
		j++;
	}
} // Finished counting haloes over threshold conditions

		pSettings[ThisTask].n_threshold=thr;
		pSettings[ThisTask].n_virialized=vir;
		pSettings[ThisTask].n_concentration=conc;
		pSettings[ThisTask].n_spin=spin;
		pSettings[ThisTask].n_all=all;
	
		mpi_determine_simulation_settings();	

fclose(h_file);
}



void mpi_set_additional_halo_properties(int n)
{
	double c = pHaloes[ThisTask][n].Rvir/pHaloes[ThisTask][n].r2;

	pHaloes[ThisTask][n].c = c;
	pHaloes[ThisTask][n].AngMom = pHaloes[ThisTask][n].lambda * sqrt(2)*pHaloes[ThisTask][n].Mvir *
			pHaloes[ThisTask][n].Rvir * sqrt(pHaloes[ThisTask][n].Mvir/pHaloes[ThisTask][n].Rvir);
	pHaloes[ThisTask][n].shape = pHaloes[ThisTask][n].cc/pHaloes[ThisTask][n].aa;
	pHaloes[ThisTask][n].triax = (pow(pHaloes[ThisTask][n].aa, 2.0) - pow(pHaloes[ThisTask][n].bb, 2.0))/
				(pow(pHaloes[ThisTask][n].aa, 2.0) - pow(pHaloes[ThisTask][n].cc, 2.0));
	pHaloes[ThisTask][n].ecc = sqrt(1 - 2*pow(pHaloes[ThisTask][n].lambdaE,2));
	pHaloes[ThisTask][n].th_vir=2*pHaloes[ThisTask][n].Ekin/pHaloes[ThisTask][n].Epot;	
	pHaloes[ThisTask][n].abs_th_vir=sqrt(pHaloes[ThisTask][n].th_vir*pHaloes[ThisTask][n].th_vir);	
	pHaloes[ThisTask][n].delta_c = (200/3)*c*c*c*(1./(log(1+c) - c/(1+c)));
#ifdef GAS
		pHaloes[ThisTask][n].N_dm = pHaloes[ThisTask][n].n_part - pHaloes[ThisTask][n].N_gas;
		pHaloes[ThisTask][n].M_dm = pHaloes[ThisTask][n].Mvir - pHaloes[ThisTask][n].M_gas;
		pHaloes[ThisTask][n].b_fraction = pHaloes[ThisTask][n].M_gas/pHaloes[ThisTask][n].Mvir;
#ifndef EXTRA_GAS
		pHaloes[ThisTask][n].T_gas = 0.0; 
#else
		pHaloes[ThisTask][n].T_gas = convert_u_to_T(pHaloes[ThisTask][n].Cum_u_gas);

		pHaloes[ThisTask][n].X_dm = 
	(pHaloes[ThisTask][n].Mvir*pHaloes[ThisTask][n].Xc - pHaloes[ThisTask][n].X_gas)/pHaloes[ThisTask][n].M_dm;

		pHaloes[ThisTask][n].Y_dm =
	(pHaloes[ThisTask][n].Mvir*pHaloes[ThisTask][n].Yc - pHaloes[ThisTask][n].Y_gas)/pHaloes[ThisTask][n].M_dm;

		pHaloes[ThisTask][n].Z_dm = 
	(pHaloes[ThisTask][n].Mvir*pHaloes[ThisTask][n].Zc - pHaloes[ThisTask][n].Z_gas)/pHaloes[ThisTask][n].M_dm;

		pHaloes[ThisTask][n].VX_dm = 
	(pHaloes[ThisTask][n].Mvir*pHaloes[ThisTask][n].VXc - pHaloes[ThisTask][n].VX_gas)/pHaloes[ThisTask][n].M_dm;

		pHaloes[ThisTask][n].VY_dm = 
	(pHaloes[ThisTask][n].Mvir*pHaloes[ThisTask][n].VYc - pHaloes[ThisTask][n].VY_gas)/pHaloes[ThisTask][n].M_dm;

		pHaloes[ThisTask][n].VZ_dm = 
	(pHaloes[ThisTask][n].Mvir*pHaloes[ThisTask][n].VZc - pHaloes[ThisTask][n].VZ_gas)/pHaloes[ThisTask][n].M_dm;
#endif
#endif
}



void mpi_read_profiles_file()
{
	int nr=0, k=0, j=0, np1=0, np2=0, halo_bins=0, counter=0, kk=0;
	double a, over, over1, over2, over3, radius1, radius2, err_p;
	char dummyline[LINE_SIZE]; 
	FILE *h_file;

		h_file = fopen(Urls_internal.profiles_file, "r");
		if(h_file==NULL) 
			fprintf(stderr, "Profiles file not found:%s\n", Urls_internal.profiles_file);
		else
			fprintf(stderr, "Found profiles file:%s\n", Urls_internal.profiles_file);

		Settings.tick=0;
		err_p=1; over3=500;

		if(Cosmo.err >0. && Cosmo.err <5.) 
			err_p = Cosmo.err;

		while(counter < Settings.n_threshold)
		{
			fgets(dummyline, LINE_SIZE, h_file);
			if(j>= Settings.halo_skip) 
			{
				sscanf(dummyline, "%lf %d %lf %lf %lf", &radius1, &np1, &a, &over1, &over);	
				halo_bins = haloes[counter].n_bins;

				if(kk==0)
				{	
					pHaloes[ThisTask][counter].radius = (double *) calloc(halo_bins,sizeof(double));
					pHaloes[ThisTask][counter].rho = (double *) calloc(halo_bins,sizeof(double));
					pHaloes[ThisTask][counter].over_rho = (double *) calloc(halo_bins,sizeof(double));
					pHaloes[ThisTask][counter].bin = (double *) calloc(halo_bins,sizeof(double));
					pHaloes[ThisTask][counter].err = (double *) calloc(halo_bins,sizeof(double));
					pHaloes[ThisTask][counter].over_err = (double *) calloc(halo_bins,sizeof(double));
				}

				pHaloes[ThisTask][counter].radius[kk] = sqrt(radius1*radius1);	
				pHaloes[ThisTask][counter].rho[kk] = over1;
				pHaloes[ThisTask][counter].over_rho[kk] = Cosmo.OmegaM*over;
				pHaloes[ThisTask][counter].bin[kk] = np1-np2;
				pHaloes[ThisTask][counter].err[kk] = err_p*over1/sqrt(np1-np2); 
				pHaloes[ThisTask][counter].over_err[kk] = err_p*over/sqrt(np1-np2); 

			if(pHaloes[ThisTask][counter].err[kk] == 0.) 
				pHaloes[ThisTask][counter].err[kk] = 2*err_p*pHaloes[ThisTask][counter].err[kk-1]; 
	
			if(pHaloes[ThisTask][counter].over_err[kk] == 0.) 
				pHaloes[ThisTask][counter].over_err[kk] = 2*err_p*pHaloes[ThisTask][counter].over_err[kk-1]; 

			if(radius1<0.) 
			{
				nr++;
			}
			np2=np1;
			kk++;
	
		if(kk == halo_bins) 
		{
			pHaloes[ThisTask][counter].neg_r_bins=nr; 
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

	fclose(h_file);
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
