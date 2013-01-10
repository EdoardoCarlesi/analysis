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
#include "../libio/halo_io.h"

#define LINE_SIZE 2048

void pread_halo_file()
{
	int n=0, j=0, thr=0, vir=0, conc=0, spin=0, skip=0;
	double a, minMass;
	char dummyline[LINE_SIZE]; 
	FILE *h_file=NULL;
	
	if(ThisTask==0)
		skip = 1;
	else
		skip=0;
	
	Settings.tick=0;
	h_file = fopen(pUrls[ThisTask].halo_file, "r");

	if(h_file==NULL) 
		fprintf(stderr, "Task=%d, halo file not found:%s\n", ThisTask, pUrls[ThisTask].halo_file);
	else
		fprintf(stderr, "Task=%d, found halo file:%s\n", ThisTask, pUrls[ThisTask].halo_file);

		pSettings[ThisTask].n_haloes = get_lines(h_file, pUrls[ThisTask].halo_file) - skip;

		pHaloes[ThisTask] = (struct halo*) calloc(pSettings[ThisTask].n_haloes, sizeof(struct halo));
		minMass = Settings.thMass;

		while(!feof(h_file))
		{
			fgets(dummyline, LINE_SIZE, h_file);
			if(j>=skip) 
			{
//fprintf(stdout, "Task=%d, j=%d, n=%d, skip=%d, n_haloes=%d, %s\n", ThisTask, j, n, skip, pSettings[ThisTask].n_haloes, dummyline);	
	sscanf(dummyline,
#ifndef GAS
	"%d  %d  %d  %lf %d  %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %d  %lf %lf %lf \
	 %lf %lf %lf \
	", 
#else // There is a gas component
#ifndef EXTRA_GAS
	"%d  %d  %d  %lf %d  %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %d  %lf %lf %lf \
	 %lf %lf %lf \
	 %d  %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	", 
#else
	"%d  %d  %d  %lf %d  %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %d  %lf %lf %lf \
	 %lf %lf %lf \
	 %d  %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf \
	", 
#endif
#endif // GAS
	&pHaloes[ThisTask][n].id,	&pHaloes[ThisTask][n].host,	&pHaloes[ThisTask][n].n_satellites, 	
	&pHaloes[ThisTask][n].Mvir, 	&pHaloes[ThisTask][n].n_part, 	&pHaloes[ThisTask][n].Xc, 	
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
	);

	pHaloes[ThisTask][n].aa = 1.0; // In the new catalogues haloes' major axis is normalized to one
#ifdef  GAS
#ifndef EXTRA_GAS
	// The cumulative u has to be read from the profile catalogues
	pHaloes[ThisTask][n].Cum_u_gas = 0.0;
#endif
#endif
//	fprintf(stdout, "\nTask=%d, setting additional halo properties.\n", ThisTask);
//	set_additional_halo_properties(n);
	//fprintf(stdout, "\nTask=%d, setting additional done.\n", ThisTask);

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
#endif 
#endif
#endif
		// Checking the various threshold conditions
	if(pHaloes[ThisTask][n].Mvir>minMass) 
	{
	thr++;

		if(sqrt(pHaloes[ThisTask][n].th_vir*pHaloes[ThisTask][n].th_vir) < Cosmo.virial) 
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
	}

		if(Settings.thNum < pHaloes[ThisTask][n].n_part)
		{
			Settings.haloes_over_thnum = pHaloes[ThisTask][n].id + 1;
	}

		//	print_counter(2000);
	
			n++;
			j++;
		} else { // Skip this line
		j++;
	}
} // Finished counting haloes over threshold conditions

		Settings.haloes_over_threshold=thr;
		Settings.virialized_haloes=vir;
		Settings.virialized_concentration=conc;
		Settings.spin_criterion=spin;

		determine_simulation_settings();	

	fclose(h_file);
}



void pset_additional_halo_properties(int n)
{
	double cc, vv, c;

	cc = pHaloes[ThisTask][n].cc;
	vv = sqrt(pHaloes[ThisTask][n].Mvir/pHaloes[ThisTask][n].Rvir*Settings.Gn);
	c = pHaloes[ThisTask][n].Rvir/pHaloes[ThisTask][n].r2;

	pHaloes[ThisTask][n].c = c;
	pHaloes[ThisTask][n].AngMom=pHaloes[ThisTask][n].lambda*sqrt(2)*pHaloes[ThisTask][n].Mvir*pHaloes[ThisTask][n].Rvir*vv;
	pHaloes[ThisTask][n].c_a = pHaloes[ThisTask][n].cc/pHaloes[ThisTask][n].aa;
	pHaloes[ThisTask][n].triax = (pow(pHaloes[ThisTask][n].aa, 2.0) - pow(pHaloes[ThisTask][n].bb, 2.0))/(pow(pHaloes[ThisTask][n].aa, 2.0) - cc);
	pHaloes[ThisTask][n].ecc = sqrt(1 - 2*pow(pHaloes[ThisTask][n].lambdaE,2));
	pHaloes[ThisTask][n].th_vir=2*pHaloes[ThisTask][n].Ekin/pHaloes[ThisTask][n].Epot;	
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



void pread_profiles_file()
{
	int nr=0, k=0, j=0, np1=0, np2=0, halo_bins=0, counter=0, kk=0;
	double a, over, over1, over2, over3, radius1, radius2, err_p;
	char dummyline[LINE_SIZE]; 
	FILE *h_file=NULL;

		h_file = fopen(Urls_internal.profiles_file, "r");
		if(h_file==NULL) 
			fprintf(stderr, "Profiles file not found:%s\n", Urls_internal.profiles_file);
		else
			fprintf(stderr, "Found profiles file:%s\n", Urls_internal.profiles_file);

		Settings.tick=0;
		Settings.Gn=1.;
		err_p=1; over3=500;

		if(Cosmo.err >0. && Cosmo.err <5.) 
			err_p = Cosmo.err;

		while(counter < Settings.haloes_over_threshold)
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



void pget_halo_files_urls()
{
	int n=0;
	char dummyline[100], *url_fc, *url_fc_pro, *url_fc_sub;
	FILE *fc=NULL, *fc_pro=NULL, *fc_sub=NULL;

	fprintf(stdout, "\nget_halo_files_urls() for %d haloes.\n", FC.numFiles);

	url_fc=pUrls[ThisTask].halo_list;
	url_fc_pro=pUrls[ThisTask].profile_list;
	url_fc_sub=pUrls[ThisTask].subhalo_list; 

	fc 	= fopen(url_fc,"r");
	fc_pro  = fopen(url_fc_pro,"r");
	fc_sub  = fopen(url_fc_sub,"r");
	
	if(fc==NULL) 
		fprintf(stderr, "%s %s\n","Could not find halo list file: ", url_fc);

	if(fc_pro==NULL) 
		fprintf(stderr, "%s %s\n","Could not find profiles list file: ", url_fc_pro);

	if(fc_sub==NULL) 
		fprintf(stderr, "%s %s\n","Could not find subhalo list file: ", url_fc_sub);

	if(pFC[ThisTask].numFiles==0) 
		pFC[ThisTask].numFiles = get_lines(fc, url_fc);

		pFC[ThisTask].urls = (char **) calloc(pFC[ThisTask].numFiles, sizeof(char *));
		pFC[ThisTask].urls_profiles = (char **) calloc(pFC[ThisTask].numFiles, sizeof(char *));
		pFC[ThisTask].urls_satellites = (char **) calloc(pFC[ThisTask].numFiles, sizeof(char *));

		for(n=0; n<pFC[ThisTask].numFiles; n++)
		{
			fgets(dummyline,100,fc);
			dummyline[strlen(dummyline)-1]='\0';
			pFC[ThisTask].urls[n] = (char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(pFC[ThisTask].urls[n],dummyline);
		}

	fclose(fc);
	
		for(n=0; n<pFC[ThisTask].numFiles; n++)
		{
			fgets(dummyline,100,fc_sub);
			dummyline[strlen(dummyline)-1]='\0';
			pFC[ThisTask].urls_satellites[n] = (char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(pFC[ThisTask].urls_satellites[n],dummyline);
		}
	fclose(fc_sub);

		for(n=0; n<pFC[ThisTask].numFiles; n++)
		{
			fgets(dummyline,100,fc_pro);
			dummyline[strlen(dummyline)-1]='\0';
			pFC[ThisTask].urls_profiles[n] = (char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(pFC[ThisTask].urls_profiles[n],dummyline);
		}
	fclose(fc_pro);
}



