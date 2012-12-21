#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "read_io.h"
#include "halo_io.h"
#include "../general_functions.h"
#include "../general_variables.h"
#include "../libcosmo/cosmological_relations.h"



void determine_simulation_settings()
{
#ifndef GAS
	Settings.pMass = haloes[0].Mvir/haloes[0].n_part;
	Settings.rho_0 = Settings.pMass*pow(Settings.nP_1D,3)*pow(Settings.box_size, -3);
	fprintf(stderr, "\nSetting particle mass: %e, rho_0: %e\n", Settings.pMass, Settings.rho_0);
#else
	Settings.dmMass  = haloes[0].M_dm/haloes[0].N_dm;
	Settings.gasMass = haloes[0].M_gas/haloes[0].N_gas; 
	Settings.rho_dm = Settings.dmMass*pow(Settings.nP_1D,3)*pow(Settings.box_size, -3);
	Settings.rho_b = Settings.gasMass*pow(Settings.nP_1D,3)*pow(Settings.box_size, -3);
	Settings.rho_0 = Settings.rho_b + Settings.rho_dm;
	fprintf(stderr, "\nSetting gasMass: %e, rho_b: %e\n", Settings.dmMass, Settings.rho_dm);
	fprintf(stderr, "\nSetting dmMass: %e, rho_dm: %e\n", Settings.gasMass, Settings.rho_b);
#endif

		Settings.rho_c = (3*haloes[0].Mvir)/(4*3.14*haloes[0].Rvir*haloes[0].Rvir*haloes[0].Rvir*200);

		fprintf(stderr, "\nThere are %d haloes over the %e mass threshold of which:\n\
		%d virialized\n\
		%d with the right concentration\n\
		%d satisfying the spin criterion\n\
		and %d haloes with more than %d particles.\n", 
		Settings.haloes_over_threshold,	Settings.thMass,
		Settings.virialized_haloes,
		Settings.virialized_concentration, 
		Settings.spin_criterion,
		Settings.haloes_over_thnum, Settings.thNum
		);

} 



	 /* Set additional halo properties derived from the basic ones */
void set_additional_halo_properties(int n)
{
	double cc, vv, c;

	cc = haloes[n].cc;
	vv = sqrt(haloes[n].Mvir/haloes[n].Rvir*Settings.Gn);
	c = haloes[n].Rvir/haloes[n].r2;

	haloes[n].c = c;
	haloes[n].AngMom=haloes[n].lambda*sqrt(2)*haloes[n].Mvir*haloes[n].Rvir*vv;
	haloes[n].c_a = haloes[n].cc/haloes[n].aa;
	haloes[n].triax = (pow(haloes[n].aa, 2.0) - pow(haloes[n].bb, 2.0))/(pow(haloes[n].aa, 2.0) - cc);
	haloes[n].ecc = sqrt(1 - 2*pow(haloes[n].lambdaE,2));
	haloes[n].th_vir=2*haloes[n].Ekin/haloes[n].Epot;	
	haloes[n].delta_c = (200/3)*c*c*c*(1./(log(1+c) - c/(1+c)));
#ifdef GAS
		haloes[n].N_dm = haloes[n].n_part - haloes[n].N_gas;
		haloes[n].M_dm = haloes[n].Mvir - haloes[n].M_gas;
		haloes[n].b_fraction = haloes[n].M_gas/haloes[n].Mvir;
#ifndef EXTRA_GAS
		haloes[n].T_gas = 0.0; 
#else
		haloes[n].T_gas = convert_u_to_T(haloes[n].Cum_u_gas);
		haloes[n].X_dm = (haloes[n].Mvir*haloes[n].Xc - haloes[n].X_gas)/haloes[n].M_dm;
		haloes[n].Y_dm = (haloes[n].Mvir*haloes[n].Yc - haloes[n].Y_gas)/haloes[n].M_dm;
		haloes[n].Z_dm = (haloes[n].Mvir*haloes[n].Zc - haloes[n].Z_gas)/haloes[n].M_dm;
		haloes[n].VX_dm = (haloes[n].Mvir*haloes[n].VXc - haloes[n].VX_gas)/haloes[n].M_dm;
		haloes[n].VY_dm = (haloes[n].Mvir*haloes[n].VYc - haloes[n].VY_gas)/haloes[n].M_dm;
		haloes[n].VZ_dm = (haloes[n].Mvir*haloes[n].VZc - haloes[n].VZ_gas)/haloes[n].M_dm;
#endif
#endif
}



void read_halo_file()
{
	int n=0, j=0, thr=0, vir=0, conc=0, spin=0;
	double a, minMass;
	char dummyline[10000]; 
	FILE *h_file=NULL;

	Settings.tick=0;
	h_file = fopen(Urls_internal.halo_file, "r");

	if(h_file==NULL) 
		fprintf(stderr, "%s %s \n", "Halo file not found: ", Urls_internal.halo_file);
	else
		fprintf(stderr, "%s %s \n", "Found halo file: ", Urls_internal.halo_file);

		Settings.n_haloes = get_lines(h_file, Urls_internal.halo_file) - Settings.halo_skip;

		haloes = (struct halo*) calloc(Settings.n_haloes, sizeof(struct halo));
		minMass = Settings.thMass;

		while(!feof(h_file))
		{
			fgets(dummyline, 10000, h_file);
				if(j>=Settings.halo_skip) 
			{	
#ifdef AHF_v1
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
	&haloes[n].id, 		&haloes[n].host, 	&haloes[n].n_satellites, 	&haloes[n].Mvir, 
	&haloes[n].n_part, 	&haloes[n].Xc, 		&haloes[n].Yc, 			&haloes[n].Zc, 
	&haloes[n].VXc, 	&haloes[n].VYc, 	&haloes[n].VZc, 		&haloes[n].Rvir, 
	&haloes[n].Rmax, 	&haloes[n].r2, 		&a, 				&a,
	&haloes[n].Vmax, 	&a, 			&a, 				&haloes[n].lambda, // First 20 entries
	&haloes[n].lambdaE, 	&haloes[n].Lx, 		&haloes[n].Ly, 			&haloes[n].Lz, 
	&haloes[n].bb,		&haloes[n].cc, 		&haloes[n].Eax, 		&haloes[n].Eay, 
	&haloes[n].Eaz,		&a, 			&a, 				&a, 
	&a, 			&a,			&a, 				&a,
	&haloes[n].n_bins,	&a, 			&haloes[n].Ekin, 		&haloes[n].Epot, // First 40 entries
	&a, 			&a, 			&haloes[n].c_nfw
#ifdef GAS
	, &haloes[n].N_gas, 	&haloes[n].M_gas, 	&haloes[n].lambda_gas, 		&haloes[n].lambdaE_gas, 
	&a, 			&a, 			&a, 				&haloes[n].b_gas, 
	&haloes[n].c_gas, 	&haloes[n].Eax_gas, 	&haloes[n].Eay_gas, 		&haloes[n].Eaz_gas, 
	&a, 			&a, 			&a,				&a, 
	&a, 			&a, 			&haloes[n].Ekin_gas, 		&haloes[n].Epot_gas
#ifdef EXTRA_GAS
	, &haloes[n].X_gas,	&haloes[n].Y_gas, 	&haloes[n].Z_gas, 		&haloes[n].VX_gas, 
	&haloes[n].VY_gas, 	&haloes[n].VZ_gas,	&haloes[n].Cum_u_gas
#endif
#endif
	);
	haloes[n].aa = 1.0; // In the new catalogues haloes' major axis is normalized to one
#ifdef  GAS
#ifndef EXTRA_GAS
	// The cumulative u has to be read from the profile catalogues
	haloes[n].Cum_u_gas = 0.0;
#endif
#endif

#else // If using old AHF catalogues
	sscanf(dummyline, 
	"%d  %d  %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %d %lf %lf %lf %lf %lf", 
	&haloes[n].n_part, &haloes[n].nv_part, 
	&haloes[n].Xc, &haloes[n].Yc, 
	&haloes[n].Zc, &haloes[n].VXc, 
	&haloes[n].VYc, &haloes[n].VZc, 
	&haloes[n].Mvir, &haloes[n].Rvir,
	&haloes[n].Vmax, &a, &a, &haloes[n].lambda, 
	&haloes[n].Lx, &haloes[n].Ly, &haloes[n].Lz, 
	&a, &haloes[n].Eax, &haloes[n].Eay, // First 20 entries
	&haloes[n].Eaz, &haloes[n].bb, &a, &a, &a, &haloes[n].cc, &a, &a, &a, &a,
	&a, &haloes[n].n_bins, &haloes[n].Ekin, &haloes[n].Epot, &a, &a, &a
	);	
#endif // Use AHF-v1.0

	set_additional_halo_properties(n);

#ifdef USE_UNIT_MPC
	haloes[n].Rvir *= 1.e-3;
	haloes[n].Xc *= 1.e-3;
	haloes[n].Yc *= 1.e-3;
	haloes[n].Zc *= 1.e-3;
	haloes[n].VXc *= 1.e-3;
	haloes[n].VYc *= 1.e-3;
	haloes[n].VZc *= 1.e-3;
#ifdef GAS
#ifdef EXTRA_GAS
	haloes[n].X_gas *= 1.e-3;
	haloes[n].Y_gas *= 1.e-3;
	haloes[n].Z_gas *= 1.e-3;
	haloes[n].VX_gas *= 1.e-3;
	haloes[n].VY_gas *= 1.e-3;
	haloes[n].VZ_gas *= 1.e-3;
	haloes[n].X_dm *= 1.e-3;
	haloes[n].Y_dm *= 1.e-3;
	haloes[n].Z_dm *= 1.e-3;
	haloes[n].VX_dm *= 1.e-3;
	haloes[n].VY_dm *= 1.e-3;
	haloes[n].VZ_dm *= 1.e-3;
#endif 
#endif
#endif
		// Checking the various threshold conditions
	if(haloes[n].Mvir>minMass) 
	{
	thr++;

		if(sqrt(haloes[n].th_vir*haloes[n].th_vir) < Cosmo.virial) 
		{

		if(haloes[n].c_nfw == -1) 
		{
			haloes[n].conc=0;		
		} else {
			haloes[n].conc=1;
			conc++;
		}

		if(haloes[n].lambda > Cosmo.spin)
		{
			haloes[n].spin=0;
		} else {
			haloes[n].spin=1;
			spin++;	
		}

			vir++;

		haloes[n].virial=1;
		} else {
			haloes[n].conc=0;
			haloes[n].virial=0;
	 	}
	}

		if(Settings.thNum < haloes[n].n_part)
		{
			Settings.haloes_over_thnum = haloes[n].id + 1;
	}

			print_counter(2000);
	
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



void read_profiles_file()
{
	int nr=0, k=0, j=0, np1=0, np2=0, halo_bins=0, counter=0, kk=0;
	double a, over, over1, over2, over3, radius1, radius2, err_p;
	char dummyline[4096]; 
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
			fgets(dummyline, 4096, h_file);
			if(j>= Settings.halo_skip) 
			{
				sscanf(dummyline, "%lf %d %lf %lf %lf", &radius1, &np1, &a, &over1, &over);	
				halo_bins = haloes[counter].n_bins;

				if(kk==0)
				{	
					haloes[counter].radius = (double *) calloc(halo_bins,sizeof(double));
					haloes[counter].rho = (double *) calloc(halo_bins,sizeof(double));
					haloes[counter].over_rho = (double *) calloc(halo_bins,sizeof(double));
					haloes[counter].bin = (double *) calloc(halo_bins,sizeof(double));
					haloes[counter].err = (double *) calloc(halo_bins,sizeof(double));
					haloes[counter].over_err = (double *) calloc(halo_bins,sizeof(double));
				}

				haloes[counter].radius[kk] = sqrt(radius1*radius1);	
				haloes[counter].rho[kk] = over1;
				haloes[counter].over_rho[kk] = Cosmo.OmegaM*over;
				haloes[counter].bin[kk] = np1-np2;
				haloes[counter].err[kk] = err_p*over1/sqrt(np1-np2); 
				haloes[counter].over_err[kk] = err_p*over/sqrt(np1-np2); 

			if(haloes[counter].err[kk] == 0.) 
				haloes[counter].err[kk] = 2*err_p*haloes[counter].err[kk-1]; 
	
			if(haloes[counter].over_err[kk] == 0.) 
				haloes[counter].over_err[kk] = 2*err_p*haloes[counter].over_err[kk-1]; 

			if(radius1<0.) 
			{
				nr++;
			}
			np2=np1;
			kk++;
	
		if(kk == halo_bins) 
		{
			haloes[counter].neg_r_bins=nr; 
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



void get_halo_files_urls()
{
	int n=0;
	char dummyline[100], *url_fc, *url_fc_pro, *url_fc_sub;
	FILE *fc=NULL, *fc_pro=NULL, *fc_sub=NULL;

	fprintf(stdout, "\nget_halo_files_urls() for %d haloes.\n", FC.numFiles);

	url_fc=Urls_internal.halo_list;
	url_fc_pro=Urls_internal.profile_list;
	url_fc_sub=Urls_internal.subhalo_list; 

	fc 	= fopen(url_fc,"r");
	fc_pro  = fopen(url_fc_pro,"r");
	fc_sub  = fopen(url_fc_sub,"r");
	
	if(fc==NULL) 
		fprintf(stderr, "%s %s\n","Could not find halo list file: ", url_fc);

	if(fc_pro==NULL) 
		fprintf(stderr, "%s %s\n","Could not find profiles list file: ", url_fc_pro);

	if(fc_sub==NULL) 
		fprintf(stderr, "%s %s\n","Could not find subhalo list file: ", url_fc_sub);

	if(FC.numFiles==0) 
		FC.numFiles = get_lines(fc, url_fc);

		FC.urls = (char **) calloc(FC.numFiles, sizeof(char *));
		FC.urls_profiles = (char **) calloc(FC.numFiles, sizeof(char *));
		FC.urls_satellites = (char **) calloc(FC.numFiles, sizeof(char *));

		for(n=0; n<FC.numFiles; n++)
		{
			fgets(dummyline,100,fc);
			dummyline[strlen(dummyline)-1]='\0';
			FC.urls[n] = (char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(FC.urls[n],dummyline);
		}

	fclose(fc);
	
		for(n=0; n<FC.numFiles; n++)
		{
			fgets(dummyline,100,fc_sub);
			dummyline[strlen(dummyline)-1]='\0';
			FC.urls_satellites[n] = (char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(FC.urls_satellites[n],dummyline);
		}
	fclose(fc_sub);

		for(n=0; n<FC.numFiles; n++)
		{
			fgets(dummyline,100,fc_pro);
			dummyline[strlen(dummyline)-1]='\0';
			FC.urls_profiles[n] = (char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(FC.urls_profiles[n],dummyline);
		}
	fclose(fc_pro);
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

		if(sqrt(haloes[bb].th_vir*haloes[bb].th_vir)<Cosmo.virial)
		{
			cc_ids[hh] = bb;
			hh++;
		}
	} while(hh<n_halos_comp);

	return cc_ids;
}
