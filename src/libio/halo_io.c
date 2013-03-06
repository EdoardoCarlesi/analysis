#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../libcosmo/cosmo.h"
#include "../libmath/math.h"
#include "../general_def.h"

#include "io.h"

#ifdef WITH_MPI
#include "../libparallel/general.h"
#endif

#define LINE_SIZE 4096

/*
 *  Function definition
 */
void set_additional_halo_properties(int);

void determine_simulation_settings(void);

void set_box(int);

/*
 *  Function implementation
 */
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

	if(ThisTask==0)
	{
#ifndef GAS
		Settings.pMass = HALO[0].Mvir/HALO[0].n_part;
		Settings.rho_0 = Settings.pMass*pow3(Settings.n_part_1D)*(1./pow3(Settings.box_size));
		fprintf(stdout, "\nSetting particle mass: %e, rho_0: %e\n", Settings.pMass, Settings.rho_0);
#else
		Settings.dmMass  = HALO[1].dm.M/(double)HALO[1].dm.N;
		Settings.gasMass = HALO[1].gas.M/(double)HALO[1].gas.N; 

		Settings.rho_dm = Settings.dmMass*pow3(Settings.n_part_1D)*(1./pow3(Settings.box_size));
		Settings.rho_b = Settings.gasMass*pow3(Settings.n_part_1D)*(1./pow3(Settings.box_size));
		Settings.rho_0 = Settings.rho_b + Settings.rho_dm;
		fprintf(stdout, "\nSetting gasMass: %e, rho_b: %e\n", Settings.dmMass, Settings.rho_dm);
		fprintf(stdout, "\nSetting dmMass: %e, rho_dm: %e\n", Settings.gasMass, Settings.rho_b);
#endif
		Settings.rho_c = (3*HALO[0].Mvir)/(4*3.14*pow3(HALO[0].Rvir)*200);
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
	
		fprintf(stdout, "\nTask=%d has %d haloes over the %e mass threshold of which:\n\
					- %d virialized\n\
					- %d with the right concentration\n\
					- %d satisfying the spin criterion\n\
					and %d haloes complying with all criteria.\n", 
					ThisTask, 
					SETTINGS->n_threshold, Settings.mass_min,
					SETTINGS->n_virialized,
					SETTINGS->n_concentration, 
					SETTINGS->n_spin,
					SETTINGS->n_all);

		fprintf(stdout, "\nTask=%d has box edges:\n \t\tX=%lf  |  %lf\n \t\tY=%lf  |  %lf\n \t\tZ=%lf  |  %lf\n", 
					ThisTask,
					SETTINGS->box.X[0][0], SETTINGS->box.X[0][1], 
					SETTINGS->box.X[1][0], SETTINGS->box.X[1][1], 
					SETTINGS->box.X[2][0], SETTINGS->box.X[2][1]);

}



void set_box(int i)
{
	int j=0;
	double x[3];
	struct halo *HALO;
	struct general_settings *SETTINGS;

#ifdef WITH_MPI
	HALO = pHaloes[ThisTask];
	SETTINGS = &pSettings[ThisTask];
#else
	HALO = Haloes;
	SETTINGS = &Settings;
#endif

	for(j=0; j<3; j++)
		x[j] = HALO[i].X[j];

		for(j=0; j<3; j++)
		{		
			if( x[j] < SETTINGS->box.X[j][0])
				SETTINGS->box.X[j][0] = x[j];

			if( x[j] > SETTINGS->box.X[j][1])
				SETTINGS->box.X[j][1] = x[j];
		}
}


	 /* Set additional halo properties derived from the basic ones */
void set_additional_halo_properties(int n)
{
	int i=0;
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

	HALO[n].shape = HALO[n].a[2]/HALO[n].a[0];
	HALO[n].triax = (pow2(HALO[n].a[0]) - pow2(HALO[n].a[1]))/(pow2(HALO[n].a[0]) - pow2(HALO[n].a[2]));

	HALO[n].ecc = sqrt(1 - 2*pow2(HALO[n].lambdaE));
	HALO[n].th_vir=2*HALO[n].Ekin/HALO[n].Epot;	
	HALO[n].abs_th_vir=sqrt(pow2(HALO[n].th_vir)); 

#ifdef GAS
		HALO[n].dm.N = HALO[n].n_part - HALO[n].gas.N;
		HALO[n].dm.M = HALO[n].Mvir - HALO[n].gas.M;
		HALO[n].gas_only.b_fraction = HALO[n].gas.M/HALO[n].Mvir;
#ifdef EXTRA_GAS
		HALO[n].gas_only.T = convert_u_to_T(HALO[n].gas_only.Cum_u);
		
		for(i=0; i<3; i++)
		{
			HALO[n].dm.X[i] = (HALO[n].Mvir*HALO[n].X[i] - HALO[n].gas.X[i])/HALO[n].dm.M;
			HALO[n].dm.V[i] = (HALO[n].Mvir*HALO[n].V[i] - HALO[n].gas.V[i])/HALO[n].dm.M;
		}
#else
		HALO[n].gas_only.T = 0.0; 
#endif // EXTRA_GAS
#endif // GAS
}



void get_halo_files_urls()
{
	int n=0, lin=0, lin_pro=0, lin_sub=0;
	char dummyline[200], url_fc[200], url_fc_pro[200], url_fc_sub[200];
	FILE *fc=NULL, *fc_pro=NULL, *fc_sub=NULL;

	struct internal_urls *URLS;


#ifdef WITH_MPI
	TASK_INFO_MSG(ThisTask, "reading halo file urls...");
	URLS = &pUrls[ThisTask];
#else
	INFO_MSG("Reading halo file urls...");
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

		URLS->urls = (char **) calloc(URLS->nCatalogueFiles, sizeof(char *));
		URLS->urls_profiles = (char **) calloc(URLS->nCatalogueFiles, sizeof(char *));
		URLS->urls_satellites = (char **) calloc(URLS->nCatalogueFiles, sizeof(char *));

	if(lin > 0)
	{
		for(n=0; n<URLS->nCatalogueFiles; n++)
		{
			fgets(dummyline,200,fc);
			URLS->urls[n] = (char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(URLS->urls[n], dummyline);
			URLS->urls[n][strlen(dummyline)-1]='\0';
		}
	fclose(fc);

		} else {
			if(ThisTask==0)
				WARNING("Trying to read empty file", url_fc);
		}

	if(lin_sub > 0)
	{
		for(n=0; n<URLS->nCatalogueFiles; n++)
		{
			fgets(dummyline,200,fc_sub);
			URLS->urls_satellites[n] = (char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(URLS->urls_satellites[n], dummyline);
			URLS->urls_satellites[n][strlen(dummyline)-1]='\0';
		}
	fclose(fc_sub);

		} else {
			if(ThisTask==0)
				WARNING("Trying to read empty file", url_fc_sub);
		}

	if(lin_pro > 0)
	{
		for(n=0; n<URLS->nCatalogueFiles; n++)
		{
			fgets(dummyline,200,fc_pro);
			URLS->urls_profiles[n] = (char*) calloc(strlen(dummyline), sizeof(char));
			strcpy(URLS->urls_profiles[n], dummyline);
			URLS->urls_profiles[n][strlen(dummyline)-1]='\0';
		}
	fclose(fc_pro);

		} else {
			if(ThisTask==0)
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

	URLS->halo_file = (char*) calloc(strlen(URLS->urls[n])+1, sizeof(char));
	strcpy(URLS->halo_file, URLS->urls[n]);

	URLS->profiles_file = (char*) calloc(strlen(URLS->urls_profiles[n])+1, sizeof(char));
	strcpy(URLS->profiles_file, URLS->urls_profiles[n]);
}



void read_halo_file()
{
	int n=0, i=0, j=0, thr=0, vir=0, conc=0, spin=0, skip=0, all=0, condition=0;
	double a=0; // Dummy variable to read useless columns
	char dummyline[LINE_SIZE]; 
	FILE *h_file;

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

		if(ThisTask==0)
			skip = 1;
		else
			skip = 0;

		h_file = fopen(URLS->halo_file, "r");
		SETTINGS->n_haloes = get_lines(h_file, URLS->halo_file) - skip;

#ifdef WITH_MPI
		pHaloes[ThisTask] = (struct halo*) calloc(SETTINGS->n_haloes, sizeof(struct halo));
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
	"%llu %llu %d  %lf %d  %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %d  %lf %lf %lf \
	 %lf %lf %lf \
	",
#else // There is a gas component
#ifndef EXTRA_GAS
	"%llu %llu %d  %lf %d  %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %d  %lf %lf %lf \
	 %lf %lf %lf \
	 %d  %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
	", 
#else	// Use extra gas columns
	"%llu %llu %d  %lf %d  %lf %lf %lf %lf %lf \
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
	&HALO[n].Mvir, 		&HALO[n].n_part,	&HALO[n].X[0], 	
	&HALO[n].X[1], 		&HALO[n].X[2], 		&HALO[n].V[0], 	
	&HALO[n].V[1], 		&HALO[n].V[2],		&HALO[n].Rvir, 
	&HALO[n].Rmax, 		&HALO[n].r2, 		&a, 				
	&a,			&HALO[n].Vmax, 		&a, 			
	&a, 			&HALO[n].lambda, 	// First 20 entries
	&HALO[n].lambdaE,	&HALO[n].L[0], 		&HALO[n].L[1], 			
	&HALO[n].L[2], 		&HALO[n].a[1],		&HALO[n].a[2], 		
	&HALO[n].Ea[0], 	&HALO[n].Ea[1],		&HALO[n].Ea[2],		
	&a, 			&a, 			&a, 
	&a, 			&a,			&a, 				
	&a,			&HALO[n].n_bins,	&a, 			
	&HALO[n].Ekin,	 	&HALO[n].Epot, 		// First 40 entries
	&a, 			&a,			&HALO[n].c_nfw
#ifdef GAS
	, &HALO[n].gas.N, 	&HALO[n].gas.M, 	&HALO[n].gas_only.lambda, 		
	&HALO[n].gas_only.lambdaE, 	&a,			&a, 			
	&a, 			&HALO[n].gas_only.a[1], 
	&HALO[n].gas_only.a[2], 	&HALO[n].gas_only.Ea[0],&HALO[n].gas_only.Ea[1], 		
	&HALO[n].gas_only.Ea[2], 	&a,			&a, 			
	&a,			&a,			&a, 			
	&a, 			&HALO[n].gas_only.Ekin, 	&HALO[n].gas_only.Epot
#ifdef EXTRA_GAS
	, &HALO[n].gas.X[0],	&HALO[n].gas.X[1], 	&HALO[n].gas.X[2], 		
	&HALO[n].gas.V[0], 	&HALO[n].gas.V[1], 	&HALO[n].gas.V[2],	
	&HALO[n].gas_only.Cum_u
#endif
#endif
	); 
#ifdef GAS
#ifndef EXTRA_GAS
	// The cumulative u has to be read from the profile catalogues
	HALO[n].Cum_u_gas = 0.0;
#endif
#endif
	// In the new catalogues haloes' major axis is normalized to one
	HALO[n].a[0] = 1.0; 

	set_additional_halo_properties(n);

#ifdef USE_UNIT_MPC
	HALO[n].Rvir *= 1.e-3;

	for(i=0; i<3; i++)
	{
		HALO[n].X[i] *= 1.e-3;
		HALO[n].V[i] *= 1.e-3;
#ifdef GAS
#ifdef EXTRA_GAS
		HALO[n].gas.X[i] *= 1.e-3;
		HALO[n].gas.V[i] *= 1.e-3;
		HALO[n].dm.X[i] *= 1.e-3;
		HALO[n].dm.V[i] *= 1.e-3;
#endif // Extra Gas
#endif // Gas
	}
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
	int skip=0, i=0, j=0, k=0; 
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

		if(p_file==NULL)
			ERROR("Profiles file not found", URLS->profiles_file);
		else
			fprintf(stdout, "Found profiles file:%s\n", 
				URLS->profiles_file);

		if(ThisTask==0)
			skip = 1;
		else
			skip = 0;

#ifdef WITH_MPI
		HALO = pHaloes[ThisTask];
#else
		HALO = Haloes;
#endif

		while(counter < SETTINGS->n_threshold)
		{
			fgets(dummyline, LINE_SIZE, p_file);

		if(j >= skip) 
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
		&a, 	 &a,     &a, &a,     &m_gas,&a,	     &u_gas  // 27
#endif
		);
				halo_bins = HALO[counter].n_bins;

				if(i == 0)
				{	
					HALO[counter].radius = (double *) calloc(halo_bins,sizeof(double));
					HALO[counter].rho = (double *) calloc(halo_bins,sizeof(double));
					HALO[counter].err = (double *) calloc(halo_bins,sizeof(double));
#ifdef GAS
					HALO[counter].gas_only.u = (double *) calloc(halo_bins,sizeof(double));
					HALO[counter].gas_only.m = (double *) calloc(halo_bins,sizeof(double));
#endif
				}

					HALO[counter].radius[i] = sqrt(pow2(radius));	
					HALO[counter].rho[i] = dens;
					HALO[counter].err[i] = dens/sqrt(npart-npart_old); 
#ifdef GAS
					HALO[counter].gas_only.u[i] = u_gas;
					HALO[counter].gas_only.m[i] = m_gas;
#endif
			//fprintf(stderr,dummyline);
			//fprintf(stderr, "\nHaloR  [%d][%d]=%e ", counter, i, HALO[counter].radius[i]);
			//fprintf(stderr, "HaloGas[%d][%d]=%f ", counter, i, HALO[counter].u_gas[i]);
			//fprintf(stderr, "HaloRho[%d][%d]=%f\n", counter, i, HALO[counter].rho[i]);

			if(HALO[counter].err[i] == 0.) 
				HALO[counter].err[i] = HALO[counter].err[i-1]; 

			if(radius<0.) 
				neg_r++;

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
