#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>

#include "power.h"
#include "tinker.h"
#include "mass_function.h"
#include "cosmological_relations.h"

#include "../libmath/mathtools.h"
#include "../general_variables.h"
#include "../general_functions.h"

#ifdef _OPENMP
#include <omp.h>
#endif


int PK_INDEX;


double M_r(double r)
{
	if (Settings.rho_c == 0) 
		fprintf(stdout, "\nWARNING **** Rho_c = 0.000 **** WARNING\n");

	return Settings.rho_c * 4 * PI * (1./3) * r * r * r;
}



double R_m(double m)
{
	if (Settings.rho_0 == 0) 
		fprintf(stdout, "\nWARNING **** Rho_0 = 0.000 **** WARNING\n");

	return pow((3 * m)/(Settings.rho_0 * 4 * PI), 1./3.);
}



double W_r(double k, double R)
{
	return 3*((sin(k*R) - k * R * cos(k*R))) * pow(k*R,-3);
}



double sigmaM(double m)
{
	return get_interpolated_value(ThMassFunc.mass, ThMassFunc.sigma, ThMassFunc.bins, m);
}



double ln_sigma(double log_M, void *p)
{
	return get_interpolated_value(ThMassFunc.ln_mass, ThMassFunc.ln_sigma, ThMassFunc.bins, log_M);
}



double dln_sigma_dln_mass(double M)
{
	return get_interpolated_value(ThMassFunc.ln_mass, ThMassFunc.dln_sigma_dln_mass, 
		ThMassFunc.bins, inv_ln_10*log(M));
}



double sigma_integ(double k, void *params)
{ 
	double r, *pp;

		pp = (double*) params;
		r = pp[0];

	return k * k * pow(W_r(k,r),2.) * power_k(k, PK_INDEX);
}



double integrate_sigma2(double M)
{
	int WORKSP = 5000; 
	double k_max, k_min, result, error, par[2]; 

	par[0] = R_m(M);

	k_min = Pks[PK_INDEX].k[0];
	k_max = Pks[PK_INDEX].k[Pks[0].npts-1];

		gsl_function F;
		F.function = &sigma_integ;
		F.params = &par;

		gsl_integration_workspace *w = gsl_integration_workspace_alloc(WORKSP);
		gsl_integration_qags(&F, k_min, k_max, 0, 1e-6, WORKSP, w, &result, &error);
		gsl_integration_workspace_free (w);

		result*=(1./(2*PI*PI));

	return result;
}



void normalize_sigma8()
{
	int WORKSP = 5000; 
	double s8, k_min, k_max, result, error, par[2]; 
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(WORKSP);

	s8 = Cosmo.sigma8;
	par[0]=8.; 
	PK_INDEX=0;

	k_min = Pks[0].k[0]; 
	k_max = Pks[0].k[Pks[0].npts-1];

		gsl_function F;
		F.function = &sigma_integ;
		F.params = &par;

		gsl_integration_qags(&F, k_min, k_max, 0, 1e-6, WORKSP, w, &result, &error);
		result*=(1./(2*PI*PI));

		Cosmo.norm_sigma8 = s8*s8/result;

			fprintf(stdout, "Original sigma8:%lf, required sigma8: %lf, normalization factor:%lf\n", 
				sqrt(result), s8, Cosmo.norm_sigma8);

		gsl_integration_workspace_free (w);

	normalize_all_power_spectra_to_sigma8();
}



void init_MassFunc()
{
	int nBins = MassFunc.bins;

	fprintf(stdout, "\nAllocating memory for mass function structures...\n");

	MassFunc.mass = (double*) calloc(nBins, sizeof(double));
	MassFunc.mass_halfstep = (double*) calloc(nBins-1, sizeof(double));
	MassFunc.n = (double*) calloc(nBins-1, sizeof(double));
	MassFunc.n_tot = (int*) calloc(nBins-1, sizeof(int));
	MassFunc.n_bin = (int*) calloc(nBins-1, sizeof(int));
	MassFunc.dn  = (double*) calloc(nBins-1, sizeof(double));
	MassFunc.err = (double*) calloc(nBins-1, sizeof(double));
	MassFunc.err_dn = (double*) calloc(nBins-1, sizeof(double));
}



void init_ThMassFunc()
{
	int nTot=ThMassFunc.bins;

	fprintf(stdout, "\nAllocating memory for theoretical mass function structures...\n");

	ThMassFunc.mass = (double *) calloc(nTot, sizeof(double)); 
	ThMassFunc.sigma = (double *) calloc(nTot, sizeof(double)); 
	ThMassFunc.ln_mass = (double *) calloc(nTot, sizeof(double)); 
	ThMassFunc.ln_sigma = (double *) calloc(nTot, sizeof(double)); 
	ThMassFunc.dln_sigma_dln_mass = (double *) calloc(nTot, sizeof(double)); 
	ThMassFunc.n = (double*) calloc(nTot, sizeof(double));
	ThMassFunc.dn = (double*) calloc(nTot, sizeof(double));
}


void init_sigmas()
{
	int k=0, nTot=0; 
	double M=0, Mmax=0, Mmin=0, sig=0; 
	double *m_steps=NULL; 

	nTot = ThMassFunc.bins;
	m_steps = (double*) calloc(nTot, sizeof(double));

	fprintf(stdout, "\nInitializing sigmas to calculate the theoretical Tinker mass function...");

	// ThMassFunc Min and Max are initialized when the program starts in the values read in the exec.sh script
	Mmin=ThMassFunc.Mmin;
	Mmax=ThMassFunc.Mmax;

	m_steps = log_stepper(Mmin, Mmax, nTot); 

#ifdef _OPENMP
		omp_set_num_threads(OMP_THREADS);
#endif

#		pragma omp parallel for 		\
		shared(ThMassFunc) private(k, M, sig)
		for(k=0; k<nTot; k++)
		{
			M=m_steps[k];
			sig=sqrt(integrate_sigma2(M));
			ThMassFunc.mass[k]=M;
			ThMassFunc.sigma[k]=sig;
			ThMassFunc.ln_mass[k]=inv_ln_10*log(M);
			ThMassFunc.ln_sigma[k]=inv_ln_10*log(sig);
		}

	fprintf(stdout, "done.\n");
	free(m_steps);
}



void init_dln_sigma_dln_mass()
{
	int k=0, nTot=0;
	double result, abserr, mass;

	nTot = ThMassFunc.bins;

		gsl_function F;
		F.function = &ln_sigma;
		F.params = 0;
	
#ifdef _OPENMP
		omp_set_num_threads(OMP_THREADS);
#endif

#		pragma omp parallel for 		\
		shared(ThMassFunc) private(mass, k, result)
		for(k=0; k<nTot; k++)
		{
			mass = ThMassFunc.ln_mass[k];

			if	(k==0)
				gsl_deriv_forward  (&F, mass, 1e-4, &result, &abserr);
			else if (k==nTot-1)
				gsl_deriv_backward (&F, mass, 1e-4, &result, &abserr);
			else
				gsl_deriv_central  (&F, mass, 1e-4, &result, &abserr);
			
			ThMassFunc.dln_sigma_dln_mass[k]=result;
		}
}



double mf_normalization(double M)
{
	return sqrt(2./PI) * Settings.rho_0/(M*M) * sqrt(pow(dln_sigma_dln_mass(M),2));
}



void initialize_mass_function_datastore()
{
	int numPkFiles=0, k=0, nTot=0; 
	
	numPkFiles = Urls.nPkFiles;
	MassFuncZ = (struct mass_function*) calloc(numPkFiles, sizeof(struct mass_function));

		for(k=0; k<numPkFiles; k++)
		{
			nTot = ThMassFunc.bins;
			ThMassFuncZ[k].bins = nTot;
			ThMassFuncZ[k].mass = (double*) calloc(nTot, sizeof(double));
			ThMassFuncZ[k].n = (double*) calloc(nTot, sizeof(double));
			ThMassFuncZ[k].dn = (double*) calloc(nTot, sizeof(double));
#ifndef TH_ONLY
			nTot = MassFunc.bins;
			MassFuncZ[k].bins = nTot;	
			MassFuncZ[k].mass = (double*) calloc(nTot, sizeof(double));
			MassFuncZ[k].n = (double*) calloc(nTot, sizeof(double));
			MassFuncZ[k].dn = (double*) calloc(nTot, sizeof(double));
#endif
		}
}



void store_mf(int m)
{
	int k=0, nTot=0;

	fprintf(stdout, "\nstore_mf(%d). Saving the (cumulative) mass functions into mf structure.\n",m);

	nTot = ThMassFunc.bins;

	MassFuncZ[m].bins = nTot;

	for(k=0; k<nTot; k++)
	{
		ThMassFuncZ[m].mass[k] = ThMassFunc.mass[k];
		ThMassFuncZ[m].n[k] = ThMassFunc.n[k];
		ThMassFuncZ[m].dn[k] = ThMassFunc.dn[k];
	}

#ifndef TH_ONLY
	nTot = ThMassFunc.bins;
	ThMassFuncZ[m].bins = nTot;

	for(k=0; k<nTot; k++)
	{
		MassFuncZ[m].mass[k] = MassFunc.mass[k];
		MassFuncZ[m].n[k] = MassFunc.n[k];
		MassFuncZ[m].dn[k] = MassFunc.dn[k];
	}
#endif
}



void compute_theoretical_mass_function()
{
	int k=0, nTot=0;
	void *p=0; 
	double Mmin, Mmax;
//	double tstart, tend, dt;	
	
	nTot = ThMassFunc.bins;

	fprintf(stdout, "\ncompute_theoretical_mass_function()\n");
	
		PK_INDEX = Settings.use_cat;

//		tstart = time(NULL); 
		init_ThMassFunc();
		init_sigmas();
		init_dln_sigma_dln_mass();

		if(Settings.fit==1) 
			best_fit_mf_tinker(MassFunc.mass, MassFunc.dn, MassFunc.err_dn, MassFunc.bins);

		Mmax = ThMassFunc.mass[nTot-1];

#ifdef _OPENMP
		omp_set_num_threads(OMP_THREADS);
#endif

#		pragma omp parallel for 		\
		shared(ThMassFunc) private(Mmin, k)
		for(k=0; k<nTot; k++)
		{
			Mmin=ThMassFunc.mass[k];
			ThMassFunc.dn[k]=tinker(Mmin, p);
			ThMassFunc.n[k]=integral_tinker(Mmin, Mmax);
		}
		
//		tend=time(NULL);
//		dt = tend - tstart;
//	fprintf(stdout, "\nTime Of Run=%lf\n", dt);
}



void compute_numerical_mass_function(void)
{
	int nBins=Settings.n_bins, nHaloes=Settings.n_haloes, i=0; 
	double dn_norm=1., volume=0, mMin=0, mMax=0, halfstep=0, dM=0; 
	double *mass, *mass_bin; 
	int *n_mass, *cum_n_mass;

	fprintf(stdout, "\nSorting mass function for %d halos in %d bins\n", nHaloes, nBins);
	
		Settings.tick=0;
	
		mass = (double*) calloc(nHaloes, sizeof(double));
		mass_bin = (double*) calloc(nBins, sizeof(double));
		n_mass = (int*) calloc(nBins-1, sizeof(int));
		cum_n_mass = (int*) calloc(nBins-1, sizeof(int));

		init_MassFunc();

#		pragma omp parallel for 	\
		shared(Haloes, mass) private(i)			
		for(i=0; i<nHaloes; i++)
			mass[i] = Haloes[i].Mvir;			
	
		mMin=minimum(mass, nHaloes)*0.999;
		mMax=maximum(mass, nHaloes)*1.001;
		mass_bin = log_stepper(mMin, mMax, nBins);
	
		lin_bin(mass, mass_bin, nBins, nHaloes, n_mass);	
		
		cum_bin(n_mass, cum_n_mass, nBins-1);

			volume=Settings.box_size*Settings.box_size*Settings.box_size;

		for(i=0; i<nBins-1; i++)
		{
			halfstep = 0.5*(mass_bin[i+1]-mass_bin[i]);
			dn_norm = 2*halfstep/nHaloes;
			dM = mass_bin[i+1]-mass_bin[i];
			MassFunc.mass[i]=mass_bin[i];
			MassFunc.mass_halfstep[i]=mass_bin[i]+halfstep;
			MassFunc.dn[i]=n_mass[i]/(volume*dM);
			MassFunc.n[i]=cum_n_mass[i]/volume;
			MassFunc.err_dn[i]=sqrt(n_mass[i])/(volume*dM);
			MassFunc.err[i]=sqrt(cum_n_mass[i])/volume;

			MassFunc.n_bin[i]=n_mass[i];
			MassFunc.n_tot[i]=cum_n_mass[i];
		}
	
		free(cum_n_mass);
		free(mass_bin); 
		free(n_mass); 
		free(mass); 

	fprintf(stdout, "\nSorting done.\n");
}
