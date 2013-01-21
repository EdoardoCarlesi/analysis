#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>

#include "tinker.h"
#include "cosmological_relations.h"
#include "mass_function.h"
#include "power.h"
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
		fprintf(stderr, "\nWARNING **** Rho_c = 0.000 **** WARNING\n");

	return Settings.rho_c * 4 * PI * (1./3) * r * r * r;
}



double R_m(double m)
{
	if (Settings.rho_0 == 0) 
		fprintf(stderr, "\nWARNING **** Rho_0 = 0.000 **** WARNING\n");

	return pow((3 * m)/(Settings.rho_0 * 4 * PI), 1./3.);
}



double W_r(double k, double R)
{
	return 3*((sin(k*R) - k * R * cos(k*R))) * pow(k*R,-3);
}



double sigma_M(double m)
{
	return get_interpolated_value(ThMassFunc.th_masses, ThMassFunc.sigmas, ThMassFunc.bins, m);
}



double ln_sigma_M(double log_M, void *p)
{
	return get_interpolated_value(ThMassFunc.ln_masses, ThMassFunc.ln_sigmas, ThMassFunc.bins, log_M);
}



double d_ln_sigma_d_ln_M(double M)
{
	return get_interpolated_value(ThMassFunc.ln_masses, ThMassFunc.dlnSdlnM, ThMassFunc.bins, inv_ln_10*log(M));
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
	k_max = Pks[PK_INDEX].k[Pks[0].n_pk_entries-1];

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

	s8 = Cosmo.s8;
	par[0]=8.; 
	PK_INDEX=0;

	k_min = Pks[0].k[0]; 
	k_max = Pks[0].k[Pks[0].n_pk_entries-1];

		gsl_function F;
		F.function = &sigma_integ;
		F.params = &par;

		gsl_integration_qags(&F, k_min, k_max, 0, 1e-6, WORKSP, w, &result, &error);
		result*=(1./(2*PI*PI));

		Cosmo.n_s8 = s8*s8/result;

			fprintf(stdout, "Original sigma8:%lf, required sigma8: %lf, normalization factor:%lf\n", 
				sqrt(result), s8, Cosmo.n_s8);

		gsl_integration_workspace_free (w);

	normalize_all_power_spectra_to_sigma8();
}



void init_MassFunc()
{
	fprintf(stdout, "\nAllocating memory for mass function structures...\n");

	int nBins = Settings.n_bins;
	// We need N-1 points, since the value saved is the one at half length of every bin
	MassFunc.bins = nBins-1;
	MassFunc.num_masses = (double*) calloc(nBins-1, sizeof(double));
	MassFunc.dn_num_masses = (double*) calloc(nBins-1, sizeof(double));
	MassFunc.n = (double*) calloc(nBins-1, sizeof(double));
	MassFunc.n_tot = (int*) calloc(nBins-1, sizeof(int));
	MassFunc.n_bin = (int*) calloc(nBins-1, sizeof(int));
	MassFunc.dn  = (double*) calloc(nBins-1, sizeof(double));
	MassFunc.err = (double*) calloc(nBins-1, sizeof(double));
	MassFunc.err_dn = (double*) calloc(nBins-1, sizeof(double));
}



void init_ThMassFunc()
{
	fprintf(stdout, "\nAllocating memory for theoretical mass function structures...\n");

	int pts=ThMassFunc.bins;

	ThMassFunc.th_masses = (double *) calloc(pts, sizeof(double)); 
	ThMassFunc.sigmas = (double *) calloc(pts, sizeof(double)); 
	ThMassFunc.ln_masses = (double *) calloc(pts, sizeof(double)); 
	ThMassFunc.ln_sigmas = (double *) calloc(pts, sizeof(double)); 
	ThMassFunc.dlnSdlnM = (double *) calloc(pts, sizeof(double)); 
	ThMassFunc.tin = (double*) calloc(pts, sizeof(double));
	ThMassFunc.diff_tin = (double*) calloc(pts, sizeof(double));
}



void init_sigmas()
{
	int k=0, pts=ThMassFunc.bins;
	double M, Mmax, Mmin, sig, *m_steps; 

	m_steps = (double*) calloc(pts, sizeof(double));

	fprintf(stdout, "\nInitializing sigmas to calculate the theoretical Tinker mass function...");
	// ThMassFunc Min and Max are initialized when the program starts in the values read in the exec.sh script
	Mmin=ThMassFunc.Mmin;
	Mmax=ThMassFunc.Mmax;

	m_steps = log_stepper(Mmin, Mmax, pts); 

#ifdef _OPENMP
#endif
		for(k=0; k<ThMassFunc.bins; k++)
		{
			M=m_steps[k];
			sig = sqrt(integrate_sigma2(M));
			ThMassFunc.th_masses[k]=M;
			ThMassFunc.ln_masses[k]=inv_ln_10*log(M);
			ThMassFunc.sigmas[k]=sig;
			ThMassFunc.ln_sigmas[k]=inv_ln_10*log(sig);
		}

	fprintf(stdout, "done.\n");
	free(m_steps);
}



void init_dlnsdlnm()
{
	int k=0;
	double result, abserr, mm;

		gsl_function F;
		F.function = &ln_sigma_M;
		F.params = 0;

		for(k=0; k<ThMassFunc.bins; k++)
		{
			mm=ThMassFunc.ln_masses[k];
			gsl_deriv_central(&F, mm, 1e-4, &result, &abserr);
			ThMassFunc.dlnSdlnM[k]=result;
		}
}




double mf_normalization(double M)
{
	double lnslnm, norm;

		lnslnm = sqrt(pow(d_ln_sigma_d_ln_M(M),2)); 
		norm = Settings.rho_0/(M*M);

	return norm*lnslnm*sqrt(2./PI);
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
			MassFunc.num_masses[i]=mass_bin[i];
			MassFunc.dn_num_masses[i]=mass_bin[i]+halfstep;
			MassFunc.dn[i]=n_mass[i]/(volume*dM);
			MassFunc.n[i]=cum_n_mass[i]/volume;
			MassFunc.err_dn[i]=sqrt(n_mass[i])/(volume*dM);
			MassFunc.err[i]=sqrt(cum_n_mass[i])/volume;

			MassFunc.n_bin[i]=n_mass[i];
			MassFunc.n_tot[i]=cum_n_mass[i];
		}
	
		free(mass); 
		free(mass_bin); 
		free(n_mass); 
		free(cum_n_mass);

	fprintf(stdout, "\nSorting done.\n");
}



void compute_theoretical_mass_function(int index)
{
	int k=0;
	void *p=0; 
	double Mmin, Mmax;
		
	fprintf(stdout, "\ncompute_theoretical_mass_function()\n");
	
		PK_INDEX = index;
	
		init_ThMassFunc();

		init_sigmas();

		init_dlnsdlnm();

		if(Settings.fit==1) 
			best_fit_mf_tinker(MassFunc.num_masses, MassFunc.dn, MassFunc.err_dn, MassFunc.bins);

		Mmin = ThMassFunc.th_masses[0];
		Mmax = ThMassFunc.th_masses[ThMassFunc.bins-1];

		for(k=0; k<ThMassFunc.bins-1; k++)
		{
			Mmin=ThMassFunc.th_masses[k];
			ThMassFunc.diff_tin[k]=tinker(Mmin, p);
			ThMassFunc.tin[k]=integral_tinker(Mmin, Mmax);
		}
}



void initialize_mass_function_datastore()
{
	int numPkFiles=Settings.n_pk_files-1, k=0, dim = MassFunc.bins;

	MassFuncZ = 
		(struct mass_function*) calloc(numPkFiles, sizeof(struct mass_function));

		for(k=0; k<numPkFiles; k++)
		{
			MassFuncZ[k].bins = dim;
			MassFuncZ[k].th_masses  = (double*) calloc(dim, sizeof(double));
			MassFuncZ[k].tin = (double*) calloc(dim, sizeof(double));
			MassFuncZ[k].diff_tin = (double*) calloc(dim, sizeof(double));
#ifndef TH_ONLY
			MassFuncZ[k].num_masses  = (double*) calloc(dim, sizeof(double));
			MassFuncZ[k].n = (double*) calloc(dim, sizeof(double));
			MassFuncZ[k].dn = (double*) calloc(dim, sizeof(double));
#endif
		}
}



void store_mf(int m)
{
	int k=0;

	fprintf(stdout, "\nstore_mf(%d). Saving the (cumulative) mass functions into mf structure.\n",m);

	MassFuncZ[m].bins = ThMassFunc.bins;

	for(k=0; k<ThMassFunc.bins; k++)
	{
		MassFuncZ[m].th_masses[k] = ThMassFunc.th_masses[k];
		MassFuncZ[m].tin[k] = ThMassFunc.tin[k];
		MassFuncZ[m].diff_tin[k] = ThMassFunc.diff_tin[k];
#ifndef TH_ONLY
		MassFuncZ[m].num_masses[k] = MassFunc.num_masses[k];
		MassFuncZ[m].n[k] = MassFunc.n[k];
		MassFuncZ[m].dn[k] = MassFunc.dn[k];
#endif
	}
}
