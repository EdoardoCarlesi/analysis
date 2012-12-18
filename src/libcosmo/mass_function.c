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

#ifdef KURAC
#include "tareador_hooks.h"
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
	if (Settings.rho_0 == 0) fprintf(stderr, "\nWARNING **** Rho_0 = 0.000 **** WARNING\n");

	return pow((3 * m)/(Settings.rho_0 * 4 * PI), 1./3.);
}



double W_r(double k, double R)
{
	return 3*((sin(k*R) - k * R * cos(k*R))) * pow(k*R,-3);
}



double sigma_M(double m)
{
	return get_interpolated_value(AMF.th_masses, AMF.sigmas, AMF.bins, m);
}



double ln_sigma_M(double log_M, void *p)
{
	return get_interpolated_value(AMF.ln_masses, AMF.ln_sigmas, AMF.bins, log_M);
}



double d_ln_sigma_d_ln_M(double M)
{
	return get_interpolated_value(AMF.ln_masses, AMF.dlnSdlnM, AMF.bins, inv_ln_10*log(M));
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

			fprintf(stderr, "Original sigma8:%lf, required sigma8: %lf, normalization factor:%lf\n", 
				sqrt(result), s8, Cosmo.n_s8);

		gsl_integration_workspace_free (w);

	normalize_all_power_spectra_to_sigma8();
}



void init_MF()
{
	int nBins = Settings.n_bins;
	// We need N-1 points, since the value saved is the one at half length of every bin
	MF.bins = nBins-1;
	MF.num_masses = (double*) calloc(nBins-1, sizeof(double));
	MF.n = (double*) calloc(nBins-1, sizeof(double));
	MF.n_tot = (int*) calloc(nBins-1, sizeof(int));
	MF.n_bin = (int*) calloc(nBins-1, sizeof(int));
	MF.dn  = (double*) calloc(nBins-1, sizeof(double));
	MF.err = (double*) calloc(nBins-1, sizeof(double));
	MF.err_dn = (double*) calloc(nBins-1, sizeof(double));
}



void init_AMF()
{
	int pts=AMF.bins;

	AMF.th_masses = (double *) calloc(pts, sizeof(double)); 
	AMF.sigmas = (double *) calloc(pts, sizeof(double)); 
	AMF.ln_masses = (double *) calloc(pts, sizeof(double)); 
	AMF.ln_sigmas = (double *) calloc(pts, sizeof(double)); 
	AMF.dlnSdlnM = (double *) calloc(pts, sizeof(double)); 
	AMF.tin = (double*) calloc(pts, sizeof(double));
	AMF.diff_tin = (double*) calloc(pts, sizeof(double));
}



void init_sigmas()
{
	int k=0, pts=AMF.bins;
	double M, Mmax, Mmin, sig, *m_steps; 

	fprintf(stderr, "\nInitializing sigmas to calculate the theoretical Tinker mass function.\n");
	// AMF Min and Max are initialized when the program starts in the values read in the exec.sh script
	Mmin=AMF.Mmin;
	Mmax=AMF.Mmax;

	m_steps = log_stepper(Mmin, Mmax, pts); 

		for(k=0; k<AMF.bins; k++)
		{
			M=m_steps[k];
			sig = sqrt(integrate_sigma2(M));
			AMF.th_masses[k]=M;
			AMF.ln_masses[k]=inv_ln_10*log(M);
			AMF.sigmas[k]=sig;
			AMF.ln_sigmas[k]=inv_ln_10*log(sig);
		}

	free(m_steps);
}


void init_dlnsdlnm()
{
	int k=0;
	double result, abserr, mm;

		gsl_function F;
		F.function = &ln_sigma_M;
		F.params = 0;

		for(k=0; k<AMF.bins; k++)
		{
			mm=AMF.ln_masses[k];
			gsl_deriv_central(&F, mm, 1e-4, &result, &abserr);
			AMF.dlnSdlnM[k]=result;
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
	int nBins=Settings.n_bins, nHaloes=Settings.n_haloes, n=0, i=0, k=0; 
	double cumHalo=0, volume=0, mMin=0, mMax=0, halfstep=0, *mass_bins;

	fprintf(stderr, "\nSorting mass function for %d halos in %d bins\n", nHaloes, nBins);
	
		Settings.tick=0;
	
		mass = (double*) calloc(totSub, sizeof(double));
		mass_bin = (double*) calloc(nBins, sizeof(double));
		n_mass = (int*) calloc(nBins-1, sizeof(int));
		cum_n_mass = (int*) calloc(nBins-1, sizeof(int));

		init_MF();
	
		mMin=haloes[nHaloes-1].Mvir;
		mMax=haloes[0].Mvir;
		mass_bin = log_stepper(mMin, mMax, nBins);
		volume=Settings.box_size*Settings.box_size*Settings.box_size;

		lin_bin(mass, mass_bin, nBins, totSub, n_mass);	


			for(i=0; i<nBins-1; i++)
			{
				k = nBins - i - 1;
				cumHalo += n_mass[k];
				cum_n_mass[k] = cumHalo;
			}


		for(i=0; i<nBins-1; i++)
		{
				halfstep = 0.5*(mass_bin[i+1]-mass_bin[i]);
				MF.num_masses[i]=mass_bin[i]+halfstep;
				MF.n[i]=n_mass[i];
				MF.n_tot[i]=cum_n_mass[i];
				MF.err[i]=sqrt(n_mass[i])/volume;
				MF.dn[i]=n_mass[i]/volume;
				MF.n_bin[i]=cum_n_mass[i]/volume;
				MF.err_dn[i]=sqrt(cum_n_mass[i])/volume;
		}	

		free(mass); 
		free(mass_bin); 
		free(n_mass); 
		free(cum_n_mass);
}



void compute_theoretical_mass_function(int index)
{
	int k=0;
	void *p=0; 
	double Mmin, Mmax;
		
	fprintf(stdout, "\ncompute_theoretical_mass_function()\n");
	
		PK_INDEX = index;
	
		init_AMF();

		init_sigmas();

		init_dlnsdlnm();

#ifndef TH_ONLY
			if(Settings.fit==1) 
				best_fit_mf_tinker(MF.num_masses, MF.dn, MF.err_dn, MF.bins);
#endif

		Mmin = AMF.th_masses[0];
		Mmax = AMF.th_masses[AMF.bins-1];

		for(k=0; k<AMF.bins-1; k++)
		{
			Mmin=AMF.th_masses[k];
			AMF.diff_tin[k]=tinker(Mmin, p);
			AMF.tin[k]=integral_tinker(Mmin, Mmax);
		}
}



void initialize_mass_function_datastore()
{
	int numPkFiles=Settings.n_pk_files-1, k=0, dim = MF.bins;

	mf = (struct mass_function*) calloc(numPkFiles, sizeof(struct mass_function));

		for(k=0; k<numPkFiles; k++)
		{
			mf[k].bins = dim;
			mf[k].th_masses  = (double*) calloc(dim, sizeof(double));
			mf[k].tin = (double*) calloc(dim, sizeof(double));
			mf[k].diff_tin = (double*) calloc(dim, sizeof(double));
#ifndef TH_ONLY
			mf[k].num_masses  = (double*) calloc(dim, sizeof(double));
			mf[k].n = (double*) calloc(dim, sizeof(double));
			mf[k].dn = (double*) calloc(dim, sizeof(double));
#endif
		}
}



void store_mf(int m)
{
	int k=0;

	fprintf(stderr, "\nstore_mf(%d). Saving the (cumulative) mass functions into mf structure.\n",m);

	mf[m].bins = AMF.bins;

	for(k=0; k<AMF.bins; k++)
	{
		mf[m].th_masses[k] = AMF.th_masses[k];
		mf[m].tin[k] = AMF.tin[k];
		mf[m].diff_tin[k] = AMF.diff_tin[k];
#ifndef TH_ONLY
		mf[m].num_masses[k] = MF.num_masses[k];
		mf[m].n[k] = MF.n[k];
		mf[m].dn[k] = MF.dn[k];
#endif
	}
}
