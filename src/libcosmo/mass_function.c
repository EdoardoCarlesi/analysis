#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>

#include "../libmath/math.h"
#include "../general_def.h"

#include "cosmo.h"

#ifdef _OPENMP
#include <omp.h>
#endif


/*
 * Declare functions
 */
int PK_INDEX;
int MF_INDEX;

double M_r(double);
double W_r(double, double);
double R_m(double);

double sigma_integ(double, void*);
double integrate_sigma2(double);
void normalize_sigma8(void);
void store_mf(int);
double mf_normalization(double);

void init_ThMassFunc(void);
void init_sigmas(void);
void init_dln_sigma_dln_mass(void);

double sigmaM(double);
double ln_sigma(double, void*);
double dln_sigma_dln_mass(double);


/*
 * Initialize functions
 */ 
double M_r(double r)
{
	if (Settings.rho_c == 0) 
		WARNING("RhoCritical set to zero", "M_r()");

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
	return get_interpolated_value(ThMassFunc[MF_INDEX].mass, ThMassFunc[MF_INDEX].sigma, ThMassFunc[MF_INDEX].bins, m);
}



double ln_sigma(double log_M, void *p)
{
	return get_interpolated_value(ThMassFunc[MF_INDEX].ln_mass, ThMassFunc[MF_INDEX].ln_sigma, ThMassFunc[MF_INDEX].bins, log_M);
}



double dln_sigma_dln_mass(double M)
{
	return get_interpolated_value(ThMassFunc[MF_INDEX].ln_mass, ThMassFunc[MF_INDEX].dln_sigma_dln_mass, 
		ThMassFunc[MF_INDEX].bins, inv_ln_10*log(M));
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



void init_ThMassFunc()
{
	int nTot=ThMassFunc[MF_INDEX].bins;

	fprintf(stdout, "\nAllocating memory for theoretical mass function structures...\n");

	ThMassFunc[MF_INDEX].mass = (double *) calloc(nTot, sizeof(double)); 
	ThMassFunc[MF_INDEX].sigma = (double *) calloc(nTot, sizeof(double)); 
	ThMassFunc[MF_INDEX].ln_mass = (double *) calloc(nTot, sizeof(double)); 
	ThMassFunc[MF_INDEX].ln_sigma = (double *) calloc(nTot, sizeof(double)); 
	ThMassFunc[MF_INDEX].dln_sigma_dln_mass = (double *) calloc(nTot, sizeof(double)); 
	ThMassFunc[MF_INDEX].n = (double*) calloc(nTot, sizeof(double));
	ThMassFunc[MF_INDEX].dn = (double*) calloc(nTot, sizeof(double));
}


void init_sigmas()
{
	int k=0, nTot=0; 
	double M=0, Mmax=0, Mmin=0, sig=0; 
	double *m_steps=NULL; 

	nTot = ThMassFunc[MF_INDEX].bins;
	m_steps = (double*) calloc(nTot, sizeof(double));

	fprintf(stdout, "\nInitializing sigmas to calculate the theoretical Tinker mass function...");

	// ThMassFunc[MF_INDEX] Min and Max are initialized when the program starts in the values read in the exec.sh script
	Mmin=ThMassFunc[MF_INDEX].Mmin;
	Mmax=ThMassFunc[MF_INDEX].Mmax;

	m_steps = log_stepper(Mmin, Mmax, nTot); 

#		pragma omp parallel for 		\
		shared(ThMassFunc) private(k, M, sig)
		for(k=0; k<nTot; k++)
		{
			M=m_steps[k];
			sig=sqrt(integrate_sigma2(M));
			ThMassFunc[MF_INDEX].mass[k]=M;
			ThMassFunc[MF_INDEX].sigma[k]=sig;
			ThMassFunc[MF_INDEX].ln_mass[k]=inv_ln_10*log(M);
			ThMassFunc[MF_INDEX].ln_sigma[k]=inv_ln_10*log(sig);
		}

	fprintf(stdout, "done.\n");
	free(m_steps);
}



void init_dln_sigma_dln_mass()
{
	int k=0, nTot=0;
	double result, abserr, mass;

	nTot = ThMassFunc[MF_INDEX].bins;

		gsl_function F;
		F.function = &ln_sigma;
		F.params = 0;

#		pragma omp parallel for 		\
		shared(ThMassFunc) private(mass, k, result)
		for(k=0; k<nTot; k++)
		{
			mass = ThMassFunc[MF_INDEX].ln_mass[k];

			if	(k==0)
				gsl_deriv_forward  (&F, mass, 1e-4, &result, &abserr);
			else if (k==nTot-1)
				gsl_deriv_backward (&F, mass, 1e-4, &result, &abserr);
			else
				gsl_deriv_central  (&F, mass, 1e-4, &result, &abserr);
			
			ThMassFunc[MF_INDEX].dln_sigma_dln_mass[k]=result;
		}
}



double mf_normalization(double M)
{
	return sqrt(2./PI) * Settings.rho_0/(M*M) * sqrt(pow(dln_sigma_dln_mass(M),2));
}



void initialize_mass_function_structs()
{
	int numPkFiles=0, k=0, nTot=0; 
	
	numPkFiles = Urls.nPkFiles;
	ThMassFunc = (struct mass_function*) calloc(numPkFiles, sizeof(struct mass_function));
	MassFunc = (struct mass_function*) calloc(numPkFiles, sizeof(struct mass_function));

}



void compute_theoretical_mass_function()
{
	int k=0, nTot=0;
	void *p=0; 
	double Mmin, Mmax;
//	double tstart, tend, dt;	
	
	nTot = ThMassFunc[MF_INDEX].bins;

	fprintf(stdout, "\ncompute_theoretical_mass_function()\n");
	
		PK_INDEX = Settings.use_cat;

//		tstart = time(NULL); 
		init_ThMassFunc();
		init_sigmas();
		init_dln_sigma_dln_mass();

		if(Settings.fit==1) 
		best_fit_mf_tinker(MassFunc[MF_INDEX].mass, MassFunc[MF_INDEX].dn, MassFunc[MF_INDEX].err_dn, MassFunc[MF_INDEX].bins);

		Mmax = ThMassFunc[MF_INDEX].mass[nTot-1];

#		pragma omp parallel for 		\
		shared(ThMassFunc) private(Mmin, k)
		for(k=0; k<nTot; k++)
		{
			Mmin=ThMassFunc[MF_INDEX].mass[k];
			ThMassFunc[MF_INDEX].dn[k]=tinker(Mmin, p);
			ThMassFunc[MF_INDEX].n[k]=integral_tinker(Mmin, Mmax);
		}
		
//		tend=time(NULL);
//		dt = tend - tstart;
//	fprintf(stdout, "\nTime Of Run=%lf\n", dt);
}
