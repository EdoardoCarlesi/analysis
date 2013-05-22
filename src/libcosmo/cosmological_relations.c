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


/*
 * Declare functions
 */
double comoving_distance(double, double);
double comoving_vol(double, void*);
double integrate_comoving_volume(double, double);
double H_z(double, void*);
double inv_H_z(double, void*);


/*
 * Initialize functions
 */ 
double mass_temperature(double M)
{
	double alpha, log_M0;
	
		alpha=1.62;
		log_M0 = 14.56;

		fprintf(stdout, "M: %e, M0: %e \n", M, exp(log_M0));

	return 3*pow(M/pow(10,log_M0),1./alpha);
}



double w_z(double z)
{
	return get_interpolated_value(Cosmo.z_hub, Cosmo.w, Cosmo.npts, z);
}



double H_z(double z, void *p)
{
	double hz;

#ifdef USE_LCDM
		hz = sqrt(Cosmo.OmegaL*pow(z+1,0) + Cosmo.OmegaM*pow(z+1,3));
#else
		hz = get_interpolated_value(Cosmo.z_hub, Cosmo.Hubble, Cosmo.npts, z);
#endif

	return hz;
}



double inv_H_z(double z, void *p)
{
	return 1./H_z(z,p);
}



	/* Get the rho0 in SM / Mpc^3 getting rid of h */
double default_rho0()
{
	double rho;

		rho = 6.979e+7 * pow(512, 3.0) * 0.7 * 0.7 * pow(50, -3.0);
		fprintf(stdout, "\nSetting default rho_0() density to:%e in SolarMasses * h^-1 * (Mpc/h)^-3 units\n", rho);

	return rho;
}



double comoving_distance(double z0, double z1)
{
	double cd, fac, result, error;

	fac = D_H * pow(Cosmo.h,-1);
	gsl_function H;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	H.function=&inv_H_z;
	H.params=0;

		gsl_integration_qags(&H, z0, z1, 0, 1e-3, 1000, w, &result, &error);
		cd = result*fac;
		gsl_integration_workspace_free (w);

	return cd;
}



double comoving_vol(double z, void *param)
{
	double dV, z0, *pp;

		pp = (double*) param;
		z0 = pp[0];
		dV = pow(Cosmo.h,-1) * D_H * inv_H_z(z, param) * comoving_distance(z0,z)*comoving_distance(z0,z);

	return dV;
}



double integrate_comoving_volume(double z0, double z1)
{
	double result, error;

	gsl_function F;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	F.function=&comoving_vol;
	F.params=&z0;

		gsl_integration_qags(&F, z0, z1, 0, 1e-5, 1000, w, &result, &error);
		gsl_integration_workspace_free (w);

	return result;
}
