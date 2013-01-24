#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../general_variables.h"
#include "../general_functions.h"

// Define some useful macros
#define FILE_HEADER(fout, info, count) fprintf(fout, info "(%i)\t",count++)

#define DUMP_MSG(name, url) fprintf(stdout, "Printing %s to file: %s\n", name, url)


// These variables are used throughout all functions
static int i;
static int count;
static int nTot;
static double z;
static double M;
static char out_url[200];
static FILE *out_file;


void print_number_densities()
{
	count = 1;
	nTot = NumDen.npts; 
	M = Settings.mass_min;
	sprintf(out_url, "%sM_%.2e%s", Urls.output_prefix, M, "_number_density.dat");
	out_file = fopen(out_url, "w");

	DUMP_MSG("number density", out_url);

		fprintf(out_file,"#");
		FILE_HEADER(out_file, "z          ", count);
		FILE_HEADER(out_file, "Tinker n   ", count);
		FILE_HEADER(out_file, "numerical n", count);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++) 
			{
				fprintf(out_file,"%.4lf", NumDen.z[i]);
				fprintf(out_file,"\t%.3e", NumDen.n_th[i]);
				fprintf(out_file,"\t%.3e", NumDen.n_num[i]);
				fprintf(out_file, "\n");
			}

	fclose(out_file);
}



void print_all_subhalo_properties_to_one_file()
{
	count = 1;
	nTot = SubHaloZ.n_bins;
	z = GrowthFac.z[Settings.use_cat];
 	sprintf(out_url, "%sz%.2f%s", Urls.output_prefix, z, "_all_sub_statistics.dat");
	out_file = fopen(out_url,"w");

	DUMP_MSG("subhalo properties", out_url);

		fprintf(out_file,"#");
		FILE_HEADER(out_file, "lambda ", count);
		FILE_HEADER(out_file, "P(l)   ", count);
		FILE_HEADER(out_file, "shape  ", count);
		FILE_HEADER(out_file, "P(s)   ", count);
		FILE_HEADER(out_file, "triax  ", count);
		FILE_HEADER(out_file, "P(t)   ", count);
		FILE_HEADER(out_file, "Cos(th)", count);
		FILE_HEADER(out_file, "P(c(th))", count);
		FILE_HEADER(out_file, "Cos(ph) ", count);
		FILE_HEADER(out_file, "P(c(ph))", count);
		FILE_HEADER(out_file, "sub(r)/Rv", count);
		FILE_HEADER(out_file, "n(r)   ", count);
		FILE_HEADER(out_file, "n(>r)  ", count);
		FILE_HEADER(out_file, "V_sub  ", count);
		FILE_HEADER(out_file, "P(V_sub)", count);
		FILE_HEADER(out_file, "mass   ", count);
		FILE_HEADER(out_file, "n(>M)  ", count);
		FILE_HEADER(out_file, "subset(r)", count);
		FILE_HEADER(out_file, "subset(n(>r))", count);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)	
			{
				fprintf(out_file, "%lf"  , SubHaloZ.l[i]);
				fprintf(out_file, "\t%lf", SubHaloZ.p_l[i]);
				fprintf(out_file, "\t%lf", SubHaloZ.shape[i]);
				fprintf(out_file, "\t%lf", SubHaloZ.triax[i]);
				fprintf(out_file, "\t%lf", SubHaloZ.costh[i]);
				fprintf(out_file, "\t%lf", SubHaloZ.costh_count[i]);
				fprintf(out_file, "\t%lf", SubHaloZ.cosphi[i]);
				fprintf(out_file, "\t%lf", SubHaloZ.cosphi_count[i]);
				fprintf(out_file, "\t%lf", SubHaloZ.r_sub[i]);
				fprintf(out_file, "\t%d" , SubHaloZ.n_r_sub[i]);
				fprintf(out_file, "\t%d" , SubHaloZ.cum_n_r_sub[i]);
				fprintf(out_file, "\t%lf", SubHaloZ.vel_sub[i]);
				fprintf(out_file, "\t%lf", SubHaloZ.p_vel_sub[i]);
				fprintf(out_file, "\t%e" , SubHaloZ.mass_sub[i]);
				fprintf(out_file, "\t%d" , SubHaloZ.cum_n_sub[i]);
				fprintf(out_file, "\t%lf", SubHaloZ.r_sub_subset[i]);
				fprintf(out_file, "\t%d" , SubHaloZ.cum_n_r_sub_subset[i]);
				fprintf(out_file, "\n");
			}

	fclose(out_file);	
}



void print_theoretical_mass_function()
{
	count = 1;
	nTot = ThMassFunc.bins-1;
	z = GrowthFac.z[Settings.use_cat];
	sprintf(out_url, "%sz%.2f%s", Urls.output_prefix, z, "_theoretical_mass_function.dat");
	out_file = fopen(out_url, "w");

	DUMP_MSG("theoretical mass function", out_url);

		fprintf(out_file,"#");
		FILE_HEADER(out_file, "Mass  ", count);
		FILE_HEADER(out_file, "n(>M) ", count);
		FILE_HEADER(out_file, "dn    ", count);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)
			{
				fprintf(out_file,"%e", ThMassFunc.mass[i]);
				fprintf(out_file,"\t%e", ThMassFunc.n[i]);
				fprintf(out_file,"\t%e", ThMassFunc.dn[i]);
				fprintf(out_file,"\n");
			}

	fclose(out_file);
}



void print_correlation_function()
{
	count = 1;
	nTot = Xi.npts;
	z = GrowthFac.z[Settings.use_cat];
	sprintf(out_url, "%sz%.2f%s", Urls.output_prefix, z, "_correlation_function.dat");
	out_file = fopen(out_url,"w");

	DUMP_MSG("correlation function", out_url);

		fprintf(out_file,"#");
		FILE_HEADER(out_file, "R", count);
		FILE_HEADER(out_file, "Xi", count);
		FILE_HEADER(out_file, "Xi fit", count);
		fprintf(out_file, "\n");
	
			for(i=0; i<nTot; i++) 
			{
				fprintf(out_file, "%lf\n", Xi.r[i]); 
				fprintf(out_file, "\t%lf\n", Xi.xi_r[i]);
				fprintf(out_file, "\t%lf\n", Xi.xi_fit[i]);
				fprintf(out_file, "\n");
			}

	fclose(out_file);
}



void print_growth_factor()
{
	count = 1;
	nTot = GrowthFac.npts;
	sprintf(out_url, "%sk%.3f%s", Urls.output_prefix, GrowthFac.scale_k, "_growth_factor.dat");	
	out_file = fopen(out_url,"w");

	DUMP_MSG("growth factor", out_url);
	
		fprintf(out_file, "#");
		FILE_HEADER(out_file, "z     ", count);
		FILE_HEADER(out_file, "gf    ", count);
		FILE_HEADER(out_file, "gf/a  ", count);
		fprintf(out_file, "\n");
	
			for (i=0; i<nTot; i++) 
			{
				fprintf(out_file,"%lf", GrowthFac.z[i]);
				fprintf(out_file,"\t%lf", GrowthFac.gf[i]);
				fprintf(out_file,"\t%lf", GrowthFac.gf[i]/GrowthFac.a[i]);
				fprintf(out_file,"\n");
			}

	fclose(out_file);
}



void print_evolution_to_file()
{
	count = 1;
	nTot = Urls.nCatalogueFiles-1;
	sprintf(out_url, "%s%s", Urls.output_prefix, "halo_subhalo_evolution.dat");
	out_file = fopen(out_url, "w");

	DUMP_MSG("evolution", out_url);
	
		fprintf(out_file,"#");
		FILE_HEADER(out_file, "z         ", count);
		FILE_HEADER(out_file, "conc      ", count);
		FILE_HEADER(out_file, "lambda    ", count);
		FILE_HEADER(out_file, "shape     ", count);
		FILE_HEADER(out_file, "triax     ", count);
		FILE_HEADER(out_file, "avg_sub   ", count);
		FILE_HEADER(out_file, "sub_conc  ", count);
		FILE_HEADER(out_file, "sub_lambda", count);
		FILE_HEADER(out_file, "sub_shape ", count);
		FILE_HEADER(out_file, "sub_triax ", count);
		FILE_HEADER(out_file, "sub_costh ", count);
		FILE_HEADER(out_file, "sub_cosph ", count);
		FILE_HEADER(out_file, "sub_avg_v ", count);
		FILE_HEADER(out_file, "sub_avg_m ", count);
		fprintf(out_file,"\n");

			for(i=0; i<nTot; i++)
			{ 
				fprintf(out_file, "%lf", HaloProperties[i].z);   
				fprintf(out_file, "\t%lf", HaloProperties[i].c_0);   
				fprintf(out_file, "\t%lf", HaloProperties[i].l_0);   
				fprintf(out_file, "\t%lf", HaloProperties[i].s0);   
				fprintf(out_file, "\t%lf", HaloProperties[i].t0);   
				fprintf(out_file, "\t%lf", HaloProperties[i].avgSub);   
				fprintf(out_file, "\t%lf", SubHaloProperties[i].c_0);   	
				fprintf(out_file, "\t%lf", SubHaloProperties[i].l_0);   	
				fprintf(out_file, "\t%lf", SubHaloProperties[i].s0);   	
				fprintf(out_file, "\t%lf", SubHaloProperties[i].t0);   	
				fprintf(out_file, "\t%lf", SubHaloProperties[i].costh0);   	
				fprintf(out_file, "\t%lf", SubHaloProperties[i].cosphi0);   	
				fprintf(out_file, "\t%lf", SubHaloProperties[i].vel_0);   	
				fprintf(out_file, "\t%e", SubHaloProperties[i].avgMass);   	
				fprintf(out_file,"\n");
			}

	fclose(out_file);
}



void print_all_halo_properties_to_one_file()
{
	count = 1;
	nTot = HaloZ.n_bins-1;
	z = GrowthFac.z[Settings.use_cat];
	sprintf(out_url, "%sz%.2f%s", Urls.output_prefix, z, "_all_halo_statistical_properties.dat");
	out_file = fopen(out_url,"w");

	DUMP_MSG("halo properties", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "Mass  ", count);
		FILE_HEADER(out_file, "avg_c ", count);
		FILE_HEADER(out_file, "conc  ", count);
		FILE_HEADER(out_file, "P(c)  ", count);
		FILE_HEADER(out_file, "lambda", count);
		FILE_HEADER(out_file, "P(l)  ", count);
		FILE_HEADER(out_file, "triax ", count);
		FILE_HEADER(out_file, "P(t)  ", count);
		FILE_HEADER(out_file, "shape ", count);
		FILE_HEADER(out_file, "P(s)  ", count);
#ifdef GAS
		FILE_HEADER(out_file, "gas_temp", count);
		FILE_HEADER(out_file, "gas_frac", count);
#endif
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)	
			{

				fprintf(out_file, "%e", MassFunc.mass[i]);
				fprintf(out_file, "\t%lf", HaloZ.c_avg[i]);
				fprintf(out_file, "\t%lf", HaloZ.c[i]);
				fprintf(out_file, "\t%lf", HaloZ.p_c[i]);
				fprintf(out_file, "\t%lf", HaloZ.l[i]);
				fprintf(out_file, "\t%lf", HaloZ.p_l[i]);
				fprintf(out_file, "\t%lf", HaloZ.triax[i]);
				fprintf(out_file, "\t%lf", HaloZ.p_triax[i]);
				fprintf(out_file, "\t%lf", HaloZ.shape[i]);
				fprintf(out_file, "\t%lf", HaloZ.p_shape[i]);
#ifdef GAS
				fprintf(out_file, "\t%lf", HaloZ.gas_T[i]);
				fprintf(out_file, "\t%lf", HaloZ.gas_fraction[i]);

#endif
				fprintf(out_file, "\n");
			}

	fclose(out_file);
}



void print_axis_alignment()
{
	count = 1;
	nTot = HaloZ.r_bins-1;
	z = GrowthFac.z[Settings.use_cat];
	sprintf(out_url, "%sz%.2f%s", Urls.output_prefix, z, "_axis_alignement.dat");
	out_file = fopen(out_url, "w");

	DUMP_MSG("axis alignment", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "r      ", count);
		FILE_HEADER(out_file, "c(Tc)^2", count);
		FILE_HEADER(out_file, "c(Tp)^2", count);
		FILE_HEADER(out_file, "N      ", count);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++) 
			{
				fprintf(out_file, "%lf", HaloZ.R[i]);
				fprintf(out_file, "\t%lf", HaloZ.Th_c[i]/HaloZ.N_pairs[i]);
				fprintf(out_file, "\t%lf", HaloZ.Th_p[i]/HaloZ.N_pairs[i]);
				fprintf(out_file, "\t%d", HaloZ.N_pairs[i]);
				fprintf(out_file, "\n"); 
			}

	fclose(out_file);
}



void print_numerical_mass_function()
{
	count = 1;
	nTot = MassFunc.bins-1;
	sprintf(out_url, "%sz%.2f%s", Urls.output_prefix, z, "_numerical_mass_function.dat");
	out_file = fopen(out_url, "w");

	DUMP_MSG("numerical mass function", out_url);

		fprintf(out_file,"#");
		FILE_HEADER(out_file, "M    ", count);
		FILE_HEADER(out_file, "n    ", count);
		FILE_HEADER(out_file, "n_tot", count);
		FILE_HEADER(out_file, "n_err", count);
		FILE_HEADER(out_file, "M_step", count);
		FILE_HEADER(out_file, "dn   ", count);
		FILE_HEADER(out_file, "n_bin", count);
		FILE_HEADER(out_file, "dn_err", count);
		fprintf(out_file,"\n");

			for(i=0; i<nTot; i++) 
			{
				fprintf(out_file, "%e", MassFunc.mass[i]);
				fprintf(out_file, "\t%e", MassFunc.n[i]);
				fprintf(out_file, "\t%d\t", MassFunc.n_tot[i]);
				fprintf(out_file, "\t%e", MassFunc.err[i]);
				fprintf(out_file, "%e", MassFunc.mass_halfstep[i]);
				fprintf(out_file, "\t%e", MassFunc.dn[i]);
				fprintf(out_file, "\t%d\t", MassFunc.n_bin[i]);
				fprintf(out_file, "\t%e", MassFunc.err_dn[i]);
				fprintf(out_file,"\n");
			}

	fclose(out_file);
}



void print_nfw()
{
// TODO //FIXME
	int k=0, count=1, nTot=0;
	char *out = merge_strings(Urls.output_prefix, "nfw_test.dat");
	FILE *out_nfw = fopen(out, "w");

	nTot = NFW.bins;		

	DUMP_MSG("Navarro Frenk White profile", out);

		fprintf(out_nfw, "#");
		FILE_HEADER(out_nfw, "r", count);
		FILE_HEADER(out_nfw, "overd", count);
		FILE_HEADER(out_nfw, "rho_NFW", count);
		FILE_HEADER(out_nfw, "err", count);
		fprintf(out_nfw, "\n");

			for(k=0; k<nTot; k++)
			{
				fprintf(out_nfw, "%lf", NFW.radius[k]);
				fprintf(out_nfw, "\t%lf", NFW.overd[k]);
				fprintf(out_nfw, "\t%lf", NFW.profile[k]);
				fprintf(out_nfw, "\t%lf", NFW.err[k]);
				fprintf(out_nfw, "\n");
			}

	fclose(out_nfw);
}



void print_best_fit_results(){
// TODO //FIXME
	int i=0, count=1, nBins=HaloZ.n_bins;
	char* out_url=merge_strings(Urls.output_prefix,"all_halo_best_fit_distributions.dat");
	FILE* out_file=fopen(out_url,"w");
	DUMP_MSG("best fit", out_url);
	
		fprintf(out_file, "#Concentration distribution best fit values, c_0=%lf, c_sig  :%lf \n", 
			HaloZ.c_0,HaloZ.c_sig );
		fprintf(out_file, "#Lambda parameter distribution best fit values, l_0  :%lf, l_sig  :%lf \n", 
			HaloZ.l_0,HaloZ.l_sig );
		fprintf(out_file, "#Average shape parameter s_0  :%lf, average triaxiality t_0    :%lf \n", 
			HaloZ.s0,HaloZ.t0);

		fprintf(out_file, "#lambda_0  = %lf\n", SubHaloZ.l_0);
		fprintf(out_file, "#cos(th_0) = %lf\n", SubHaloZ.costh0);
		fprintf(out_file, "#cos(phi_0)= %lf\n", SubHaloZ.cosphi0);
		fprintf(out_file, "#t0 = %lf\n", SubHaloZ.t0);
		fprintf(out_file, "#s0 = %lf\n", SubHaloZ.s0);
		fprintf(out_file, "\n");

#ifdef GAS
			fprintf(out_file, "#Mass-Tx relation best fit values, ln_M0: %lf, alpha :%lf  \n", 
				HaloZ.ln_M0, HaloZ.T_alpha);
#endif
}
