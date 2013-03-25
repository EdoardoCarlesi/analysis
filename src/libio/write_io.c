#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../libhalo/halo.h"
#include "../libcosmo/cosmo.h"
#include "../general_def.h"

#ifdef WITH_MPI
#include "../libparallel/general.h"
#endif

#include "io.h"

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
	nTot = SubHaloProperties[HALO_INDEX].n_bins;
	z = GrowthFac.z[Settings.use_cat];
 	sprintf(out_url, "%sz%.3f%s", Urls.output_prefix, z, "_all_sub_statistics.dat");
	out_file = fopen(out_url,"w");

	DUMP_MSG("subhalo properties", out_url);

		FILE_HEADER(out_file, "Mass  ", count);
		FILE_HEADER(out_file, "Nbin  ", count);
		FILE_HEADER(out_file, "vel   ", count);
		FILE_HEADER(out_file, "conc  ", count);
		FILE_HEADER(out_file, "virial", count);
		FILE_HEADER(out_file, "lambda", count);
		FILE_HEADER(out_file, "shape ", count);
		FILE_HEADER(out_file, "triax ", count);
		FILE_HEADER(out_file, "gofnfw", count);
		FILE_HEADER(out_file, "g_nfw ", count);
		FILE_HEADER(out_file, "pg_nfw", count);
		FILE_HEADER(out_file, "c     ", count);
		FILE_HEADER(out_file, "P(c)  ", count);
		FILE_HEADER(out_file, "l     ", count);
		FILE_HEADER(out_file, "P(l)  ", count);
		FILE_HEADER(out_file, "t     ", count);
		FILE_HEADER(out_file, "P(t)  ", count);
		FILE_HEADER(out_file, "s     ", count);
		FILE_HEADER(out_file, "P(s)  ", count);
#ifdef GAS
		FILE_HEADER(out_file, "gas_t ", count);
		FILE_HEADER(out_file, "P(g_t)", count);
		FILE_HEADER(out_file, "gas_s ", count);
		FILE_HEADER(out_file, "P(g_s)", count);
		FILE_HEADER(out_file, "diff_t", count);
		FILE_HEADER(out_file, "P(d_t)", count);
		FILE_HEADER(out_file, "diff_s", count);
		FILE_HEADER(out_file, "P(d_s)", count);
		FILE_HEADER(out_file, "del_cm", count);
		FILE_HEADER(out_file, "P(cm) ", count);
		FILE_HEADER(out_file, "cos_th", count);
		FILE_HEADER(out_file, "P(cth)", count);
		FILE_HEADER(out_file, "g_spin", count);
		FILE_HEADER(out_file, "g_shap", count);
		FILE_HEADER(out_file, "g_tria", count);
		FILE_HEADER(out_file, "g_temp", count);
		FILE_HEADER(out_file, "g_frac", count);
		FILE_HEADER(out_file, "g_ekin", count);
		FILE_HEADER(out_file, "g_vir ", count);
		FILE_HEADER(out_file, "dmekin", count);
		FILE_HEADER(out_file, "dmvir ", count);
		FILE_HEADER(out_file, "g_coth", count);
		FILE_HEADER(out_file, "g_diff", count);
#endif

		FILE_HEADER(out_file, "Cos(th) ", count);
		FILE_HEADER(out_file, "P(c(th))", count);
		FILE_HEADER(out_file, "Cos(ph) ", count);
		FILE_HEADER(out_file, "P(c(ph))", count);
		FILE_HEADER(out_file, "sub(r)/R", count);
		FILE_HEADER(out_file, "n(r)    ", count);
		FILE_HEADER(out_file, "n(>r)   ", count);
		FILE_HEADER(out_file, "V_sub   ", count);
		FILE_HEADER(out_file, "P(V_sub)", count);
		FILE_HEADER(out_file, "n(>M)   ", count);
//		FILE_HEADER(out_file, "subset(r)", count);
//		FILE_HEADER(out_file, "subset(n(>r))", count);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)	
			{
				fprintf(out_file, "%e",    SubHaloProperties[HALO_INDEX].mass[i]);
				fprintf(out_file, "\t%d\t", SubHaloProperties[HALO_INDEX].n_entry[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].vel[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].conc[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].halo.virial[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].halo.lambda[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].halo.shape[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].halo.triax[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].fit_nfw.gof[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].p_fit_nfw.gof[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].p_fit_nfw.p_gof[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].c[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].p_c[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].halo.l[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].halo.p_l[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].halo.t[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].halo.p_t[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].halo.s[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].halo.p_s[i]);
#ifdef GAS
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas.t[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas.p_t[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas.s[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas.p_s[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].diff.t[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].diff.p_t[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].diff.s[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].diff.p_s[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].cm[i]);
				fprintf(out_file, "\t%e", SubHaloProperties[HALO_INDEX].p_cm[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas_dm_cth[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].p_gas_dm_cth[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas.lambda[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas.shape[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas.triax[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas_T[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas_fraction[i]);
				fprintf(out_file, "\t%e", SubHaloProperties[HALO_INDEX].gas.ekin[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas.virial[i]);
				fprintf(out_file, "\t%e", SubHaloProperties[HALO_INDEX].dm.ekin[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].dm.virial[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas_dm_costh[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].gas_diff_cm[i]);

#endif
	
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].costh[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].costh_count[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].cosphi[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].cosphi_count[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].r_sub[i]);
				fprintf(out_file, "\t%d" , SubHaloProperties[HALO_INDEX].n_r_sub[i]);
				fprintf(out_file, "\t%d" , SubHaloProperties[HALO_INDEX].cum_n_r_sub[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].vel_sub[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].p_vel_sub[i]);
				fprintf(out_file, "\t%d" , SubHaloProperties[HALO_INDEX].cum_n_sub[i]);
				//fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].r_sub_subset[i]);
				//fprintf(out_file, "\t%d" , SubHaloProperties[HALO_INDEX].cum_n_r_sub_subset[i]);
				fprintf(out_file, "\n");
			}

	fclose(out_file);	
}



void print_theoretical_mass_function()
{
	count = 1;
	nTot = ThMassFunc[MF_INDEX].bins-1;
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
				fprintf(out_file,"%e", ThMassFunc[MF_INDEX].mass[i]);
				fprintf(out_file,"\t%e", ThMassFunc[MF_INDEX].n[i]);
				fprintf(out_file,"\t%e", ThMassFunc[MF_INDEX].dn[i]);
				fprintf(out_file,"\n");
			}

	fclose(out_file);
}



void print_average_profiles()
{
	count = 1;
	nTot = BIN_PROFILE;
	z = GrowthFac.z[Settings.use_cat];
	sprintf(out_url, "%sz%.3f%s", Urls.output_prefix, z, "_avg_profiles.dat");
	out_file = fopen(out_url,"w");

	DUMP_MSG("average profiles", out_url);

		fprintf(out_file,"#");
		FILE_HEADER(out_file, "R/Rv    ", count);
		FILE_HEADER(out_file, "rho/rho0", count);
#ifdef GAS
		FILE_HEADER(out_file, "rho_gas ", count);
		FILE_HEADER(out_file, "Ix/Ix0  ", count);
		//FILE_HEADER(out_file, "R/Rv_fra", count);
		FILE_HEADER(out_file, "frac_gas", count);
#endif
		fprintf(out_file, "\n");
	
			for(i=0; i<nTot; i++) 
			{
				fprintf(out_file, "%f", HaloProperties[HALO_INDEX].nfw.x[i]); 
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].nfw.y[i]); 
#ifdef GAS
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].rho_gas.y[i]); 
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].i_x.y[i]); 
				//fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].f_gas.x[i]); 
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].f_gas.y[i]); 
				fprintf(out_file, "\n");
#endif
			}

	fclose(out_file);

}



void print_correlation_function()
{
	count = 1;
	nTot = Xi.npts;
	z = GrowthFac.z[Settings.use_cat];
	sprintf(out_url, "%sz%.3f%s", Urls.output_prefix, z, "_correlation_function.dat");
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
				fprintf(out_file, "\t%lf", HaloProperties[i].l_0);   
				fprintf(out_file, "\t%lf", HaloProperties[i].halo.s0);   
				fprintf(out_file, "\t%lf", HaloProperties[i].halo.t0);   
				fprintf(out_file, "\t%lf", HaloProperties[i].avgSub);   
				fprintf(out_file, "\t%lf", SubHaloProperties[i].l_0);   	
				fprintf(out_file, "\t%lf", SubHaloProperties[i].halo.s0);   	
				fprintf(out_file, "\t%lf", SubHaloProperties[i].halo.t0);   	
				fprintf(out_file, "\t%lf", SubHaloProperties[i].costh0);   	
				fprintf(out_file, "\t%lf", SubHaloProperties[i].cosphi0);   	
				fprintf(out_file, "\t%lf", SubHaloProperties[i].vel_0);   	
				fprintf(out_file, "\t%e", SubHaloProperties[i].avgMass);   	
				fprintf(out_file,"\n");
			}

	fclose(out_file);
}



void print_all_haloes()
{
	count = 1;
	nTot = Settings.n_haloes;
	z = GrowthFac.z[Settings.use_cat];
	sprintf(out_url, "%sz%.3f%s", Urls.output_prefix, z, "_tot_haloes.dat");
	out_file = fopen(out_url,"w");

	DUMP_MSG("halo", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "Mass  ", count);
		FILE_HEADER(out_file, "conc  ", count);
		FILE_HEADER(out_file, "virial", count);
		FILE_HEADER(out_file, "lambda", count);
		FILE_HEADER(out_file, "shape ", count);
		FILE_HEADER(out_file, "triax ", count);
#ifdef GAS
		FILE_HEADER(out_file, "gas_T ", count);
		FILE_HEADER(out_file, "gas_fr", count);
#endif
		fprintf(out_file, "\n");

		for(i=0; i<nTot; i++)	
		{

			if(halo_condition(i)==1)
			{
				fprintf(out_file, "%e", Haloes[i].Mvir); 
				fprintf(out_file, "\t%f", Haloes[i].c);
				fprintf(out_file, "\t%f", Haloes[i].abs_th_vir);
				fprintf(out_file, "\t%f", Haloes[i].lambda);
				fprintf(out_file, "\t%f", Haloes[i].shape);
				fprintf(out_file, "\t%f", Haloes[i].triax);
#ifdef GAS
				fprintf(out_file, "\t%f", Haloes[i].gas_only.T_mw);
				fprintf(out_file, "\t%f", Haloes[i].gas_only.b_fraction);
				fprintf(out_file, "\n");
#endif
			}

		}

	fclose(out_file);
}




void print_all_halo_properties_to_one_file()
{
	count = 1;
	nTot = HaloProperties[HALO_INDEX].n_bins;
	z = GrowthFac.z[Settings.use_cat];
	sprintf(out_url, "%sz%.3f%s", Urls.output_prefix, z, "_all_halo_statistical_properties.dat");
	out_file = fopen(out_url,"w");

	DUMP_MSG("halo properties", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "Mass  ", count);
		FILE_HEADER(out_file, "Nbin  ", count);
		FILE_HEADER(out_file, "vel   ", count);
		FILE_HEADER(out_file, "conc  ", count);
		FILE_HEADER(out_file, "virial", count);
		FILE_HEADER(out_file, "lambda", count);
		FILE_HEADER(out_file, "shape ", count);
		FILE_HEADER(out_file, "triax ", count);
		FILE_HEADER(out_file, "gofnfw", count);
		FILE_HEADER(out_file, "g_nfw ", count);
		FILE_HEADER(out_file, "pg_nfw", count);
		FILE_HEADER(out_file, "T     ", count);
		FILE_HEADER(out_file, "n(>T) ", count);
		FILE_HEADER(out_file, "c     ", count);
		FILE_HEADER(out_file, "P(c)  ", count);
		FILE_HEADER(out_file, "l     ", count);
		FILE_HEADER(out_file, "P(l)  ", count);
		FILE_HEADER(out_file, "t     ", count);
		FILE_HEADER(out_file, "P(t)  ", count);
		FILE_HEADER(out_file, "s     ", count);
		FILE_HEADER(out_file, "P(s)  ", count);
#ifdef GAS
		FILE_HEADER(out_file, "gas_t ", count);
		FILE_HEADER(out_file, "P(g_t)", count);
		FILE_HEADER(out_file, "gas_s ", count);
		FILE_HEADER(out_file, "P(g_s)", count);
		FILE_HEADER(out_file, "diff_t", count);
		FILE_HEADER(out_file, "P(d_t)", count);
		FILE_HEADER(out_file, "diff_s", count);
		FILE_HEADER(out_file, "P(d_s)", count);
		FILE_HEADER(out_file, "del_cm", count);
		FILE_HEADER(out_file, "P(cm) ", count);
		FILE_HEADER(out_file, "cos_th", count);
		FILE_HEADER(out_file, "P(cth)", count);
		FILE_HEADER(out_file, "g_spin", count);
		FILE_HEADER(out_file, "g_shap", count);
		FILE_HEADER(out_file, "g_tria", count);
		FILE_HEADER(out_file, "g_temp", count);
		FILE_HEADER(out_file, "g_beta", count);
		FILE_HEADER(out_file, "g_frac", count);
		FILE_HEADER(out_file, "g_ekin", count);
		FILE_HEADER(out_file, "g_vir ", count);
		FILE_HEADER(out_file, "dmekin", count);
		FILE_HEADER(out_file, "dmvir ", count);
		FILE_HEADER(out_file, "g_coth", count);
		FILE_HEADER(out_file, "g_diff", count);
#endif
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)	
			{

				fprintf(out_file, "%e",    HaloProperties[HALO_INDEX].mass[i]);
				fprintf(out_file, "\t%d\t",  HaloProperties[HALO_INDEX].n_entry[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].vel[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].conc[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.virial[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.lambda[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.shape[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.triax[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].fit_nfw.gof[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_fit_nfw.gof[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_fit_nfw.p_gof[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].T[i]);
				fprintf(out_file, "\t%e", HaloProperties[HALO_INDEX].n_T[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].c[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_c[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.l[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.p_l[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.p_t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.s[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.p_s[i]);
#ifdef GAS
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.p_t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.s[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.p_s[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].diff.t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].diff.p_t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].diff.s[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].diff.p_s[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].cm[i]);
				fprintf(out_file, "\t%e",  HaloProperties[HALO_INDEX].p_cm[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas_dm_cth[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_gas_dm_cth[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.lambda[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.shape[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.triax[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas_T[i]);
				fprintf(out_file, "\t%f",  HaloProperties[HALO_INDEX].gas.beta[i]);
				fprintf(out_file, "\t%f",  HaloProperties[HALO_INDEX].gas_fraction[i]);
				fprintf(out_file, "\t%e",  HaloProperties[HALO_INDEX].gas.ekin[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.virial[i]);
				fprintf(out_file, "\t%e",  HaloProperties[HALO_INDEX].dm.ekin[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].dm.virial[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas_dm_costh[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas_diff_cm[i]);

#endif
				fprintf(out_file, "\n");
			}

	fclose(out_file);
}



void print_axis_alignment()
{
	count = 1;
	nTot = HaloProperties[HALO_INDEX].r_bins-1;
	z = GrowthFac.z[Settings.use_cat];
	sprintf(out_url, "%sz%.3f%s", Urls.output_prefix, z, "_axis_alignement.dat");
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
				fprintf(out_file, "%lf", HaloProperties[HALO_INDEX].R[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].Th_c[i]
					/ HaloProperties[HALO_INDEX].N_pairs[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].Th_p[i]
					/ HaloProperties[HALO_INDEX].N_pairs[i]);
				fprintf(out_file, "\t%d", HaloProperties[HALO_INDEX].N_pairs[i]);
				fprintf(out_file, "\n"); 
			}

	fclose(out_file);
}



void print_numerical_mass_function()
{
	count = 1;
	nTot = MassFunc[MF_INDEX].bins-1;
	sprintf(out_url, "%sz%.3f%s", Urls.output_prefix, z, "_numerical_mass_function.dat");
	out_file = fopen(out_url, "w");

	DUMP_MSG("numerical mass function", out_url);

		fprintf(out_file,"#");
		FILE_HEADER(out_file, "M    ", count);
		FILE_HEADER(out_file, "M_step", count);
		FILE_HEADER(out_file, "n    ", count);
		FILE_HEADER(out_file, "n_tot", count);
		FILE_HEADER(out_file, "n_err", count);
		FILE_HEADER(out_file, "dn   ", count);
		FILE_HEADER(out_file, "n_bin", count);
		FILE_HEADER(out_file, "dn_err", count);
		fprintf(out_file,"\n");

			for(i=0; i<nTot; i++) 
			{
				fprintf(out_file, "%e", MassFunc[MF_INDEX].mass[i]);
				fprintf(out_file, "\t%e", MassFunc[MF_INDEX].mass_halfstep[i]);
				fprintf(out_file, "\t%e", MassFunc[MF_INDEX].n[i]);
				fprintf(out_file, "\t%d\t", MassFunc[MF_INDEX].n_tot[i]);
				fprintf(out_file, "\t%e", MassFunc[MF_INDEX].err[i]);
				fprintf(out_file, "\t%e", MassFunc[MF_INDEX].dn[i]);
				fprintf(out_file, "\t%d\t", MassFunc[MF_INDEX].n_bin[i]);
				fprintf(out_file, "\t%e", MassFunc[MF_INDEX].err_dn[i]);
				fprintf(out_file,"\n");
			}

	fclose(out_file);
}



void print_halo_profile(int m)
{
	sprintf(out_url, "%s.%04d.%03d_%s", Urls.output_prefix, ThisTask, i, "profiles.dat");
	out_file = fopen(out_url, "w");

	struct halo * HALO;
	nTot = HALO[m].n_bins - HALO[m].neg_r_bins;

#ifdef WITH_MPI
	HALO = pHaloes[ThisTask];
#else 
	HALO = Haloes;
#endif

	if(ThisTask==0)

	DUMP_MSG("Dumping halo density profile", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "r      ", count);
		FILE_HEADER(out_file, "rho_DM ", count);
#ifdef GAS
		FILE_HEADER(out_file, "rho_gas", count);
		FILE_HEADER(out_file, "gas_fra", count);
		FILE_HEADER(out_file, "I_X    ", count);
#endif
		fprintf(out_file, "\n");
		fprintf(out_file, "#");
		fprintf(out_file, "M=%e ", HALO[m].Mvir);
		fprintf(out_file, "X=%f ", HALO[m].X[0]);
		fprintf(out_file, "Y=%f ", HALO[m].X[1]);
		fprintf(out_file, "Z=%f ", HALO[m].X[2]);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)
			{
				fprintf(out_file, "%lf",   HALO[m].radius[i]);
				fprintf(out_file, "\t%lf", HALO[m].rho[i]);
#ifdef GAS
				fprintf(out_file, "\t%lf", HALO[m].gas_only.frac[i]);
				fprintf(out_file, "\t%lf", HALO[m].gas_only.rho[i]);
				fprintf(out_file, "\t%lf", HALO[m].gas_only.i_x[i]);
#endif
				fprintf(out_file, "\n");
			}

	fclose(out_file);
}


void print_halo_best_fit_results()
{
	sprintf(out_url, "%sz%.3f%s", Urls.output_prefix, z, "_halo_best_fit.dat");
	DUMP_MSG("best fit", out_url);
	out_file=fopen(out_url,"w");
	
		fprintf(out_file, "#Mass-concentration best fit values, c_0=%f, c_beta  :%f \n", 
			HaloProperties[HALO_INDEX].c_0,HaloProperties[HALO_INDEX].c_beta );
		fprintf(out_file, "#Mass-velocity best fit values, vel_0=%f, vel_beta  :%f \n", 
			HaloProperties[HALO_INDEX].vel_0,HaloProperties[HALO_INDEX].vel_beta );
		fprintf(out_file, "#Lambda parameter distribution best fit values, l_0  :%f, l_sig  :%f \n", 
			HaloProperties[HALO_INDEX].l_0,HaloProperties[HALO_INDEX].l_sig );
		fprintf(out_file, "#Average shape parameter s_0  :%lf, average triaxiality t_0    :%lf \n", 
			HaloProperties[HALO_INDEX].halo.s0, HaloProperties[HALO_INDEX].halo.t0);
		fprintf(out_file, "\n");

#ifdef GAS
			fprintf(out_file, "#Mass-T best fit values, T0: %lf, alpha :%lf  \n", 
				HaloProperties[HALO_INDEX].T0, HaloProperties[HALO_INDEX].alpha);
			fprintf(out_file, "#Average beta:%lf\n", HaloProperties[HALO_INDEX].beta0);
#endif
	fclose(out_file);
}
