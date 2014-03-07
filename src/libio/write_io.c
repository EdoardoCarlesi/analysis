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



void print_subhalo_only_properties()
{
	count = 1;
	nTot = (int) (F_SUB * HaloProperties[HALO_INDEX].n_bins) - 1;
 	sprintf(out_url, "%s%s", Urls.output_prefix, "_sub_only_statistics.dat");
	out_file = fopen(out_url,"w");

	DUMP_MSG("subhalo properties", out_url);

		FILE_HEADER(out_file, "r/R     ", count);
		FILE_HEADER(out_file, "n(r)    ", count);
		FILE_HEADER(out_file, "n(>r)   ", count);
		FILE_HEADER(out_file, "V_sub(r)", count);
		FILE_HEADER(out_file, "V_sub   ", count);
		FILE_HEADER(out_file, "P(V_sub)", count);
		FILE_HEADER(out_file, "Cos(th) ", count);
		FILE_HEADER(out_file, "P(c(th))", count);
		FILE_HEADER(out_file, "Cos(ph) ", count);
		FILE_HEADER(out_file, "P(c(ph))", count);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)	
			{
				fprintf(out_file, "%f", HaloProperties[HALO_INDEX].r_sub[i]);
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].n_r_sub[i]);
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].cum_n_r_sub[i]);
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].vel_sub_r[i]);
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].vel_sub[i]);
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].p_vel_sub[i]);
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].costh[i]);
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].costh_count[i]);
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].cosphi[i]);
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].cosphi_count[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].r_sub_subset[i]);
				//fprintf(out_file, "\t%d" , HaloProperties[HALO_INDEX].cum_n_r_sub_subset[i]);
				fprintf(out_file, "\n");
			}

	fclose(out_file);	
}



void print_theoretical_mass_function()
{
	count = 1;
	nTot = ThMassFunc[MF_INDEX].bins-1;
	sprintf(out_url, "%s%s", Urls.output_prefix, "_theoretical_mass_function.dat");
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
	sprintf(out_url, "%s%s", Urls.output_prefix, "_avg_profiles.dat");
	out_file = fopen(out_url,"w");

	double T0 = HaloProperties[HALO_INDEX].temp.y[BIN_PROFILE-1];
	double rho0 = HaloProperties[HALO_INDEX].rho_gas.y[BIN_PROFILE-1];
			//get_interpolated_value(HaloProperties[HALO_INDEX].temp.x, 
			//HaloProperties[HALO_INDEX].temp.y, 0.25, BIN_PROFILE);

	DUMP_MSG("average profiles", out_url);

		fprintf(out_file,"#");
		FILE_HEADER(out_file, "R/Rv    ", count);
		FILE_HEADER(out_file, "rho/rho0", count);
		FILE_HEADER(out_file, "N_bin   ", count);
#ifdef GAS
		FILE_HEADER(out_file, "rho_gas ", count);
		FILE_HEADER(out_file, "Ix/Ix0  ", count);
		FILE_HEADER(out_file, "frac_gas", count);
		FILE_HEADER(out_file, "T/T0    ", count);
	//	FILE_HEADER(out_file, "hydro_M ", count);
	//	FILE_HEADER(out_file, "pressure", count);
#endif
		fprintf(out_file, "\n");

			fprintf(out_file,"#Average Mvir = %e SM/h, Median Mvir = %e SM/h\n",
				HaloProperties[HALO_INDEX].avgMvir, 
				HaloProperties[HALO_INDEX].medMvir);

			fprintf(out_file,"#Average Msub = %e SM/h, Median Msub = %e SM/h\n",
				HaloProperties[HALO_INDEX].avgMsub, 
				HaloProperties[HALO_INDEX].medMsub);

			fprintf(out_file,"#Average Nsub = %f, Median Nsub = %f\n",
				HaloProperties[HALO_INDEX].avgNsub, 
				HaloProperties[HALO_INDEX].medNsub);

			fprintf(out_file,"#Average Rvir = %e Mpc/h, Median Rvir = %e Mpc/h\n",
				HaloProperties[HALO_INDEX].avgRvir, 
				HaloProperties[HALO_INDEX].medRvir);
	
			for(i=2; i<nTot; i++) 
			{
				fprintf(out_file, "%lf",  HaloProperties[HALO_INDEX].nfw.x[i]); 
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].nfw.y[i]); 
				fprintf(out_file, "\t%d", HaloProperties[HALO_INDEX].nfw.n[i]); 
#ifdef GAS
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].rho_gas.y[i]/rho0); 
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].i_x.y[i]); 
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].f_gas.y[i]); 
				fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].temp.y[i]); 
			//	fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].hydro_m.y[i]); 
			//	fprintf(out_file, "\t%f", HaloProperties[HALO_INDEX].pressure.y[i]); 
#endif
				fprintf(out_file, "\n");
			}

	fclose(out_file);

}



void print_correlation_function()
{
	count = 1;
	nTot = Xi.npts;
	sprintf(out_url, "%s%s", Urls.output_prefix, "_correlation_function.dat");
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

/*
			for(i=0; i<nTot; i++)
			{ 
				fprintf(out_file, "%lf  ", HaloProperties[i].z);   
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
				fprintf(out_file, "\t%e ", SubHaloProperties[i].avgMass);   	
				fprintf(out_file,"\n");
			}
*/

	fclose(out_file);
}



void print_all_haloes()
{
	count = 1;
	nTot = Settings.n_haloes;
	sprintf(out_url, "%s%s", Urls.output_prefix, "_tot_haloes.dat");
	out_file = fopen(out_url,"w");

	DUMP_MSG("halo", out_url);

		fprintf(out_file, "#");
#ifndef WEB_ONLY
		FILE_HEADER(out_file, "Mass  ", count);
		FILE_HEADER(out_file, "Rvir  ", count);
		FILE_HEADER(out_file, "Msub  ", count);
#endif
		FILE_HEADER(out_file, "n_part", count);
		FILE_HEADER(out_file, "X     ", count);
		FILE_HEADER(out_file, "Y     ", count);
		FILE_HEADER(out_file, "Z     ", count);
#ifndef NO_PROFILES
		FILE_HEADER(out_file, "conc  ", count);
#endif

#ifdef WEB_ONLY
		FILE_HEADER(out_file, "Web   ", count);
		FILE_HEADER(out_file, "CatNum", count);
		FILE_HEADER(out_file, "Line  ", count);

#else
		FILE_HEADER(out_file, "lambda", count);
		FILE_HEADER(out_file, "shape ", count);
		FILE_HEADER(out_file, "triax ", count);
#ifdef GAS
		FILE_HEADER(out_file, "gas_fr", count);
	      	FILE_HEADER(out_file, "gasTmw", count);
		FILE_HEADER(out_file, "gasTew", count);
		FILE_HEADER(out_file, "gasTsl", count);
#endif
#endif
		FILE_HEADER(out_file, "vir_dm", count);
		fprintf(out_file, "\n");

		for(i=0; i<nTot; i++)	
		{
 
			if(halo_condition(i)==1)
			{
#ifndef WEB_ONLY
				fprintf(out_file, "%e", Haloes[i].Mvir); 
				fprintf(out_file, "\t%lf", Haloes[i].Rvir); 
				fprintf(out_file, "\t%lf", Haloes[i].Msub/Haloes[i].Mvir); 

#ifndef NO_PROFILES
				fprintf(out_file, "\t%f", Haloes[i].fit_nfw.c);
#endif
				fprintf(out_file, "%e", Haloes[i].Mvir); 
				fprintf(out_file, "\t%lf", Haloes[i].Rvir); 
				fprintf(out_file, "\t%lf", Haloes[i].Msub/Haloes[i].Mvir); 
#endif
				fprintf(out_file, "%d\t", Haloes[i].n_part); 
				fprintf(out_file, "\t%f", Haloes[i].X[0]);
				fprintf(out_file, "\t%f", Haloes[i].X[1]);
				fprintf(out_file, "\t%f", Haloes[i].X[2]);
#ifndef NO_WEB
				fprintf(out_file, "\t%d\t", Haloes[i].c_web);
#endif

#ifdef WEB_ONLY
				fprintf(out_file, "\t%d\t", Haloes[i].cat_numb);
				fprintf(out_file, "\t%d\t", Haloes[i].cat_line);
#else
				fprintf(out_file, "\t%f", Haloes[i].fit_nfw.c);
				fprintf(out_file, "\t%f", Haloes[i].lambda);
				fprintf(out_file, "\t%f", Haloes[i].shape);
				fprintf(out_file, "\t%f", Haloes[i].triax);
#ifdef GAS
				fprintf(out_file, "\t%f", -2.*Haloes[i].dm.Ekin/Haloes[i].dm.Epot);
				fprintf(out_file, "\t%f", Haloes[i].gas_only.b_fraction);
				fprintf(out_file, "\t%f", Haloes[i].gas_only.T_mw);
				fprintf(out_file, "\t%f", Haloes[i].gas_only.T_ew);
				fprintf(out_file, "\t%f", Haloes[i].gas_only.T_sl);
#else
				fprintf(out_file, "\t%f", -2.*Haloes[i].Ekin/Haloes[i].Epot);
#endif
#endif
				fprintf(out_file, "\n");
			}

		}

	fclose(out_file);
}




void print_all_halo_properties_to_one_file()
{
	count = 1;
	nTot = HaloProperties[HALO_INDEX].n_bins;
	sprintf(out_url, "%s%s", Urls.output_prefix, "_all_halo_statistical_properties.dat");
	out_file = fopen(out_url,"w");

	DUMP_MSG("halo properties", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "Mass  ", count);
		FILE_HEADER(out_file, "Nbin  ", count);
		FILE_HEADER(out_file, "vel   ", count);
		FILE_HEADER(out_file, "conc  ", count);
		//FILE_HEADER(out_file, "virial", count);
		FILE_HEADER(out_file, "lambda", count);
		//FILE_HEADER(out_file, "shape ", count);
		//FILE_HEADER(out_file, "triax ", count);
		//FILE_HEADER(out_file, "gofnfw", count);
#ifndef NO_PROFILES
		FILE_HEADER(out_file, "g_nfw ", count);
		FILE_HEADER(out_file, "pg_nfw", count);
#endif
		FILE_HEADER(out_file, "c     ", count);
		FILE_HEADER(out_file, "P(c)  ", count);
		FILE_HEADER(out_file, "l     ", count);
		FILE_HEADER(out_file, "P(l)  ", count);
		FILE_HEADER(out_file, "t     ", count);
		FILE_HEADER(out_file, "P(t)  ", count);
		FILE_HEADER(out_file, "s     ", count);
		FILE_HEADER(out_file, "P(s)  ", count);
#ifdef GAS
		//FILE_HEADER(out_file, "gas_t ", count);
		//FILE_HEADER(out_file, "P(g_t)", count);
		//FILE_HEADER(out_file, "gas_s ", count);
		//FILE_HEADER(out_file, "P(g_s)", count);
		//FILE_HEADER(out_file, "diff_t", count);
		//FILE_HEADER(out_file, "P(d_t)", count);
		//FILE_HEADER(out_file, "diff_s", count);
		//FILE_HEADER(out_file, "P(d_s)", count);
		//FILE_HEADER(out_file, "del_cm", count);
		//FILE_HEADER(out_file, "P(cm) ", count);
		FILE_HEADER(out_file, "cos_th", count);
		FILE_HEADER(out_file, "P(cth)", count);
		//FILE_HEADER(out_file, "g_spin", count);
		//FILE_HEADER(out_file, "g_shap", count);
		//FILE_HEADER(out_file, "g_tria", count);
		FILE_HEADER(out_file, "g_temp", count);
		//FILE_HEADER(out_file, "g_beta", count);
		FILE_HEADER(out_file, "g_frac", count);
		//FILE_HEADER(out_file, "g_ekin", count);
		//FILE_HEADER(out_file, "g_vir ", count);
		//FILE_HEADER(out_file, "dmekin", count);
		FILE_HEADER(out_file, "dmvir ", count);
		//FILE_HEADER(out_file, "g_coth", count);
		//FILE_HEADER(out_file, "g_diff", count);
#ifndef NO_PROFILES
		FILE_HEADER(out_file, "mhydro", count);
		FILE_HEADER(out_file, "n_mhyd", count);
		FILE_HEADER(out_file, "s_mhyd", count);
		FILE_HEADER(out_file, "t_mhyd", count);
		FILE_HEADER(out_file, "gamma ", count);
		FILE_HEADER(out_file, "n_gamm", count);
		FILE_HEADER(out_file, "s_gamm", count);
		FILE_HEADER(out_file, "t_gamm", count);
#endif
		//FILE_HEADER(out_file, "n_s_mh", count);
		//FILE_HEADER(out_file, "n_t_mh", count);
		//FILE_HEADER(out_file, "sub_mh", count);
		//FILE_HEADER(out_file, "n_subm", count);
		//FILE_HEADER(out_file, "n_t_ga", count);
		//FILE_HEADER(out_file, "n_s_ga", count);
		//FILE_HEADER(out_file, "sub_ga", count);
		//FILE_HEADER(out_file, "n_subg", count);
#endif
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)	
			{

				fprintf(out_file, "%e",    HaloProperties[HALO_INDEX].mass[i]);
				fprintf(out_file, "\t%d\t",HaloProperties[HALO_INDEX].n_entry[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].vel[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].conc[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.virial[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.lambda[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.shape[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.triax[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].fit_nfw.gof[i]);
#ifndef NO_PROFILES
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_fit_nfw.gof[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_fit_nfw.p_gof[i]);
#endif
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].c[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_c[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.l[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.p_l[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.p_t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.s[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].halo.p_s[i]);
#ifdef GAS
				/*
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.p_t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.s[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.p_s[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].diff.t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].diff.p_t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].diff.s[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].diff.p_s[i]);
				fprintf(out_file, "\t%e", HaloProperties[HALO_INDEX].cm[i]);
				fprintf(out_file, "\t%e",  HaloProperties[HALO_INDEX].p_cm[i]);
				*/
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas_dm_cth[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_gas_dm_cth[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.lambda[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.shape[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.triax[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas_T[i]);
				//fprintf(out_file, "\t%f",  HaloProperties[HALO_INDEX].gas.beta[i]);
				fprintf(out_file, "\t%f",  HaloProperties[HALO_INDEX].gas_fraction[i]);
				//fprintf(out_file, "\t%e",  HaloProperties[HALO_INDEX].gas.ekin[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas.virial[i]);
				//fprintf(out_file, "\t%e",  HaloProperties[HALO_INDEX].dm.ekin[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].dm.virial[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas_dm_costh[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas_diff_cm[i]);
#ifndef NO_PROFILES
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].dM_hydro[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].dM_hydro_bin[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].shape_dM[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].triax_dM[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gamma[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gamma_bin[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].shape_g[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].triax_g[i]);
#endif
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].shape_dM_hydro[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].triax_dM_hydro[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].sub_dM[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].sub_dM_hydro[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].shape_gamma[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].triax_gamma[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].sub_g[i]);
				//fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].sub_gamma[i]);
#endif
				fprintf(out_file, "\n");
			}

	fclose(out_file);
}



void print_axis_alignment()
{
	count = 1;
	nTot = HaloProperties[HALO_INDEX].r_bins-1;
	sprintf(out_url, "%s%s", Urls.output_prefix, "_axis_alignement.dat");
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
	sprintf(out_url, "%s%s", Urls.output_prefix, "_numerical_mass_function.dat");
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
		FILE_HEADER(out_file, "V    ", count);
		FILE_HEADER(out_file, "n    ", count);
		FILE_HEADER(out_file, "dn   ", count);
		FILE_HEADER(out_file, "n_bin", count);
		FILE_HEADER(out_file, "n_tot", count);
#ifdef GAS
		FILE_HEADER(out_file, "Mgas ", count);
		FILE_HEADER(out_file, "n    ", count);
		FILE_HEADER(out_file, "n_tot", count);
		FILE_HEADER(out_file, "Nogas", count);
		FILE_HEADER(out_file, "n    ", count);
		FILE_HEADER(out_file, "n_tot", count);
		FILE_HEADER(out_file, "Dark ", count);
		FILE_HEADER(out_file, "n    ", count);
		FILE_HEADER(out_file, "n_tot", count);
		FILE_HEADER(out_file, "Temp ", count);
		FILE_HEADER(out_file, "n(>T)", count);
		FILE_HEADER(out_file, "N(>T)", count);
#endif
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
				fprintf(out_file, "\t%f", VelFunc[MF_INDEX].mass[i]);
				fprintf(out_file, "\t%e", VelFunc[MF_INDEX].n[i]);
				fprintf(out_file, "\t%e", VelFunc[MF_INDEX].dn[i]);
				fprintf(out_file, "\t%d\t", VelFunc[MF_INDEX].n_bin[i]);
				fprintf(out_file, "\t%d\t", VelFunc[MF_INDEX].n_tot[i]);
#ifdef GAS
				fprintf(out_file, "\t%e", GasFunc[MF_INDEX].mass[i]);
				fprintf(out_file, "\t%e", GasFunc[MF_INDEX].n[i]);
				fprintf(out_file, "\t%d\t", GasFunc[MF_INDEX].n_tot[i]);
				fprintf(out_file, "\t%e", NoGasFunc[MF_INDEX].mass[i]);
				fprintf(out_file, "\t%e", NoGasFunc[MF_INDEX].n[i]);
				fprintf(out_file, "\t%d\t", NoGasFunc[MF_INDEX].n_tot[i]);
				fprintf(out_file, "\t%e", DarkFunc[MF_INDEX].mass[i]);
				fprintf(out_file, "\t%e", DarkFunc[MF_INDEX].n[i]);
				fprintf(out_file, "\t%d\t", DarkFunc[MF_INDEX].n_tot[i]);
				fprintf(out_file, "\t%e", TempFunc[MF_INDEX].mass[i]);
				fprintf(out_file, "\t%e", TempFunc[MF_INDEX].n[i]);
				fprintf(out_file, "\t%d\t", TempFunc[MF_INDEX].n_tot[i]);
#endif
				fprintf(out_file,"\n");
			}

	fclose(out_file);
}


#ifndef NO_PROFILES
void print_halo_profile(int m)
{
	count = 1;
	sprintf(out_url, "%shalo.%04d.%03d_%s", Urls.output_prefix, ThisTask, m, "profile.dat");
	out_file = fopen(out_url, "w");
	struct halo * HALO;

#ifdef WITH_MPI
	HALO = pHaloes[ThisTask];
#else 
	HALO = Haloes;
#endif

	int skip = HALO[m].neg_r_bins;
	nTot = HALO[m].n_bins - HALO[m].neg_r_bins;

	DUMP_MSG("halo density profile", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "r      ", count);
		FILE_HEADER(out_file, "r/Rv   ", count);
		FILE_HEADER(out_file, "rho_DM ", count);
#ifdef GAS
		FILE_HEADER(out_file, "gas_fra", count);
		//FILE_HEADER(out_file, "rho_gas", count);
		//FILE_HEADER(out_file, "I_X    ", count);
		FILE_HEADER(out_file, "M(r)   ", count);
		FILE_HEADER(out_file, "hydro_M", count);
		FILE_HEADER(out_file, "T      ", count);
#endif
		fprintf(out_file, "\n");
		fprintf(out_file, "#");
		fprintf(out_file, "M=%e\t", HALO[m].Mvir);
		fprintf(out_file, "R=%f\t", HALO[m].Rvir);
		fprintf(out_file, "X=%f\t", HALO[m].X[0]);
		fprintf(out_file, "Y=%f\t", HALO[m].X[1]);
		fprintf(out_file, "Z=%f\t", HALO[m].X[2]);
		fprintf(out_file, "\n");
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)
			{
				fprintf(out_file, "%f",   HALO[m].radius[i+skip]);
				fprintf(out_file, "\t%f",  HALO[m].radius[i+skip]/HALO[m].Rvir);
				fprintf(out_file, "\t%lf", HALO[m].rho[i+skip]);
#ifdef GAS
				fprintf(out_file, "\t%lf", HALO[m].gas_only.frac[i+skip]);
				//fprintf(out_file, "\t%lf", HALO[m].gas_only.rho[i+skip]);
				//fprintf(out_file, "\t%lf", HALO[m].gas_only.i_x[i+skip]);
				fprintf(out_file, "\t%e", HALO[m].mass_r[i+skip]);
				fprintf(out_file, "\t%e", HALO[m].gas_only.hydro_m[i+skip]);
				fprintf(out_file, "\t%e", HALO[m].gas_only.T[i+skip]);
#endif
				fprintf(out_file, "\n");
			}

	fclose(out_file);
}
#else
void print_halo_profile(int m)
{
	// DO nothing
}
#endif
			/* Averaged distribution of all subhaloes */
void print_all_sub_per_host(int bins, int NsubTot, int NsubTh, int *n, double *r, double *r_n, double *r_n_c, 
	double *r_v, double *r_m, double *r_m_c, double *costh, int *p_costh, double *cosphi, int *p_cosphi)
{
	count = 1;
	sprintf(out_url, "%shalo.bin_%s", Urls.output_prefix, "distribution.dat");
	out_file = fopen(out_url, "w");
	nTot = bins-1;

	DUMP_MSG("satellite distribution in host halo", out_url);
/*
		fprintf(out_file, "#");
		FILE_HEADER(out_file, "r/Rv   ", count);
		FILE_HEADER(out_file, "n_bin  ", count);
		FILE_HEADER(out_file, "N_sub  ", count);
		FILE_HEADER(out_file, "N_sub_c", count);
		FILE_HEADER(out_file, "M_sub  ", count);
		FILE_HEADER(out_file, "M_sub_c", count);
		FILE_HEADER(out_file, "V_sub  ", count);
		fprintf(out_file, "\n");
		fprintf(out_file, "#Nsub=%d\tNsubTh=%d", NsubTot, NsubTh);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)
			{
				fprintf(out_file, "%lf", r[i]);
				fprintf(out_file, "\t%d\t",n[i]); 
				fprintf(out_file, "\t%lf", r_n[i]); 
				fprintf(out_file, "\t%lf", r_n_c[i]); 
				fprintf(out_file, "\t%lf", r_m[i]);
				fprintf(out_file, "\t%lf", r_m_c[i]);
				fprintf(out_file, "\t%lf", r_v[i]);
				fprintf(out_file, "\n");
			}
	fclose(out_file);
*/
	count = 1;
	sprintf(out_url, "%shalo.bin_%s", Urls.output_prefix, "angle_distribution.dat");
	out_file = fopen(out_url, "w");
	nTot = bins-1;

	fprintf(out_file, "#");
		FILE_HEADER(out_file, "costh  ", count);
		FILE_HEADER(out_file, "P(cth) ", count);
		FILE_HEADER(out_file, "cosphi ", count);
		FILE_HEADER(out_file, "P(cphi)", count);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)
			{
				fprintf(out_file, "\t%lf", costh[i]);
				fprintf(out_file, "\t%lf", (float)p_costh[i]/(float)NsubTh); 
				fprintf(out_file, "\t%lf", cosphi[i]);
				fprintf(out_file, "\t%lf", (float)p_cosphi[i]/(float)NsubTot); 
				fprintf(out_file, "\n");
			}


	fclose(out_file);

	fprintf(out_file, "#");
		FILE_HEADER(out_file, "costh  ", count);
		FILE_HEADER(out_file, "P(cth) ", count);
		FILE_HEADER(out_file, "cosphi ", count);
		FILE_HEADER(out_file, "P(cphi)", count);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)
			{
				fprintf(out_file, "\t%lf", costh[i]);
				fprintf(out_file, "\t%lf", (float)p_costh[i]/(float)NsubTh); 
				fprintf(out_file, "\t%lf", cosphi[i]);
				fprintf(out_file, "\t%lf", (float)p_cosphi[i]/(float)NsubTot); 
				fprintf(out_file, "\n");
			}


	fclose(out_file);
}

	/* Subhaloes for individual host */ 
void print_sub_per_host(int index, int bins, int NsubTh, int host, 
	double *r, int *r_n, int *r_n_c, double *r_v, double *r_m, 
	double *r_m_c, double *costh, int *p_costh, double *cosphi, int *p_cosphi, 
	double *all_r, double *all_m, double *all_v, double *all_cosphi)
{
	count = 1;
	sprintf(out_url, "%shalo.%02d_%s", Urls.output_prefix, host, "distribution.dat");
	out_file = fopen(out_url, "w");
	nTot = bins-1;
	struct halo * HALO;

	HALO = Haloes;
	int Nsub = HALO[index].n_satellites;

	DUMP_MSG("satellite distribution in host halo", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "r/Rv   ", count);
		FILE_HEADER(out_file, "r Mpc  ", count);
		FILE_HEADER(out_file, "N_sub  ", count);
		FILE_HEADER(out_file, "N_sub_c", count);
		FILE_HEADER(out_file, "V_sub  ", count);
		FILE_HEADER(out_file, "M_sub  ", count);
		FILE_HEADER(out_file, "M_sub_c", count);
		FILE_HEADER(out_file, "costh  ", count);
		FILE_HEADER(out_file, "P(cth) ", count);
		FILE_HEADER(out_file, "cosphi ", count);
		FILE_HEADER(out_file, "P(cphi)", count);
		fprintf(out_file, "\n");

		fprintf(out_file, "#Host\t\t");
		fprintf(out_file, "NsubTh=%d\t", NsubTh);
		fprintf(out_file, "Nsub=%d\t", HALO[index].n_satellites);
		fprintf(out_file, "M=%e\t", HALO[index].Mvir);
		fprintf(out_file, "R=%f\t", HALO[index].Rvir);
		fprintf(out_file, "X=%f\t", HALO[index].X[0]);
		fprintf(out_file, "Y=%f\t", HALO[index].X[1]);
		fprintf(out_file, "Z=%f\t", HALO[index].X[2]);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)
			{
				fprintf(out_file, "%lf", r[i]);
				fprintf(out_file, "\t%lf", r[i]*HALO[index].Rvir);
				fprintf(out_file, "\t%d\t", r_n[i]); 
				fprintf(out_file, "\t%d\t", r_n_c[i]); 
				fprintf(out_file, "\t%lf", r_v[i]);
				fprintf(out_file, "\t%lf", r_m[i]);
				fprintf(out_file, "\t%lf", r_m_c[i]);
				fprintf(out_file, "\t%lf", costh[i]);
				fprintf(out_file, "\t%lf", (float)p_costh[i]/(float)NsubTh);
				fprintf(out_file, "\t%e ", cosphi[i]);
				fprintf(out_file, "\t%e ", (float)p_cosphi[i]/(float)r_n_c[0]); 
				fprintf(out_file, "\n");
			}
	fclose(out_file);

	sprintf(out_url, "%shalo.%02d_all_%s", Urls.output_prefix, host, "distribution.dat");
	out_file = fopen(out_url, "w");
	count = 1;
	nTot = Nsub;

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "r/Rv   ", count);
		FILE_HEADER(out_file, "M_sub  ", count);
		FILE_HEADER(out_file, "V_sub  ", count);
		FILE_HEADER(out_file, "costh  ", count);
		fprintf(out_file, "\n");

		fprintf(out_file, "#Host\t");
		fprintf(out_file, "Nsub=%d\t", NsubTh);
		fprintf(out_file, "Nsub=%d\t", HALO[index].n_satellites);
		fprintf(out_file, "M=%e\t", HALO[index].Mvir);
		fprintf(out_file, "R=%f\t", HALO[index].Rvir);
		fprintf(out_file, "X=%f\t", HALO[index].X[0]);
		fprintf(out_file, "Y=%f\t", HALO[index].X[1]);
		fprintf(out_file, "Z=%f\t", HALO[index].X[2]);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)
			{
				fprintf(out_file, "%f", all_r[i]);
				fprintf(out_file, "\t%lf", all_m[i]);
				fprintf(out_file, "\t%lf", all_v[i]);
				fprintf(out_file, "\t%lf", all_cosphi[i]);
				fprintf(out_file, "\n");
			}
	
	fclose(out_file);
}



void print_halo_best_fit_results()
{
	sprintf(out_url, "%s%s", Urls.output_prefix, "_halo_best_fit.dat");
	DUMP_MSG("best fit", out_url);
	out_file=fopen(out_url,"w");

	double MassFrac = Settings.totSubMass / Settings.totHaloMass;
	
		fprintf(out_file, "# Mass fraction in haloes :%e \n", Settings.totHaloMass);
		fprintf(out_file, "# Mass fraction in dark haloes :%e \n", Settings.totDarkMass);
		fprintf(out_file, "# Total gas fraction in haloes :%e \n", Settings.totGasMassInHalo);
		fprintf(out_file, "# Mass-concentration best fit values, c_0=%f, c_beta  :%f \n", 
			HaloProperties[HALO_INDEX].c_0,HaloProperties[HALO_INDEX].c_beta );
		fprintf(out_file, "# Mass-velocity best fit values, vel_0=%f, vel_beta  :%f \n", 
			HaloProperties[HALO_INDEX].vel_0,HaloProperties[HALO_INDEX].vel_beta );
		fprintf(out_file, "# Lambda parameter distribution best fit values, l_0  :%f, l_sig  :%f \n", 
			HaloProperties[HALO_INDEX].l_0,HaloProperties[HALO_INDEX].l_sig );
		fprintf(out_file, "# Average shape parameter s_0  :%lf, average triaxiality t_0    :%lf \n", 
			HaloProperties[HALO_INDEX].halo.s0, HaloProperties[HALO_INDEX].halo.t0);
		fprintf(out_file, "\n");

#ifdef GAS
			fprintf(out_file, "#Mass-T best fit values, T0: %lf, alpha :%lf  \n", 
				HaloProperties[HALO_INDEX].T0, HaloProperties[HALO_INDEX].alpha);
			fprintf(out_file, "#Average beta:%lf\n", HaloProperties[HALO_INDEX].beta0);
			fprintf(out_file, "#Average gamma:%lf\n", HaloProperties[HALO_INDEX].gamma0);
			fprintf(out_file, "#Average delta hydro mass:%lf\n", HaloProperties[HALO_INDEX].dM0);
#endif
	fclose(out_file);
}
