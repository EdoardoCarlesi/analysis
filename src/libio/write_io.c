#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../libhalo/halo.h"
#include "../libcosmo/cosmo.h"
#include "../general_def.h"

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
				fprintf(out_file, "%lf"  , SubHaloProperties[HALO_INDEX].l[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].p_l[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].shape[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].triax[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].costh[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].costh_count[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].cosphi[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].cosphi_count[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].r_sub[i]);
				fprintf(out_file, "\t%d" , SubHaloProperties[HALO_INDEX].n_r_sub[i]);
				fprintf(out_file, "\t%d" , SubHaloProperties[HALO_INDEX].cum_n_r_sub[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].vel_sub[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].p_vel_sub[i]);
				fprintf(out_file, "\t%e" , SubHaloProperties[HALO_INDEX].mass_sub[i]);
				fprintf(out_file, "\t%d" , SubHaloProperties[HALO_INDEX].cum_n_sub[i]);
				fprintf(out_file, "\t%lf", SubHaloProperties[HALO_INDEX].r_sub_subset[i]);
				fprintf(out_file, "\t%d" , SubHaloProperties[HALO_INDEX].cum_n_r_sub_subset[i]);
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
				fprintf(out_file, "\t%lf", HaloProperties[i].l_0);   
				fprintf(out_file, "\t%lf", HaloProperties[i].s0);   
				fprintf(out_file, "\t%lf", HaloProperties[i].t0);   
				fprintf(out_file, "\t%lf", HaloProperties[i].avgSub);   
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
	nTot = HaloProperties[HALO_INDEX].n_bins-1;
	z = GrowthFac.z[Settings.use_cat];
	sprintf(out_url, "%sz%.2f%s", Urls.output_prefix, z, "_all_halo_statistical_properties.dat");
	out_file = fopen(out_url,"w");

	DUMP_MSG("halo properties", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "Mass  ", count);
		FILE_HEADER(out_file, "vel   ", count);
		FILE_HEADER(out_file, "conc  ", count);
		FILE_HEADER(out_file, "lambda", count);
		FILE_HEADER(out_file, "shape ", count);
		FILE_HEADER(out_file, "triax ", count);
		FILE_HEADER(out_file, "chinfw", count);
		FILE_HEADER(out_file, "pernfw", count);
		FILE_HEADER(out_file, "gofnfw", count);
		FILE_HEADER(out_file, "c_nfw ", count);
		FILE_HEADER(out_file, "pc_nfw", count);
		FILE_HEADER(out_file, "g_nfw ", count);
		FILE_HEADER(out_file, "pg_nfw", count);
		FILE_HEADER(out_file, "p_nfw ", count);
		FILE_HEADER(out_file, "pp_nfw", count);
		FILE_HEADER(out_file, "c     ", count);
		FILE_HEADER(out_file, "P(c)  ", count);
		FILE_HEADER(out_file, "l     ", count);
		FILE_HEADER(out_file, "P(l)  ", count);
		FILE_HEADER(out_file, "t     ", count);
		FILE_HEADER(out_file, "P(t)  ", count);
		FILE_HEADER(out_file, "s     ", count);
		FILE_HEADER(out_file, "P(s)  ", count);
#ifdef GAS
		FILE_HEADER(out_file, "gas_temp", count);
		FILE_HEADER(out_file, "gas_frac", count);
#endif
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)	
			{

				fprintf(out_file, "%e",    HaloProperties[HALO_INDEX].mass[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].vel[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].conc[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].lambda[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].shape[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].triax[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].fit_nfw.chi[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].fit_nfw.per[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].fit_nfw.gof[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_fit_nfw.chi[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_fit_nfw.p_chi[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_fit_nfw.per[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_fit_nfw.p_per[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_fit_nfw.gof[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_fit_nfw.p_gof[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].c[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_c[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].l[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_l[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_t[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].s[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].p_s[i]);
#ifdef GAS
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas_T[i]);
				fprintf(out_file, "\t%lf", HaloProperties[HALO_INDEX].gas_fraction[i]);

#endif
				fprintf(out_file, "\n");
			}

	fclose(out_file);
}



void print_grid_CIC()
{
/*
	count = 1;
	nTot = Grid.N;
	sprintf(out_url, "%s%d%s", "../output/halo_cic_grid_N_", nTot, "_cells.dat");
	out_file = fopen(out_url, "w");

	DUMP_MSG("CIC assignment", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "M_CIC  ", count);
		FILE_HEADER(out_file, "rho_CIC", count);
		FILE_HEADER(out_file, "x    ", count);
		FILE_HEADER(out_file, "y    ", count);
		FILE_HEADER(out_file, "z    ", count);
		fprintf(out_file, "\n");

			for(i=0; i<nTot*nTot*nTot; i++) 
			{
				if(Node[i].M_CIC > 0)
				{
					fprintf(out_file, "%e", Node[i].M_CIC);
					fprintf(out_file, "\t%lf", Node[i].M_CIC / Grid.cell_volume);
					fprintf(out_file, "\t%lf", Node[i].X[0]);
					fprintf(out_file, "\t%lf", Node[i].X[1]);
					fprintf(out_file, "\t%lf", Node[i].X[2]);
					fprintf(out_file, "\n"); 
				}
			}

	fclose(out_file);
*/
}



void print_halo_density()
{
/*
	count = 1;
	nTot = Density.N;
	sprintf(out_url, "%sN_%d%s", "../output/halo_density_test_", nTot, "_cells_.dat");
	out_file = fopen(out_url, "w");

	DUMP_MSG("Halo density", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "r      ", count);
		FILE_HEADER(out_file, "rho    ", count);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++) 
			{
				fprintf(out_file, "%e", Density.r_n[i]);
				fprintf(out_file, "\t%e", Density.rho[i]);
				fprintf(out_file, "\n"); 
			}

	fclose(out_file);
*/
}




void print_axis_alignment()
{
	count = 1;
	nTot = HaloProperties[HALO_INDEX].r_bins-1;
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
				fprintf(out_file, "%e", MassFunc[MF_INDEX].mass[i]);
				fprintf(out_file, "\t%e", MassFunc[MF_INDEX].n[i]);
				fprintf(out_file, "\t%d\t", MassFunc[MF_INDEX].n_tot[i]);
				fprintf(out_file, "\t%e", MassFunc[MF_INDEX].err[i]);
				fprintf(out_file, "%e", MassFunc[MF_INDEX].mass_halfstep[i]);
				fprintf(out_file, "\t%e", MassFunc[MF_INDEX].dn[i]);
				fprintf(out_file, "\t%d\t", MassFunc[MF_INDEX].n_bin[i]);
				fprintf(out_file, "\t%e", MassFunc[MF_INDEX].err_dn[i]);
				fprintf(out_file,"\n");
			}

	fclose(out_file);
}


/*
void print_nfw()
{
// TODO //FIXME
	nTot = NFW.bins;		
	sprintf(out_url, "%s%s", Urls.output_prefix, "nfw_test.dat");
	out_file = fopen(out_url, "w");

	DUMP_MSG("Navarro Frenk White profile", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "r", count);
		FILE_HEADER(out_file, "overd", count);
		FILE_HEADER(out_file, "rho_NFW", count);
		FILE_HEADER(out_file, "err", count);
		fprintf(out_file, "\n");

			for(i=0; i<nTot; i++)
			{
				fprintf(out_file, "%lf", NFW.radius[i]);
				fprintf(out_file, "\t%lf", NFW.overd[i]);
				fprintf(out_file, "\t%lf", NFW.profile[i]);
				fprintf(out_file, "\t%lf", NFW.err[i]);
				fprintf(out_file, "\n");
			}

	fclose(out_file);
}
*/


void print_best_fit_results(){
// TODO //FIXME
	out_file=fopen(out_url,"w");
	sprintf(out_url, "%s%s", Urls.output_prefix,"all_halo_best_fit_distributions.dat");
	DUMP_MSG("best fit", out_url);
	
	//	fprintf(out_file, "#Concentration distribution best fit values, c_0=%lf, c_sig  :%lf \n", 
	//		HaloProperties[HALO_INDEX].c_0,HaloProperties[HALO_INDEX].c_sig );
		fprintf(out_file, "#Lambda parameter distribution best fit values, l_0  :%lf, l_sig  :%lf \n", 
			HaloProperties[HALO_INDEX].l_0,HaloProperties[HALO_INDEX].l_sig );
		fprintf(out_file, "#Average shape parameter s_0  :%lf, average triaxiality t_0    :%lf \n", 
			HaloProperties[HALO_INDEX].s0,HaloProperties[HALO_INDEX].t0);

		fprintf(out_file, "#lambda_0  = %lf\n", SubHaloProperties[HALO_INDEX].l_0);
		fprintf(out_file, "#cos(th_0) = %lf\n", SubHaloProperties[HALO_INDEX].costh0);
		fprintf(out_file, "#cos(phi_0)= %lf\n", SubHaloProperties[HALO_INDEX].cosphi0);
		fprintf(out_file, "#t0 = %lf\n", SubHaloProperties[HALO_INDEX].t0);
		fprintf(out_file, "#s0 = %lf\n", SubHaloProperties[HALO_INDEX].s0);
		fprintf(out_file, "\n");

#ifdef GAS
			fprintf(out_file, "#Mass-Tx relation best fit values, ln_M0: %lf, alpha :%lf  \n", 
				HaloProperties[HALO_INDEX].ln_M0, HaloProperties[HALO_INDEX].T_alpha);
#endif
}
