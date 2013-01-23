#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../libcosmo/nfw.h"
#include "../libmath/log_norm.h"
#include "../general_variables.h"
#include "../general_functions.h"

#define FILE_HEADER(fout, info, count) fprintf(fout, info "(%i)\t",count++)
#define DUMP_MSG(name, url) fprintf(stdout, "Printing %s to file: %s\n", name, url)


void print_number_densities()
{
	int i=0, count=1, nTot=0, N=0;
	double M=0;
	char out_url[200]; 
	FILE *out_file=NULL;

	nTot = NumDen.npts; 
	M = Settings.mass_min;
	N = Settings.n_min;

	if(Settings.use_n_min == 1)
		sprintf(out_url, "%sN_%d%s", Urls.output_prefix, N, "_number_density.dat");
	else
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



void print_all_subhalo_properties_to_one_file(int j)
{
	int i=0, count=1, nBins=0;
	double z=0;
	char out[200]; 
	FILE* fout=NULL;

	nBins = SubHaloZ.n_bins;
	z = GrowthFac.z[j];
 	sprintf(out, "%sz%.2f%s", Urls.output_prefix, z, "_all_sub_statistics.dat");
	fout = fopen(out,"w");

	DUMP_MSG("subhalo properties", out);

		fprintf(fout,"#");
		FILE_HEADER(fout, "lambda ", count);
		FILE_HEADER(fout, "P(l)", count);
		FILE_HEADER(fout, "shape", count);
		FILE_HEADER(fout, "P(s)", count);
		FILE_HEADER(fout, "triax", count);
		FILE_HEADER(fout, "P(t)", count);
		FILE_HEADER(fout, "Cos(th)", count);
		FILE_HEADER(fout, "P(c(th))", count);
		FILE_HEADER(fout, "Cos(phi)", count);
		FILE_HEADER(fout, "P(c(phi))", count);
		FILE_HEADER(fout, "sub(r)/Rv", count);
		FILE_HEADER(fout, "n(r)", count);
		FILE_HEADER(fout, "n(>r)", count);
		FILE_HEADER(fout, "V_sub", count);
		FILE_HEADER(fout, "P(V_sub)", count);
		FILE_HEADER(fout, "mass", count);
		FILE_HEADER(fout, "n(>M)", count);
		FILE_HEADER(fout, "subset(r)", count);
		FILE_HEADER(fout, "subset(n(>r))", count);
		fprintf(fout, "\n");

		for(i=0; i<nBins; i++)	
		{
			fprintf(fout, "%lf"  , SubHaloZ.l[i]);
			fprintf(fout, "\t%lf", SubHaloZ.p_l[i]);
			fprintf(fout, "\t%lf", SubHaloZ.shape[i]);
			fprintf(fout, "\t%lf", SubHaloZ.triax[i]);
			fprintf(fout, "\t%lf", SubHaloZ.costh[i]);
			fprintf(fout, "\t%lf", SubHaloZ.costh_count[i]);
			fprintf(fout, "\t%lf", SubHaloZ.cosphi[i]);
			fprintf(fout, "\t%lf", SubHaloZ.cosphi_count[i]);
			fprintf(fout, "\t%lf", SubHaloZ.r_sub[i]);
			fprintf(fout, "\t%d" , SubHaloZ.n_r_sub[i]);
			fprintf(fout, "\t%d" , SubHaloZ.cum_n_r_sub[i]);
			fprintf(fout, "\t%lf", SubHaloZ.vel_sub[i]);
			fprintf(fout, "\t%lf", SubHaloZ.p_vel_sub[i]);
			fprintf(fout, "\t%e" , SubHaloZ.mass_sub[i]);
			fprintf(fout, "\t%d" , SubHaloZ.cum_n_sub[i]);
			fprintf(fout, "\t%lf", SubHaloZ.r_sub_subset[i]);
			fprintf(fout, "\t%d" , SubHaloZ.cum_n_r_sub_subset[i]);
			fprintf(fout, "\n");
		}

	fclose(fout);	
}



void print_theoretical_mass_function(int j)
{
	double z=0; 
	int k=0, count=1, nTot=0; 
	char out_file[200];
	FILE *fout=NULL;

	nTot = ThMassFunc.bins-1;
	z = GrowthFac.z[j];
	sprintf(out_file, "%sz%.2f%s", Urls.output_prefix, z, "_theoretical_mass_function.dat");
	fout=fopen(out_file, "w");

	DUMP_MSG("theoretical mass function", out_file);

		fprintf(fout,"#");
		FILE_HEADER(fout, "Mass", count);
		FILE_HEADER(fout, "n(>M) Tinker", count);
		FILE_HEADER(fout, "dn Tinker", count);
		fprintf(fout, "\n");

			for(k=0; k<nTot; k++)
			{
				fprintf(fout,"%e", ThMassFunc.mass_halfstep[k]);
				fprintf(fout,"\t%e", ThMassFunc.n[k]);
				fprintf(fout,"\t%e", ThMassFunc.dn[k]);
				fprintf(fout,"\n");
			}

	fclose(fout);
}



void print_correlation_function(int j)
{
	int i=0, count=1, nTot=0; 
	double z=0;
	char fileName[200];
	FILE *fout=NULL;
	
	nTot = Xi.npts;
	z = GrowthFac.z[j];
	sprintf(fileName, "%sz%.2f%s", Urls.output_prefix, z, "_correlation_function.dat");
	fopen(fileName,"w");

	DUMP_MSG("correlation function", fileName);

		fprintf(fout,"#");
		FILE_HEADER(fout, "R", count);
		FILE_HEADER(fout, "Xi", count);
		FILE_HEADER(fout, "Xi fit", count);
		fprintf(fout, "\n");
	
		for(i=0; i<nTot; i++) 
		{
			fprintf(fout, "%lf\n", Xi.r[i]); 
			fprintf(fout, "\t%lf\n", Xi.xi_r[i]);
			fprintf(fout, "\t%lf\n", Xi.xi_fit[i]);
			fprintf(fout, "\n");
		}
	fclose(fout);
}



void print_growth_factor()
{
	int count=1, k=0, nTot=0; 
	double scale=0;
	char fileName[200];
	FILE *output=NULL; 

	nTot = GrowthFac.npts;
	scale = GrowthFac.scale_k;
	sprintf(fileName, "%sk%.3f%s", Urls.output_prefix, scale, "_growth_factor.dat");	
	output = fopen(fileName,"w");

	DUMP_MSG("growth factor", fileName);
	
		fprintf(output, "#");
		FILE_HEADER(output, "z", count);
		FILE_HEADER(output, "gf", count);
		FILE_HEADER(output, "gf/a", count);
		fprintf(output, "\n");
	
		for (k=0; k<nTot; k++) 
		{
			fprintf(output,"%lf", GrowthFac.z[k]);
			fprintf(output,"\t%lf", GrowthFac.gf[k]);
			fprintf(output,"\t%lf", GrowthFac.gf[k]/GrowthFac.a[k]);
			fprintf(output,"\n");
		}

	fclose(output);
}



void print_evolution_to_file()
{
	int j=0, count=1, nTot=0; 
	char out[200]; 
	FILE *f_out=NULL;

	nTot = Urls.nCatalogueFiles-1;
	sprintf(out, "%s%s", Urls.output_prefix, "halo_subhalo_evolution.dat");
	fopen(out, "w");

	DUMP_MSG("evolution", out);
	
		fprintf(f_out,"#");
		FILE_HEADER(f_out, "z", count);
		FILE_HEADER(f_out, "conc", count);
		FILE_HEADER(f_out, "lambda", count);
		FILE_HEADER(f_out, "shape", count);
		FILE_HEADER(f_out, "triax", count);
		FILE_HEADER(f_out, "avg_sub", count);
		FILE_HEADER(f_out, "sub_conc", count);
		FILE_HEADER(f_out, "sub_lambda", count);
		FILE_HEADER(f_out, "sub_shape", count);
		FILE_HEADER(f_out, "sub_triax", count);
		FILE_HEADER(f_out, "sub_costh", count);
		FILE_HEADER(f_out, "sub_cosphi", count);
		FILE_HEADER(f_out, "sub_avg_vel", count);
		FILE_HEADER(f_out, "sub_avg_mass", count);
		fprintf(f_out,"\n");

		for(j=0; j<nTot; j++)
		{ 
			fprintf(f_out, "%lf", HaloProperties[j].z);   
			fprintf(f_out, "\t%lf", HaloProperties[j].c_0);   
			fprintf(f_out, "\t%lf", HaloProperties[j].l_0);   
			fprintf(f_out, "\t%lf", HaloProperties[j].s0);   
			fprintf(f_out, "\t%lf", HaloProperties[j].t0);   
			fprintf(f_out, "\t%lf", HaloProperties[j].avgSub);   
			fprintf(f_out, "\t%lf", SubHaloProperties[j].c_0);   	
			fprintf(f_out, "\t%lf", SubHaloProperties[j].l_0);   	
			fprintf(f_out, "\t%lf", SubHaloProperties[j].s0);   	
			fprintf(f_out, "\t%lf", SubHaloProperties[j].t0);   	
			fprintf(f_out, "\t%lf", SubHaloProperties[j].costh0);   	
			fprintf(f_out, "\t%lf", SubHaloProperties[j].cosphi0);   	
			fprintf(f_out, "\t%lf", SubHaloProperties[j].vel_0);   	
			fprintf(f_out, "\t%e", SubHaloProperties[j].avgMass);   	
			fprintf(f_out,"\n");
		}

	fclose(f_out);
}



void print_all_halo_properties_to_one_file(int j)
{
	int i=0, count=1, nBins=0;
	double z=0;
	char out_url[200];
	FILE* out_file=NULL;
	
	nBins=HaloZ.n_bins-1;
	z = GrowthFac.z[j];
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

			for(i=0; i<nBins; i++)	
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



void print_axis_alignment(int j)
{
	int i=0, count=1, nBins=0;
	double z=0;
	char out_url[200]; 
	FILE *out_file=NULL;

	nBins = HaloZ.r_bins-1;
	z = GrowthFac.z[j];
	sprintf(out_url, "%sz%.2f%s", Urls.output_prefix, z, "_axis_alignement.dat");
	out_file = fopen(out_url, "w");

	DUMP_MSG("axis alignment", out_url);

		fprintf(out_file, "#");
		FILE_HEADER(out_file, "r      ", count);
		FILE_HEADER(out_file, "c(Tc)^2", count);
		FILE_HEADER(out_file, "c(Tp)^2", count);
		FILE_HEADER(out_file, "N      ", count);
		fprintf(out_file, "\n");

			for(i=0; i<nBins; i++) 
			{
				fprintf(out_file, "%lf", HaloZ.R[i]);
				fprintf(out_file, "\t%lf", HaloZ.Th_c[i]/HaloZ.N_pairs[i]);
				fprintf(out_file, "\t%lf", HaloZ.Th_p[i]/HaloZ.N_pairs[i]);
				fprintf(out_file, "\t%d", HaloZ.N_pairs[i]);
				fprintf(out_file, "\n"); 
			}

	fclose(out_file);
}



void print_numerical_mass_function(int j)
{
	int i=0, count=1, nTot=0;
	double z=0;
	char out_url[200]; 
	FILE *out_file=NULL;

	nTot = MassFunc.bins;
	z = GrowthFac.z[j];
	sprintf(out_url, "%sz%.2f%s", Urls.output_prefix, z, "_numerical_mass_function.dat");
	out_file = fopen(out_url, "w");

	DUMP_MSG("numerical mass function", out_url);

		fprintf(out_file,"#");
		FILE_HEADER(out_file, "M    ", count);
		FILE_HEADER(out_file, "n    ", count);
		FILE_HEADER(out_file, "n_tot", count);
		FILE_HEADER(out_file, "n_err", count);
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
	FILE_HEADER(out_nfw, "rho_NFW", count);
	FILE_HEADER(out_nfw, "overd", count);
	FILE_HEADER(out_nfw, "err", count);
	fprintf(out_nfw, "\n");

		for(k=0; k<nTot; k++)
		{
			fprintf(out_nfw, "%lf", NFW.radius[k]);
			fprintf(out_nfw, "\t%lf", NFW.overd[k]);
			fprintf(out_nfw, "\t%lf", nfw(NFW.radius[k], NFW.rs, NFW.rho0));
			fprintf(out_nfw, "\t%lf", NFW.err[k]);
			fprintf(out_nfw, "\n");
		}

	fclose(out_nfw);
}



void print_best_fit_results(){
// TODO //FIXME
	int i=0, count=1, nBins=HaloZ.n_bins;
	char* out_url=merge_strings(Urls.output_prefix,"all_halo_best_fit_distributions.dat");
	FILE* fout=fopen(out_url,"w");
	DUMP_MSG("best fit", out_url);
	
		fprintf(fout, "#Concentration distribution best fit values, c_0=%lf, c_sig  :%lf \n", 
			HaloZ.c_0,HaloZ.c_sig );
		fprintf(fout, "#Lambda parameter distribution best fit values, l_0  :%lf, l_sig  :%lf \n", 
			HaloZ.l_0,HaloZ.l_sig );
		fprintf(fout, "#Average shape parameter s_0  :%lf, average triaxiality t_0    :%lf \n", 
			HaloZ.s0,HaloZ.t0);

		fprintf(fout, "#lambda_0  = %lf\n", SubHaloZ.l_0);
		fprintf(fout, "#cos(th_0) = %lf\n", SubHaloZ.costh0);
		fprintf(fout, "#cos(phi_0)= %lf\n", SubHaloZ.cosphi0);
		fprintf(fout, "#t0 = %lf\n", SubHaloZ.t0);
		fprintf(fout, "#s0 = %lf\n", SubHaloZ.s0);
		fprintf(fout, "\n");

#ifdef GAS
			fprintf(fout, "#Mass-Tx relation best fit values, ln_M0: %lf, alpha :%lf  \n", 
				HaloZ.ln_M0, HaloZ.T_alpha);
#endif
}
