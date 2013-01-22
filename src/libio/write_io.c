#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../libmath/log_norm.h"
#include "../general_variables.h"
#include "../general_functions.h"

#define FILE_HEADER(fout, info, count) fprintf(fout, info "(%i)\t",count++)


void print_number_densities()
{
		int k=0, count=1, N=NumDen.npts; 
		double M=Settings.mass_min;
		char *out  = merge_strings(Urls.output_prefix,"number_density.dat");
		FILE *fout = fopen(out, "w");

		fprintf(stdout, "\nWriting number_density() output to: %s \n", out);
	
		fprintf(fout,"#");
		FILE_HEADER(fout, "z       ", count);
		FILE_HEADER(fout, "Tinker n", count);
		FILE_HEADER(fout, "numer. n", count);
		fprintf(fout, "\n");
//		fprintf(fout, "#MassCut:%e Solar Masses\n", M);

			for(k=0; k<N; k++) 
			{
			fprintf(fout,"%lf", NumDen.z[k]);
			fprintf(fout,"\t%e", NumDen.n_tin[k]);
			fprintf(fout,"\t%e", NumDen.n_num[k]);
			fprintf(fout, "\n");
			}

	fclose(fout);
}



void print_all_subhalo_properties_to_one_file()
{
	int i=0, count=1, nBins=SubHaloZ.n_bins;
	char* out  = merge_strings(Urls.output_prefix, "all_sub_statistics.dat");
	FILE* fout = fopen(out,"w");

	fprintf(stdout, "\nWriting subhalo_properties() output to:%s\n", out);

		fprintf(fout,"#");
		FILE_HEADER(fout, "lambda       ", count);
		FILE_HEADER(fout, "P(l)         ", count);
		FILE_HEADER(fout, "shape        ", count);
		FILE_HEADER(fout, "P(s)         ", count);
		FILE_HEADER(fout, "triax        ", count);
		FILE_HEADER(fout, "P(t)         ", count);
		FILE_HEADER(fout, "Cos(th)      ", count);
		FILE_HEADER(fout, "P(c(th))     ", count);
		FILE_HEADER(fout, "Cos(phi)     ", count);
		FILE_HEADER(fout, "P(c(phi))    ", count);
		FILE_HEADER(fout, "sub(r)/R_vir ", count);
		FILE_HEADER(fout, "n(r)         ", count);
		FILE_HEADER(fout, "n(>r)        ", count);
		FILE_HEADER(fout, "V_sub        ", count);
		FILE_HEADER(fout, "P(V_sub)     ", count);
		FILE_HEADER(fout, "Mass         ", count);
		FILE_HEADER(fout, "n(>M)        ", count);
		FILE_HEADER(fout, "subset(r)    ", count);
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



void print_virial_theorem()
{	
	int k=0, nTot=0, count=1; 
	char *fName=merge_strings(Urls.output_prefix,"virial_theorem.dat");
	FILE *fout = fopen(fName, "w");

		fprintf(fout,"#");
		FILE_HEADER(fout, "Mass", count);
		FILE_HEADER(fout, "VirialRatio", count);
		fprintf(fout, "\n");

		nTot = Settings.n_haloes;

			for(k=0; k<nTot; k++) 
			{
				fprintf(fout,"%e", Haloes[k].Mvir);
				fprintf(fout,"\t%lf", Haloes[k].abs_th_vir);
				fprintf(fout, "\n");
			}

	fclose(fout);
}



void print_theoretical_mass_function(float z)
{
	double diff_tin=0, tin=0, n_num=0, dn_num=0, n_mass=0, err=0, th_mass=0;
	int k=0, count=1, nTot=0; 
	char out_file[200];
	FILE *fout=NULL;

	nTot=ThMassFunc.bins-1;

	sprintf(out_file, "%s%f%s", 
		Urls.output_prefix, z, "_theoretical_mass_function.dat");

	fout=fopen(out_file, "w");

	fprintf(stderr, "Mass function output file:%s \n", out_file);

		fprintf(fout,"#");
		FILE_HEADER(fout, "Mass        ", count);
		FILE_HEADER(fout, "dn Tinker   ", count);
		FILE_HEADER(fout, "n(>M) Tinker", count);
		fprintf(fout, "\n");

			for(k=0; k<nTot; k++)
			{
				th_mass=ThMassFunc.th_masses[k];
				tin=ThMassFunc.tin[k];
				diff_tin=ThMassFunc.diff_tin[k];	
				fprintf(fout,"%e", th_mass);
				fprintf(fout,"\t%e", diff_tin);
				fprintf(fout,"\t%e", tin);
				fprintf(fout,"\n");
			}

	fclose(fout);
}



void print_correlation_function()
{
	int i=0, count=1; 
	char *fileName; 
	FILE *fout;

		fileName=merge_strings(Urls.output_prefix, "correlation_function.dat");
		fout=fopen(fileName,"w");

		fprintf(fout,"#");
		FILE_HEADER(fout, "R     ", count);
		FILE_HEADER(fout, "Xi    ", count);
		FILE_HEADER(fout, "Xi fit", count);
		fprintf(fout, "\n");
	
		for(i=0; i<Xi.n_xi_entries; i++) 
		{
			fprintf(fout, "%lf\n", Xi.r[i]); 
			fprintf(fout, "\t%lf\n", Xi.xi_r[i]);
			fprintf(fout, "\t%lf\n", Xi.xi_fit[i]);
			fprintf(fout, "\n");
		}
fclose(fout);
}



void print_nfw()
{
	int k=0, count=1, nTot=0;
	double rr, od, rs, r0, nnffww;
	FILE * out_nfw = fopen(merge_strings(Urls.output_prefix, "nfw_test.dat"), "w");

	nTot = NFW.bins;		

	fprintf(out_nfw, "#");
	FILE_HEADER(out_nfw, "r        ", count);
	FILE_HEADER(out_nfw, "rho_{NFW}", count);
	FILE_HEADER(out_nfw, "overd    ", count);
	FILE_HEADER(out_nfw, "err      ", count);
	fprintf(out_nfw, "\n");

		for(k=0; k<nTot; k++)
		{
			rr=NFW.radius[k];
			od=NFW.overd[k];
			rs=NFW.rs;
			r0=NFW.rho0;
			nnffww = nfw(rr, rs, r0);
	
			fprintf(out_nfw, "%lf", rr);
			fprintf(out_nfw, "\t%lf", nnffww);
			fprintf(out_nfw, "\t%lf", od);
			fprintf(out_nfw, "\t%lf", NFW.err[k]);
			fprintf(out_nfw, "\n");
		}

	fclose(out_nfw);
}



void print_growth_factor()
{
	int count=1, k=0, dim_eff=0; 
	char *fileName; 
	FILE *output=NULL;

	fprintf(stderr,"print_growth_factors().\n");

	dim_eff = GrowthFac.npts;
	fileName=merge_strings(Urls.output_prefix, "growth_factor.dat");
	output=fopen(fileName,"w");
	
		fprintf(output, "#");
		FILE_HEADER(output, "z   ", count);
		FILE_HEADER(output, "gf  ", count);
		FILE_HEADER(output, "gf/a", count);
		fprintf(output, "\n");
	
		for (k=0; k<dim_eff; k++) 
		{
			fprintf(output,"%lf", GrowthFac.z[k]);
			fprintf(output,"\t%lf", GrowthFac.gf_z[k]);
			fprintf(output,"\t%lf", GrowthFac.gf_over_a_z[k]);
			fprintf(output,"\n");
		}

	fclose(output);
	fprintf(stderr, "Written growth factor file to %s.\n", fileName);
}



void print_evolution_to_file()
{
	int j=0, count=1, nTot=0; 
	char *out = merge_strings(Urls.output_prefix, "halo_subhalo_evolution.dat");
	FILE *f_out = fopen(out, "w");

	nTot=FullCat.numFiles-1;
	
	fprintf(f_out,"#");
	FILE_HEADER(f_out, "z           ", count);
	FILE_HEADER(f_out, "conc        ", count);
	FILE_HEADER(f_out, "lambda      ", count);
	FILE_HEADER(f_out, "shape       ", count);
	FILE_HEADER(f_out, "triax       ", count);
	FILE_HEADER(f_out, "avg_sub     ", count);
	FILE_HEADER(f_out, "sub_conc    ", count);
	FILE_HEADER(f_out, "sub_lambda  ", count);
	FILE_HEADER(f_out, "sub_shape   ", count);
	FILE_HEADER(f_out, "sub_triax   ", count);
	FILE_HEADER(f_out, "sub_costh   ", count);
	FILE_HEADER(f_out, "sub_cosphi  ", count);
	FILE_HEADER(f_out, "sub_avg_vel ", count);
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



void print_all_halo_properties_to_one_file()
{
	int i=0, count=1, nBins=HaloZ.n_bins;
	char* out_url=merge_strings(Urls.output_prefix,"all_halo_statistical_properties.dat");
	FILE* out=fopen(out_url,"w");
	
		fprintf(out, "#");
		FILE_HEADER(out, "Mass  ", count);
		FILE_HEADER(out, "n(>M) ", count);
		FILE_HEADER(out, "avg_c ", count);
		FILE_HEADER(out, "conc  ", count);
		FILE_HEADER(out, "P(c)  ", count);
		FILE_HEADER(out, "lambda", count);
		FILE_HEADER(out, "P(l)  ", count);
		FILE_HEADER(out, "triax ", count);
		FILE_HEADER(out, "P(t)  ", count);
		FILE_HEADER(out, "shape ", count);
		FILE_HEADER(out, "P(s)  ", count);
#ifdef GAS
		FILE_HEADER(out, "gas_temp", count);
		FILE_HEADER(out, "gas_frac", count);
#endif
		fprintf(out, "\n");

			for(i=0; i<nBins; i++)	
			{
				fprintf(out, "%e", MassFunc.num_masses[i]);
				fprintf(out, "\t%e", MassFunc.n[i]);
				fprintf(out, "\t%lf", HaloZ.c_avg_mass[i]);
				fprintf(out, "\t%lf", HaloZ.c[i]);
				fprintf(out, "\t%lf", HaloZ.p_c[i]);
				fprintf(out, "\t%lf", HaloZ.l[i]);
				fprintf(out, "\t%lf", HaloZ.p_l[i]);
				fprintf(out, "\t%lf", HaloZ.triax[i]);
				fprintf(out, "\t%lf", HaloZ.p_triax[i]);
				fprintf(out, "\t%lf", HaloZ.shape[i]);
				fprintf(out, "\t%lf", HaloZ.p_shape[i]);
#ifdef GAS
				fprintf(out, "\t%lf", HaloZ.gas_T[i]);
				fprintf(out, "\t%lf", HaloZ.gas_fraction[i]);

#endif
				fprintf(out, "\n");
		}

	fclose(out);
}



void print_axis_alignment()
{
	int i=0, count=1, nBins=HaloZ.r_bins;
	char *file_out = merge_strings(Urls.output_prefix,"axis_alignement.dat");
	FILE *output = fopen(file_out, "w");

	fprintf(output, "#");
	FILE_HEADER(output, "R", count);
	FILE_HEADER(output, "cos(Th_c)^2", count);
	FILE_HEADER(output, "cos(Th_p)^2", count);
	FILE_HEADER(output, "N", count);
	fprintf(output, "\n");

		for(i=0; i<nBins; i++) 
		{
			fprintf(output, "%lf", HaloZ.R[i]);
			fprintf(output, "\t%lf", HaloZ.Th_c[i]/HaloZ.N_pairs[i]);
			fprintf(output, "\t%lf", HaloZ.Th_p[i]/HaloZ.N_pairs[i]);
			fprintf(output, "\t%d", HaloZ.N_pairs[i]);
			fprintf(output, "\n"); 
		}

	fclose(output);
}



void print_numerical_mass_function()
{
	int i=0, count=1, nBins=MassFunc.bins;
	char *file_out = merge_strings(Urls.output_prefix,"numerical_mass_function.dat");
	FILE *output = fopen(file_out, "w");

		fprintf(output,"#");
		FILE_HEADER(output, "M", count);
		FILE_HEADER(output, "n", count);
		FILE_HEADER(output, "n_tot", count);
		FILE_HEADER(output, "n_err", count);
		FILE_HEADER(output, "M_dM", count);
		FILE_HEADER(output, "dn", count);
		FILE_HEADER(output, "n_bin", count);
		FILE_HEADER(output, "dn_err", count);
		fprintf(output,"\n");

			for(i=0; i<nBins; i++) 
			{
				fprintf(output, "%e", MassFunc.num_masses[i]);
				fprintf(output, "\t%e", MassFunc.n[i]);
				fprintf(output, "\t%d", MassFunc.n_tot[i]);
				fprintf(output, "\t%e", MassFunc.err[i]);
				fprintf(output, "\t%e", MassFunc.dn_num_masses[i]);
				fprintf(output, "\t%e", MassFunc.dn[i]);
				fprintf(output, "\t%d", MassFunc.n_bin[i]);
				fprintf(output, "\t%e", MassFunc.err_dn[i]);
				fprintf(output,"\n");
			}
	fclose(output);
}



void print_best_fit_results(){
// TODO
	int i=0, count=1, nBins=HaloZ.n_bins;
	char* out_url=merge_strings(Urls.output_prefix,"all_halo_best_fit_distributions.dat");
	FILE* fout=fopen(out_url,"w");
	
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
