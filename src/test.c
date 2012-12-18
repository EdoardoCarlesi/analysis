/* Test new routines and functions */
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


#include "libio/halo_io.h"
#include "libio/read_io.h"

int main(int argc, char **argv){
fprintf(stderr, "\n*********TEST**********\n");

int *ran; int N=50; int M=30;
ran = (int*)calloc(M,sizeof(int));
generate_random_subset(N,M,ran);

/*
initialize_internal_variables(argv);

char *url = "/home/carlesi/vde1_correlation_function.dat";
FILE *file = fopen(url, "r");

int size = get_lines(file);
int k;

double *x; double *y;
double *err;
char line[256];

x   = (double*) calloc(size,sizeof(double));
err = (double*) calloc(size,sizeof(double));
y   = (double*) calloc(size,sizeof(double));

for(k=0; k<size; k++){
fgets(line,256,file);
sscanf(line,"%lf %lf", &x[k], &y[k]);
err[k]=0.2*y[k];
//fprintf(stderr,"%lf %lf %lf\n", x[k], y[k], err[k]);
}

double* a = best_fit_power_law(x,y,err,size);
char *corr = "vde1_correlation_best_fit.dat";
FILE *cbf = fopen(corr, "w");
fprintf(cbf, "#x(1)  y(2)  err(3)  y_bf(3)\n");

for(k=0; k<size; k++){
double yy = a[1]*pow(x[k],a[0]);
//fprintf(stderr, "x: %lf, y:%lf \n", x[k], yy);
fprintf(cbf, "%lf  %lf  %lf  %lf\n", x[k], y[k], err[k], yy);
}
*/
/*
Cosmo.OmegaM = 0.3;
Cosmo.OmegaL = 0.7;
Settings.model_type=1;
print_H_z("hz_lcdm1.dat");

Cosmo.OmegaM = 0.388;
Cosmo.OmegaL = 0.612;
Settings.model_type=1;
print_H_z("hz_lcdm2.dat");
*/
// TEST Mass Function and similar routines
//Urls_internal.halo_file = "/home/edoardo/UnsortedData/merged_snapshot__29_z0.AHF_halos-lcdm256";
//Urls_internal.halo_file = "/home/carlesi/data/lcdm_256/ahf/merged_snapshot__29_z0.AHF_halos";
//Urls_internal.analysis_dir = "/home/carlesi/Analysis/";
//MF.lines_skip=1;
//MF.n_bins=50;
//Settings.model_type=1;
//Settings.box_size=1;

//initialize_urls();
//fprintf(stderr, "\n URLS initialized \n");
//read_halo_file();
//fprintf(stderr, "\n Halo File Read \n");

//FILE *out;
//out = fopen("test.dat", "w");
//sort_and_print_lambda(out);

// TEST Non linear fit

//struct data *d;
//least_square_nl_fit(d);

//char *url1 = "output/mf_lcdm_500-z1.dat";
//char *url2 = "output/mf_lcdm_500-z1.dat";
//char *url3 = "output/mf_ratios_z1.dat";
//char *url1 = "hz_lcdm1.dat";
//char *url2 = "hz_lcdm2.dat";
//char *url3 = "ratios_lcdms.dat";
//calculate_files_ratio(url2, url1, url3, 3000, 1, 1);

//char *e1url1 = "/home/edoardo/UnsortedData/vde_512/Pk-512-1_snapshot_060";
//char *e1url2 = "/home/edoardo/UnsortedData/lcdm_512/Pk-512-1_snapshot_060-vde";
//char *e1url3 = "/home/edoardo/UnsortedData/vde_512/pk_ratio_vdelcdmvde_05z4.dat";
//char *e1url1 = "/home/edoardo/Analysis/output/41_rp5-scott-code-icnew-dedm_growth_factor.dat";
//char *e1url2 = "/home/edoardo/Gadget-devel/gadget-dedm/Tables/RP5c/growthFactor_Z.dat";
//char *e1url3 = "/home/edoardo/Analysis/output/ratio_41-rp5c.dat";

//char *e1url1 = "/home/edoardo/DEDM_Plots/cDE_1/32/gf.dat"; 
//char *e1url2 = "/home/edoardo/DEDM_Plots/Tables/Om_027/cDE1/growthFactor_Z.dat"; 
//char *e1url2 = "/home/edoardo/DEDM_Plots/Tables/RP3/growthFactor_NormalizedToA.dat"; 
//char *e1url3 = "/home/edoardo/DEDM_Plots/cDE_1/ratio_32-th.dat"; 

/*
char *url1 = "/home/edoardo/DEDM_Plots/Tables/Om_027/cDE3/hubble_factor_z.dat";
char *url2 = "/home/edoardo/DEDM_Plots/Tables/Om_027/cDE3/hubble_vde.dat";
normalize_to_one(url1,url2);
char *url1a = "/home/edoardo/DEDM_Plots/Tables/Om_027/cDE2/hubble_factor_z.dat";
char *url2a = "/home/edoardo/DEDM_Plots/Tables/Om_027/cDE2/hubble_vde.dat";
normalize_to_one(url1a,url2a);

//calculate_files_ratio(e1url2, e1url1, e1url3, 30, 1, 1);
//calculate_gf_ratio(e1url2, e1url1, e1url3, 30, 1, 1);

char *e2url1 = "/home/edoardo/UnsortedData/lcdm_512/Pk-512-1_snapshot_072";
char *e2url2 = "/home/edoardo/UnsortedData/lcdm_512/Pk-512-1_snapshot_072";
char *e2url3 = "../pk_ratio_z2_lcdmlcdm.dat";
calculate_files_ratio(e2url1, e2url2, e2url3, 30, 1, 1);

char *ee2url1 = "/home/edoardo/UnsortedData/lcdm_512/Pk-512-1_snapshot_060";
char *ee2url2 = "/home/edoardo/UnsortedData/lcdm_512/Pk-512-1_snapshot_060-lcdm";
char *ee2url3 = "../pk_ratio_z4_lcdmlcdm-lcdm.dat";
calculate_files_ratio(ee2url1, ee2url2, ee2url3, 30, 1, 1);
*/

// Test P(k) integration routines

//init_cosmology();
//initialize_urls();
//Urls_internal.hubble_dir="/home/edoardo/Analysis/VDE/Tables/hubble_lcdm.dat";
//read_hubble();
//char* url = "/home/edoardo/Analysis/Tables/VDE/growth_factor_from_Pk.dat";
//char* url = "/home/edoardo/Analysis/Tables/LCDM/lcdm_growth_factor.dat";
//char* url = "/home/edoardo/Analysis/Tables/LCDM/growthFactor_NormalizedToA.dat";
//char* url = "/home/edoardo/Analysis/Tables/LCDM/Cosmology.DAT";
/*
read_growth_factor(url);
Settings.model_type=3;
Cosmo.OmegaL=0.612;
Cosmo.OmegaM=0.388;
Cosmo.h=0.62;
Cosmo.H_0=Cosmo.h*100;

calculate_gamma();


double rad = PI/180;
double th1 = 0*rad;
double th2 = 5*rad;
double ph1 = 0*rad;
double ph2 = 5*rad;

double z0 = 1.4; double z1=2.2;
double h, h_inv;
void *p;
h = H_z(z1, p); h_inv = inv_H_z(z1, p);

fprintf(stderr, "z: %lf, H_z: %lf, H_inv: %lf \n", z1, h, h_inv);

double r = integrate_solid_angle(th1,th2,ph1,ph2);
double dh = comoving_distance(z0,z1);
*/
//double dh1 = comoving_distance(0,z1);
//double dh2 = comoving_distance(0,z0);
//p = (void*) calloc(1,sizeof(void));
//p[0]=z1;
//double dv = integrate_comoving_volume(z0,z1); // comoving_vol(z0,z1);
//double dh1 =integrate_number_density(0,z0); // comoving_vol(z0,z1);
//double dh2 =integrate_number_density(0,z1); // comoving_vol(z0,z1);
//fprintf(stderr, "solid angle integration res: %lf c_dist: %e c_vol: %e", r, dh, dv);

/*
double a_x[] = {80,100,123,153,189,234,290,359,444,550,681,842,1043,1290};
double a_y[] = {94075,70255,56255,43573,32829,23855,17323,12085,8221,5560,3814,2674,1883,1236};
//double a_y[] = {77363,43641,40611,29402,20826,13828,10024,6232,3903,2587,1863,1401,1000,512};
double e_y[] = {10085,7285,4392,2904,1889,1223,758,487,304,183,107,63,38,25};

int nbin = 14;
int k=0;

//for(){}

NFW.bins=nbin;
NFW.radius = (double *) calloc(nbin, sizeof(double));
NFW.overd = (double *) calloc(nbin, sizeof(double));
NFW.err = (double *) calloc(nbin, sizeof(double));

for(k=0; k<nbin; k++) {
NFW.radius[k]=a_x[k];
NFW.overd[k]=a_y[k];
}

double rb = 1; //pow(512,3)*5.55*pow(10,11)*pow(10,-12);

rb *= 1; //6.03e15/(4*PI*pow(2606,3));

double rs = 328;
double cc = 2204/rs;
double rho0 = cc*cc*cc*200/(3*(log(1+cc)-(cc/(1+cc))));

best_fit_nfw(rho0, rs, nbin, a_x, a_y, e_y);
fprintf(stderr, "Chi2: %lf \n", Chi2.chi);
print_nfw();
*/

return 0;
}



