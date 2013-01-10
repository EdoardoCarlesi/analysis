#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <malloc.h>
#include <mpi.h>
#include "general.h"


#include "../general_variables.h"


int ThisTask;
int NTask;
int NHaloes;
int TaskHaloes;

MPI_Datatype 	MPI_halo;


void copy_url(char *url)
{
	pUrls[ThisTask].halo_file = (char*) calloc(strlen(url)-1, sizeof(char));
	strcpy(pUrls[ThisTask].halo_file, url);
}


void init_pstructures()
{
	pSettings = (struct general_settings *) calloc(NTask, sizeof(struct general_settings));
	pHaloes = (struct halo **) calloc(NTask, sizeof(struct halo *));
	pUrls = (struct internal_urls *) calloc(NTask, sizeof(struct internal_urls));
	pFC = (struct full_catalogue *) calloc(NTask, sizeof(struct full_catalogue));
}



void create_mpi_halo_datatype()
{
	// Datatype tipe size:
	// 9 int
	// 1 int* 
	// 37 double
	// 7 double*
	// GAS
	// 2 int
	// 20 double
	// EXTRA GAS
	// 12 double
	int 		j=0;
	int 		items=88;
	int 		blocklen[88];

	MPI_Aint	disp[88];
	MPI_Datatype 	types[88] = 
	{
	MPI_INT,	MPI_INT,	MPI_INT,	MPI_INT,	MPI_INT,	
	MPI_INT,	MPI_INT,	MPI_INT,	MPI_INT,	
	MPI_INT,
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	 	
	MPI_INT,	MPI_INT,
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	MPI_DOUBLE, 	
	MPI_DOUBLE, 	MPI_DOUBLE, 	 	
	};

	for(j=0; j<items; j++)
		blocklen[j] = 1;
	
	disp[0] = offsetof(struct halo, n_part);
	disp[0] = offsetof(struct halo, n_satellites);
        disp[0] = offsetof(struct halo, n_bins);
	disp[0] = offsetof(struct halo, neg_r_bins);
	disp[0] = offsetof(struct halo, id);
	disp[0] = offsetof(struct halo, host);
	disp[0] = offsetof(struct halo, virial);
	disp[0] = offsetof(struct halo, spin);
	disp[0] = offsetof(struct halo, conc);

	disp[0] = offsetof(struct halo, id_satellites);

	disp[0] = offsetof(struct halo, Xc);
	disp[0] = offsetof(struct halo, Yc);
	disp[0] = offsetof(struct halo, Zc);
	disp[0] = offsetof(struct halo, VXc);
	disp[0] = offsetof(struct halo, VYc);
	disp[0] = offsetof(struct halo, VZc);
	disp[0] = offsetof(struct halo, Lx);
	disp[0] = offsetof(struct halo, Ly);
	disp[0] = offsetof(struct halo, Lz);
	disp[0] = offsetof(struct halo, Eax); 
	disp[0] = offsetof(struct halo, Eay); 
	disp[0] = offsetof(struct halo, Eaz); 
	disp[0] = offsetof(struct halo, Mvir);
	disp[0] = offsetof(struct halo, Rvir);
	disp[0] = offsetof(struct halo, Rmax);
	disp[0] = offsetof(struct halo, Vmax);
	disp[0] = offsetof(struct halo, Ekin);
	disp[0] = offsetof(struct halo, Epot);
	disp[0] = offsetof(struct halo, th_vir);
	disp[0] = offsetof(struct halo, AngMom);
	disp[0] = offsetof(struct halo, Jcirc);
	disp[0] = offsetof(struct halo, ecc);
	disp[0] = offsetof(struct halo, lambda);
	disp[0] = offsetof(struct halo, lambdaE);
	disp[0] = offsetof(struct halo, delta_c);
	disp[0] = offsetof(struct halo, c);
	disp[0] = offsetof(struct halo, r2);
	disp[0] = offsetof(struct halo, aa);
	disp[0] = offsetof(struct halo, bb);
	disp[0] = offsetof(struct halo, cc);
	disp[0] = offsetof(struct halo, c_a);
	disp[0] = offsetof(struct halo, triax);
	disp[0] = offsetof(struct halo, z_form);
	disp[0] = offsetof(struct halo, c_nfw);
	disp[0] = offsetof(struct halo, chi_nfw);
	disp[0] = offsetof(struct halo, rs_nfw);
	disp[0] = offsetof(struct halo, rho0_nfw);

	disp[0] = offsetof(struct halo, radius);
	disp[0] = offsetof(struct halo, rho);
	disp[0] = offsetof(struct halo, over_rho);
	disp[0] = offsetof(struct halo, err);
	disp[0] = offsetof(struct halo, over_err);
	disp[0] = offsetof(struct halo, err_dn);
	disp[0] = offsetof(struct halo, bin);

#ifdef GAS
	disp[0] = offsetof(struct halo, N_dm);
	disp[0] = offsetof(struct halo, N_gas);
	disp[0] = offsetof(struct halo, M_dm);
	disp[0] = offsetof(struct halo, M_gas);
	disp[0] = offsetof(struct halo, Eax_gas);
	disp[0] = offsetof(struct halo, Eay_gas);
	disp[0] = offsetof(struct halo, Eaz_gas);
	disp[0] = offsetof(struct halo, Ebx_gas);
	disp[0] = offsetof(struct halo, Eby_gas);
	disp[0] = offsetof(struct halo, Ebz_gas);
	disp[0] = offsetof(struct halo, Ecx_gas);
	disp[0] = offsetof(struct halo, Ecy_gas);
	disp[0] = offsetof(struct halo, Ecz_gas);
	disp[0] = offsetof(struct halo, lambda_gas);
	disp[0] = offsetof(struct halo, lambdaE_gas);
	disp[0] = offsetof(struct halo, b_gas);
	disp[0] = offsetof(struct halo, c_gas);
	disp[0] = offsetof(struct halo, Ekin_gas);
	disp[0] = offsetof(struct halo, Epot_gas);
	disp[0] = offsetof(struct halo, b_fraction);
	disp[0] = offsetof(struct halo, Cum_u_gas);
	disp[0] = offsetof(struct halo, T_gas);
#ifdef EXTRA_GAS
	disp[0] = offsetof(struct halo, X_dm);
	disp[0] = offsetof(struct halo, Y_dm);
	disp[0] = offsetof(struct halo, Z_dm);
	disp[0] = offsetof(struct halo, VX_dm);
	disp[0] = offsetof(struct halo, VY_dm);
	disp[0] = offsetof(struct halo, VZ_dm);
	disp[0] = offsetof(struct halo, X_gas);
	disp[0] = offsetof(struct halo, Y_gas);
	disp[0] = offsetof(struct halo, Z_gas);
	disp[0] = offsetof(struct halo, VX_gas);
	disp[0] = offsetof(struct halo, VY_gas);
	disp[0] = offsetof(struct halo, VZ_gas);
#endif
#endif

/*
	disp[0] = &pHaloes[ThisTask][0].n_part - pHaloes[ThisTask][0];
	disp[1] = &pHaloes[ThisTask][0].n_satellites - &pHaloes[ThisTask][0];
	disp[2] = &pHaloes[ThisTask][0].n_bins - &pHaloes[ThisTask][0];
	disp[3] = &pHaloes[ThisTask][0].neg_r_bins - &pHaloes[ThisTask][0];
	disp[4] = &pHaloes[ThisTask][0].id - &pHaloes[ThisTask][0];
	disp[5] = &pHaloes[ThisTask][0].host - &pHaloes[ThisTask][0];
	disp[6] = &pHaloes[ThisTask][0].virial - &pHaloes[ThisTask][0];
	disp[7] = &pHaloes[ThisTask][0].spin - &pHaloes[ThisTask][0];
	disp[8] = &pHaloes[ThisTask][0].conc - &pHaloes[ThisTask][0];

	disp[9] = &pHaloes[ThisTask][0].id_satellites - &pHaloes[ThisTask][0];

	disp[10] = &pHaloes[ThisTask][0].Xc - &pHaloes[ThisTask][0];
	disp[11] = &pHaloes[ThisTask][0].Yc - &pHaloes[ThisTask][0];
	disp[12] = &pHaloes[ThisTask][0].Zc - &pHaloes[ThisTask][0];
	disp[13] = &pHaloes[ThisTask][0].VXc - &pHaloes[ThisTask][0];
	disp[14] = &pHaloes[ThisTask][0].VYc - &pHaloes[ThisTask][0];
	disp[15] = &pHaloes[ThisTask][0].VZc - &pHaloes[ThisTask][0];
	disp[16] = &pHaloes[ThisTask][0].Lx - &pHaloes[ThisTask][0];
	disp[17] = &pHaloes[ThisTask][0].Ly - &pHaloes[ThisTask][0];
	disp[18] = &pHaloes[ThisTask][0].Lz - &pHaloes[ThisTask][0];
	disp[19] = &pHaloes[ThisTask][0].Eax - &pHaloes[ThisTask][0];
	disp[20] = &pHaloes[ThisTask][0].Eay - &pHaloes[ThisTask][0];
	disp[21] = &pHaloes[ThisTask][0].Eaz - &pHaloes[ThisTask][0];
	disp[22] = &pHaloes[ThisTask][0].Mvir - &pHaloes[ThisTask][0];
	disp[23] = &pHaloes[ThisTask][0].Rvir - &pHaloes[ThisTask][0];
	disp[24] = &pHaloes[ThisTask][0].Rmax - &pHaloes[ThisTask][0];
	disp[25] = &pHaloes[ThisTask][0].Vmax - &pHaloes[ThisTask][0];
	disp[26] = &pHaloes[ThisTask][0].Ekin - &pHaloes[ThisTask][0];
	disp[27] = &pHaloes[ThisTask][0].Epot - &pHaloes[ThisTask][0];
	disp[28] = &pHaloes[ThisTask][0].th_vir - &pHaloes[ThisTask][0];
	disp[29] = &pHaloes[ThisTask][0].AngMom - &pHaloes[ThisTask][0];
	disp[30] = &pHaloes[ThisTask][0].Jcirc - &pHaloes[ThisTask][0];
	disp[31] = &pHaloes[ThisTask][0].ecc - &pHaloes[ThisTask][0];
	disp[32] = &pHaloes[ThisTask][0].lambda - &pHaloes[ThisTask][0];
	disp[33] = &pHaloes[ThisTask][0].lambdaE - &pHaloes[ThisTask][0];
	disp[34] = &pHaloes[ThisTask][0].delta_c - &pHaloes[ThisTask][0];
	disp[35] = &pHaloes[ThisTask][0].c - &pHaloes[ThisTask][0];
	disp[36] = &pHaloes[ThisTask][0].r2 - &pHaloes[ThisTask][0];
	disp[37] = &pHaloes[ThisTask][0].aa - &pHaloes[ThisTask][0];
	disp[38] = &pHaloes[ThisTask][0].bb - &pHaloes[ThisTask][0];
	disp[39] = &pHaloes[ThisTask][0].cc - &pHaloes[ThisTask][0];
	disp[40] = &pHaloes[ThisTask][0].c_a - &pHaloes[ThisTask][0];
	disp[41] = &pHaloes[ThisTask][0].triax - &pHaloes[ThisTask][0];
	disp[42] = &pHaloes[ThisTask][0].z_form - &pHaloes[ThisTask][0];
	disp[43] = &pHaloes[ThisTask][0].c_nfw - &pHaloes[ThisTask][0];
	disp[44] = &pHaloes[ThisTask][0].chi_nfw - &pHaloes[ThisTask][0];
	disp[45] = &pHaloes[ThisTask][0].rs_nfw - &pHaloes[ThisTask][0];
	disp[46] = &pHaloes[ThisTask][0].rho0_nfw - &pHaloes[ThisTask][0];

	disp[47] = &pHaloes[ThisTask][0].radius - &pHaloes[ThisTask][0];
	disp[48] = &pHaloes[ThisTask][0].rho - &pHaloes[ThisTask][0];
	disp[49] = &pHaloes[ThisTask][0].over_rho - &pHaloes[ThisTask][0];
	disp[50] = &pHaloes[ThisTask][0].err - &pHaloes[ThisTask][0];
	disp[51] = &pHaloes[ThisTask][0].over_err - &pHaloes[ThisTask][0];
	disp[52] = &pHaloes[ThisTask][0].err_dn - &pHaloes[ThisTask][0];
	disp[53] = &pHaloes[ThisTask][0].bin - &pHaloes[ThisTask][0];

#ifdef GAS
	disp[54] = &pHaloes[ThisTask][0].N_dm - &pHaloes[ThisTask][0];
	disp[55] = &pHaloes[ThisTask][0].N_gas - &pHaloes[ThisTask][0];
	disp[56] = &pHaloes[ThisTask][0].M_dm - &pHaloes[ThisTask][0];
	disp[57] = &pHaloes[ThisTask][0].M_gas - &pHaloes[ThisTask][0];
	disp[58] = &pHaloes[ThisTask][0].Eax_gas - &pHaloes[ThisTask][0];
	disp[59] = &pHaloes[ThisTask][0].Eay_gas - &pHaloes[ThisTask][0];
	disp[60] = &pHaloes[ThisTask][0].Eaz_gas - &pHaloes[ThisTask][0];
	disp[61] = &pHaloes[ThisTask][0].Ebx_gas - &pHaloes[ThisTask][0];
	disp[62] = &pHaloes[ThisTask][0].Eby_gas - &pHaloes[ThisTask][0];
	disp[63] = &pHaloes[ThisTask][0].Ebz_gas - &pHaloes[ThisTask][0];
	disp[64] = &pHaloes[ThisTask][0].Ecx_gas - &pHaloes[ThisTask][0];
	disp[65] = &pHaloes[ThisTask][0].Ecy_gas - &pHaloes[ThisTask][0];
	disp[66] = &pHaloes[ThisTask][0].Ecz_gas - &pHaloes[ThisTask][0];
	disp[67] = &pHaloes[ThisTask][0].lambda_gas - &pHaloes[ThisTask][0];
	disp[68] = &pHaloes[ThisTask][0].lambdaE_gas - &pHaloes[ThisTask][0];
	disp[69] = &pHaloes[ThisTask][0].b_gas - &pHaloes[ThisTask][0];
	disp[70] = &pHaloes[ThisTask][0].c_gas - &pHaloes[ThisTask][0];
	disp[71] = &pHaloes[ThisTask][0].Ekin_gas - &pHaloes[ThisTask][0];
	disp[72] = &pHaloes[ThisTask][0].Epot_gas - &pHaloes[ThisTask][0];
	disp[73] = &pHaloes[ThisTask][0].b_fraction - &pHaloes[ThisTask][0];
	disp[74] = &pHaloes[ThisTask][0].Cum_u_gas - &pHaloes[ThisTask][0];
	disp[75] = &pHaloes[ThisTask][0].T_gas - &pHaloes[ThisTask][0];
#ifdef EXTRA_GAS
	disp[76] = &pHaloes[ThisTask][0].X_dm - &pHaloes[ThisTask][0];
	disp[77] = &pHaloes[ThisTask][0].Y_dm - &pHaloes[ThisTask][0];
	disp[78] = &pHaloes[ThisTask][0].Z_dm - &pHaloes[ThisTask][0];
	disp[79] = &pHaloes[ThisTask][0].VX_dm - &pHaloes[ThisTask][0];
	disp[80] = &pHaloes[ThisTask][0].VY_dm - &pHaloes[ThisTask][0];
	disp[81] = &pHaloes[ThisTask][0].VZ_dm - &pHaloes[ThisTask][0];
	disp[82] = &pHaloes[ThisTask][0].X_gas - &pHaloes[ThisTask][0];
	disp[83] = &pHaloes[ThisTask][0].Y_gas - &pHaloes[ThisTask][0];
	disp[84] = &pHaloes[ThisTask][0].Z_gas - &pHaloes[ThisTask][0];
	disp[85] = &pHaloes[ThisTask][0].VX_gas - &pHaloes[ThisTask][0];
	disp[86] = &pHaloes[ThisTask][0].VY_gas - &pHaloes[ThisTask][0];
	disp[87] = &pHaloes[ThisTask][0].VZ_gas - &pHaloes[ThisTask][0];
#endif
#endif
*/
	MPI_Type_create_struct(items,blocklen,disp,types,&MPI_halo);
	MPI_Type_commit(&MPI_halo);
}
