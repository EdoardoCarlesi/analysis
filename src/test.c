/* Test new routines and functions */
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "libio/write_io.h"
#include "libio/halo_io.h"
#include "libhalo/halo_properties.h"
#include "general_variables.h"

#include "libgrid/grid.h"
#include "libgrid/density.h"

#ifdef WITH_MPI
#include <mpi.h>
#include "libparallel/general.h"
#endif


int main(int argc, char **argv)
{
	char *halo_urls[2];

	halo_urls[0] = "/home/edoardo/snapshot_029.0000.z0.000.AHF_halos";
	halo_urls[1] = "/home/edoardo/snapshot_029.0001.z0.000.AHF_halos";

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
	MPI_Comm_size(MPI_COMM_WORLD, &NTask);

	init_comm_structures();

	init_cpu_struct();
	
	if(ThisTask == 0)
	{	// Example settings
		Settings.mass_min = 3.1e10;
		Settings.n_part_1D = 256;
		Settings.box_size = 75;
		Settings.n_min = 1e2;
		Settings.n_bins = 10;
		Cosmo.virial = 1.5;
		Cosmo.spin = 0.15;
	}
			// Broadcast settings to all Tasks
		MPI_Bcast(&Settings, sizeof(struct general_settings), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&Cosmo, sizeof(struct cosmology), MPI_BYTE, 0, MPI_COMM_WORLD);

		pUrls[ThisTask].halo_file = (char *) calloc(strlen(halo_urls[ThisTask]), sizeof(char)); 
		strcpy(pUrls[ThisTask].halo_file, halo_urls[ThisTask]);

		MPI_Barrier(MPI_COMM_WORLD);

		read_halo_file();

		// Gather all halo structures from every processor into one
	gather_halo_structures();
		
		// Get rid of the unused structures
	free_comm_structures();

	if(ThisTask==10)
	{
		// Do your (serial / OpenMP) operations	
		init_grid(128);
		fill_grid_CIC();
		print_grid_CIC();
		init_density();
		find_density_maxima_and_minima();
	}

	MPI_Finalize();

	return 0;
}
