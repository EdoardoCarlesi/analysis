#include <stdio.h>
#include <stdlib.h>

#include "libio/write_io.h"
#include "libio/halo_io.h"
#include "libhalo/halo_properties.h"

#include "general_functions.h"
#include "general_variables.h"

#ifdef WITH_MPI
#include <mpi.h>
#include "libparallel/halo_io.h"
#include "libparallel/general.h"
#endif


int main(int argc, char **argv)
{
	fprintf(stderr, "\nmain(). Computing halo statistical properties.\n");

		initialize_internal_variables(argv);

		initialize_halo_properties_structure();

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
	MPI_Comm_size(MPI_COMM_WORLD, &NTask);

	init_comm_structures();
			
	generate_urls_for_tasks();	
			// TODO init halo urls differently
		//copy_halo_url(halo_urls[ThisTask]);

		MPI_Barrier(MPI_COMM_WORLD);

		mpi_read_halo_file();

		gather_halo_structures();

		free_comm_structures();
	
	if(ThisTask==0)
	{
#endif

			// TODO Check the number of haloes to use!!!!!!!!
		Settings.n_haloes_to_use=39;

		read_halo_file();

		compute_halo_properties();

		print_all_halo_properties_to_one_file();
	
#ifdef WITH_MPI
	}

	MPI_Finalize();
#endif

	fprintf(stderr, "\nmain(). Halo statistical properties computed.\n");

return 0;
}

