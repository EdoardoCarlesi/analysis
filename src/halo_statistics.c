#include <stdio.h>
#include <stdlib.h>

#include "libio/write_io.h"
#include "libio/halo_io.h"
#include "libhalo/halo_properties.h"
#include "general_functions.h"
#include "general_variables.h"

#ifdef WITH_MPI
#include <mpi.h>
#include "libparallel/general.h"
#endif


int main(int argc, char **argv)
{
	fprintf(stderr, "\nComputing halo statistical properties.\n");

	initialize_internal_variables(argv);

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
	MPI_Comm_size(MPI_COMM_WORLD, &NTask);

	init_comm_structures();
	init_cpu_struct();	
		
	generate_url_for_tasks();	
#endif	

			get_halo_files_urls();
			
			use_halo_url(0);

			read_halo_file();

#ifdef WITH_MPI
		MPI_Barrier(MPI_COMM_WORLD);

		gather_halo_structures();

		free_comm_structures();
	
	if(ThisTask==0)
	{
#endif
	
		initialize_halo_properties_structure();

		compute_halo_properties();

		print_all_halo_properties_to_one_file();
	
#ifdef WITH_MPI
	}

	MPI_Barrier(MPI_COMM_WORLD);

		if(ThisTask == 0)
#endif
			fprintf(stderr, "\nHalo statistical properties computed.\n");

#ifdef WITH_MPI
	MPI_Finalize();
#endif

return 0;
}

