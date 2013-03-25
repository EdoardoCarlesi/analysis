#include <stdio.h>
#include <stdlib.h>

#include "libio/io.h"
#include "libhalo/halo.h"
#include "libcosmo/cosmo.h"

#include "general_def.h"

#ifdef WITH_MPI
#include <mpi.h>
#include "libparallel/general.h"
#endif


int main(int argc, char **argv)
{
	INFO_MSG("Computing halo statistical properties at fixed z");

	initialize_internal_variables(argv);
	initialize_halo_storage();

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
	MPI_Comm_size(MPI_COMM_WORLD, &NTask);

	init_comm_structures();
	init_cpu_struct();	

	MPI_Bcast(&Settings, sizeof(struct general_settings), MPI_BYTE, 0, MPI_COMM_WORLD);
		
	generate_url_for_tasks();	
#endif	
			get_halo_files_urls();
			
			set_halo_url();

			read_halo_file();

			read_profiles_file();
	
			fit_and_store_nfw_parameters();
			
#ifdef GAS
			fit_and_store_gas_parameters();
#endif
			free_halo_profiles();

#ifdef WITH_MPI

		MPI_Barrier(MPI_COMM_WORLD);

		gather_halo_structures();

	//	free_comm_structures();
	
	if(ThisTask==0)
	{
#endif
	
		omp_set_num_threads(OMP_THREADS);	
	
		HALO_INDEX=0;
		MF_INDEX=0;
		PK_INDEX=0;

		initialize_halo_properties_structure();

		find_substructure();

		compute_halo_properties();
		print_numerical_mass_function();
		print_all_halo_properties_to_one_file();
		print_all_haloes();

#ifdef GAS
		average_gas_profiles();
#endif
		average_nfw_profile();
		print_average_profiles();

		print_halo_best_fit_results();

	//	sort_axis_alignment();
	//	print_axis_alignment();

	//	compute_subhalo_properties();
	//	print_all_subhalo_properties_to_one_file();
	// Redo all of the above for subhaloes

		// FIXME
		free_halo_properties();
	
#ifdef WITH_MPI
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(ThisTask == 0)
#endif
	INFO_MSG("Computed halo statistical properties at fixed z");

#ifdef WITH_MPI
	MPI_Finalize();
#endif

return 0;
}

