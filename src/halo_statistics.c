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

		MPI_Finalize();

		//free_comm_structures();
#endif
	
	if(ThisTask==0)
	{
	
		omp_set_num_threads(OMP_THREADS);	
	
		HALO_INDEX=0;
		MF_INDEX=0;
		PK_INDEX=0;

		initialize_halo_properties_structure();

		compute_halo_properties();
#ifdef GAS
		average_gas_profiles();
#endif
		average_nfw_profile();

	//	sort_axis_alignment();
	//	print_axis_alignment();

		print_average_profiles();
		print_halo_best_fit_results();
		print_numerical_mass_function();
		print_all_halo_properties_to_one_file();
		//print_all_haloes();

		read_v_web();
		read_t_web();

		sort_web_statistics();
	
		assign_haloes_to_web();

			// Assign halo to nodes
			int i=0;
			char base_out[200]; 
			sprintf(base_out, "%s", Urls.output_prefix);
			Settings.use_web = 1;

			for(i=0; i<4; i++)
			{
				Settings.use_web_type = i;
				sprintf(Urls.output_prefix, "%s%s_%02d_", base_out, "type", i);			
				fprintf(stderr, "Sorting halo properties for web type=%d, using %d haloes.\n", 
						i, n_haloes_per_criterion());				
				//compute_halo_properties();
				average_nfw_profile();
#ifdef GAS
				average_gas_profiles();
#endif
	
				// Print files
				print_average_profiles();
				print_halo_best_fit_results();
				print_numerical_mass_function();
				print_all_halo_properties_to_one_file();
			}

		// Now select subhaloes
		find_substructure();

		Settings.use_web = 0;
		Settings.use_sub = 1;
		sprintf(Urls.output_prefix, "%s%s", base_out, "substructure_");			
		compute_subhalo_properties();
		//compute_halo_properties();
		average_nfw_profile();
#ifdef GAS
		average_gas_profiles();
#endif
	
		print_average_profiles();
		print_halo_best_fit_results();
		print_numerical_mass_function();
		print_all_halo_properties_to_one_file();

		print_subhalo_only_properties();

	//	free_halo_properties();

/*
*/
	INFO_MSG("Computed halo statistical properties at fixed z");

	}

	//MPI_Barrier(MPI_COMM_WORLD);

	//if(ThisTask == 0)
//#endif

//#ifdef WITH_MPI
//	MPI_Finalize();
//#endif

return 0;
}

