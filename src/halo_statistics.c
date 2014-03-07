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
	initialize_internal_variables(argv);
	initialize_halo_storage();

	int ifile=0, locFileNumber=0, totFilesPerTask=0, filesPerTask=0, extraFilesPerTask=0;
	char base_out[200]; 
	sprintf(base_out, "%s", Urls.output_prefix);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
	MPI_Comm_size(MPI_COMM_WORLD, &NTask);

	init_comm_structures();
	init_cpu_struct();	

	MPI_Bcast(&Settings, sizeof(struct general_settings), MPI_BYTE, 0, MPI_COMM_WORLD);
		
        if(ThisTask == 0)
	  fprintf(stderr, "Total Number Of Files = %d, Total Number Of Tasks = %d\n", Settings.tot_files, NTask);

	filesPerTask = Settings.tot_files / NTask;
	extraFilesPerTask = Settings.tot_files % NTask;

	if(ThisTask < extraFilesPerTask)
	  totFilesPerTask = filesPerTask + 1;
	else
	  totFilesPerTask = filesPerTask;

	fprintf(stderr, "Task=%d will read %d files\n", ThisTask, totFilesPerTask);

	// Loop over files' chunks
	for(ifile = 0; ifile < totFilesPerTask; ifile++)
        {
		locFileNumber = ThisTask + ifile * NTask;

		generate_url_for_tasks(locFileNumber);	

		get_halo_files_urls();
		
		set_halo_url();

		read_halo_file();

#ifndef NO_PROFILES
		read_profiles_file();

		fit_and_store_nfw_parameters();
#ifdef GAS
		fit_and_store_gas_parameters();
#endif
		free_halo_profiles();

		MPI_Barrier(MPI_COMM_WORLD);
#endif
		realloc_haloes();
	}
	
	pHaloes[ThisTask] = tempHaloes;
	pSettings[ThisTask].n_haloes = pSettings[ThisTask].n_haloes_step;

	gather_halo_structures();

	MPI_Finalize();

	free_comm_structures();

	// Files have been read in, now start the analysis
	if(ThisTask == 0)
	{
		HALO_INDEX=0;
		MF_INDEX=0;
		PK_INDEX=0;

		omp_set_num_threads(OMP_THREADS);	
	
		initialize_halo_properties_structure();

#ifndef WEB_ONLY
		compute_halo_properties();
		find_substructure();

#ifndef NO_PROFILES
		average_nfw_profile();
		print_average_profiles();
#ifdef GAS
		average_gas_profiles();
#endif

#endif

#ifdef COMPUTE_ALIGN
		sort_axis_alignment();
		print_axis_alignment();
#endif

#ifndef USE_N_HALOES
		print_numerical_mass_function();
#endif
		print_all_halo_properties_to_one_file();
		
		print_halo_best_fit_results();
#endif /* WEB_ONLY */

		/* If NO_WEB is not selected, the above analysis is repeated for every environment */ 

#ifndef NO_WEB
		read_v_web();
#ifndef WEB_ONLY
		read_t_web();
		sort_web_statistics();
		print_web_statistics();
#endif
		Settings.use_web = 0;
		assign_haloes_to_web();
#ifdef WEB_ONLY
		// This prints all haloes with informations on the region they belong to
		print_all_haloes();
#endif

#endif
			// Assign halo to nodes
			int i=0;

#ifndef WEB_ONLY
			Settings.use_sub = 0;
			Settings.use_web = 1;
#ifndef NO_WEB
			for(i=0; i<4; i++)
			{
				Settings.use_web_type = i;
				sprintf(Urls.output_prefix, "%s%s_%02d_", base_out, "type", i);			
				fprintf(stderr, "Sorting halo properties for web type=%d, using %d haloes.\n", 
						i, n_haloes_per_criterion());				
				compute_halo_properties();
#ifndef NO_PROFILES
				average_nfw_profile();
#ifdef GAS
				average_gas_profiles();
#endif
				// Print files
				print_average_profiles();
				print_halo_best_fit_results();
#endif
				print_numerical_mass_function();
				print_all_halo_properties_to_one_file();
				fprintf(stderr, "Done sorting properties for type %d.\n", i);
			}

		Settings.use_web = 0;
#endif /* NO_WEB */

#ifdef CROSS_CORRELATION
		cross_correlation("/home/carlesi/Analysis/output/CC/lcdm-250-1024-z0.000_tot_haloes.dat");
#endif

#ifdef SUBHALO
		// Now select subhaloes
	//	find_substructure();

		Settings.use_sub = 1;
		sprintf(Urls.output_prefix, "%s%s", base_out, "substructure_");			
		compute_subhalo_properties();
		compute_halo_properties();
		average_nfw_profile();
#ifdef GAS
		average_gas_profiles();
#endif

#ifndef NO_PROFILES
		print_average_profiles();
#endif
		print_halo_best_fit_results();
		print_numerical_mass_function();
		print_all_halo_properties_to_one_file();
		print_subhalo_only_properties();
		free_halo_properties();
#endif /* SUBHALO */

#endif /* WEB ONLY */
		free(Haloes);

	  INFO_MSG("Computed halo statistical properties at fixed z");
	}

  return 0;
}

