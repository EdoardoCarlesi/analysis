#include <stdio.h>

#include "libio/write_io.h"
#include "libio/halo_io.h"
#include "libio/power_io.h"
#include "libcosmo/tinker.h"
#include "libcosmo/mass_function.h"
#include "libcosmo/cosmological_relations.h"
#include "general_functions.h"
#include "general_variables.h"

#ifdef WITH_MPI
#include <mpi.h>
#include "libparallel/general.h"
#endif


int main(int arg, char **argv)
{
	fprintf(stdout, "\nComputing the mass function, starting the main().\n");

	initialize_internal_variables(argv);

#ifdef WITH_MPI
	MPI_Init(&arg, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
	MPI_Comm_size(MPI_COMM_WORLD, &NTask);

	init_comm_structures();
	init_cpu_struct();	
		
	generate_url_for_tasks();	
#endif

	Settings.rho_0 = default_rho0();
	Settings.use_one_pk = 1;
		// First guess parameters and fits for the different mass functions
	T_mf.A = 0.186; 
	T_mf.a = 1.47; 
	T_mf.b = 2.57; 
	T_mf.c = 1.19;

#ifdef WITH_MPI
	if(ThisTask == 0)
#endif
		init_pks();

		get_halo_files_urls();

		set_halo_url();

		read_halo_file();

#ifdef WITH_MPI
		MPI_Barrier(MPI_COMM_WORLD);

		gather_halo_structures();

		free_comm_structures();
	
	if(ThisTask==0)
	{
#endif
		compute_numerical_mass_function();
		compute_theoretical_mass_function();
		print_numerical_mass_function();	
		print_theoretical_mass_function();	
	
#ifdef WITH_MPI
	}

	MPI_Barrier(MPI_COMM_WORLD);

		if(ThisTask == 0)
#endif

	fprintf(stdout, "\nMass function computed.\n");

#ifdef WITH_MPI
	MPI_Finalize();
#endif
	return 0;
}
