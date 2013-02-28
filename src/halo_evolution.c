#include <stdio.h>
#include <stdlib.h>

#include "libio/io.h"
#include "libhalo/halo.h"

#include "general_def.h"


int main(int argc, char **argv)
{
	int i=0; 

	INFO_MSG("Computing halo properties as a function of redshift");
	
	initialize_internal_variables(argv);

		read_redshift_file();

		get_halo_files_urls();

		initialize_halo_storage();

			for (i=0; i<Urls.nCatalogueFiles-1; i++) 
				compute_halo_and_subhalo_statistics(i);

		print_evolution_to_file();

	INFO_MSG("Done halo properties computation");

return 0;
}
