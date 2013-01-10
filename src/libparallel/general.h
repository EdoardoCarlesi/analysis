#include <mpi.h>

extern int ThisTask;
extern int NTask;
extern int TaskHaloes;
extern int NHaloes;

extern MPI_Datatype MPI_halo;

void init_pstructures();
void copy_url(char*);
void create_mpi_halo_datatype();
