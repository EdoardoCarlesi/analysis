#include <mpi.h>

extern int ThisTask;
extern int NTask;

extern int *SizeDisplStructHalo;
extern int *SizeHaloesStructHalo;

extern int *SizeDispl;
extern int *SizeHaloes;
extern int *SizeDisplThreshold;
extern int *SizeHaloesThreshold;

void mpi_bcast_size();

void set_halo_displacement(void);
void init_pstructures(void);
void copy_url(char*);
