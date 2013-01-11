#include <mpi.h>

extern int ThisTask;
extern int NTask;

extern int *SizeDisplStructHalo;
extern int *SizeHaloesStructHalo;

extern int *SizeDispl;
extern int *SizeHaloes;
extern int *SizeDisplThreshold;
extern int *SizeHaloesThreshold;

void gather_halo_structures(void);
void init_comm_structures(void);
void free_comm_structures(void);
void copy_halo_url(char*);
