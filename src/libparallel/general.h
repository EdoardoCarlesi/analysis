extern int ThisTask;
extern int NTask;

extern struct Cpu
{
  char name[5];

} *cpu;

extern int *SizeDisplStructHalo;
extern int *SizeHaloesStructHalo;

void init_cpu_struct(void);
void generate_url_for_tasks(void);
void gather_halo_structures(void);
void init_comm_structures(void);
void free_comm_structures(void);
