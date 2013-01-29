#define WARNING(str1, str2) fprintf(stderr, "\n\t\tWARNING! %s: %s.\n", str1, str2)

#define ERROR(str1, str2) fprintf(stderr, "\n\t\tERROR! %s: %s.\n", str1, str2)

#define INFO_MSG(str) fprintf(stdout, "\n\t%s.\n", str)

void initialize_internal_variables(char **);

void default_init(void);

void set_halo_selection_criterion(void);
void check_condition_consistency(void);

int n_haloes_per_criterion(void);
int halo_condition(int);
