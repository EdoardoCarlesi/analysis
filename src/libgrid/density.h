extern struct density
{
	int N;

	double *r_n;
	double *rho;
	double *mass;

} Density;



void init_density(void);

void free_density(void);

void find_density_maxima_and_minima(void);
