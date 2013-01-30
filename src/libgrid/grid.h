extern struct node
{
	struct halo *HALO;	

	int n_haloes;

	double M;
	double M_CIC;
	double X[3];
	double X_cm[3];

} *Node;


extern struct grid
{
	int N;
	int N3;
	
	double L;
	double grid_size;
	double half_grid;	
	double cell_volume;

} Grid;


void init_grid(int);

void fill_grid_NGP(void);
void fill_grid_CIC(void);

void NGP_assignment(int);
void CIC_assignment(int);

void find_nodes_CIC(int*, int);
void set_fac_CIC(double*, double*);

void free_nodes(void);
