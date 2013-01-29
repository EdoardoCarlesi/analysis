
extern struct node
{
	struct halo *HALO;	

	int n_haloes;

	double M;
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

} Grid;


void init_grid(int);

void fill_grid(void);

void assign_halo_to_node(int);
