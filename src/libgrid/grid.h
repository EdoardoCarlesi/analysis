#include "../general_variables.h"

extern struct node
{

	struct halo *HALO;	
	
	int ID; // Node number in 1 d
	int n_haloes;

	double M;
	double X;
	double Y;
	double Z;

} *Node;


void init_grid(int);

void assign_halo_to_node(int);
