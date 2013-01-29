#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../general_variables.h"
#include "grid.h"

#ifdef WITH_MPI
#include <mpi.h>
#include "../libparallel/general.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


struct node *Node;


struct grid Grid;



void init_grid(int N)
{
	int i=0, j=0, k=0, n=0; 
	int N3;

	N3 = N*N*N;

	Grid.N = N;
	Grid.N3 = N3;
	Grid.L = Settings.box_size;
	Grid.grid_size = Settings.box_size / N;
	Grid.half_grid = 0.5 * Grid.grid_size; 

	Node = (struct node *) calloc(N3, sizeof(struct node));

#ifdef _OPENMP
		omp_set_num_threads(OMP_THREADS);
#endif

#		pragma omp parallel for 		\
		shared(Node) private(i,j,k)
		for(i=0; i<N; i++)
		{
			for(j=0; j<N; j++)
			{			
				for(k=0; k<N; k++)
				{
					Node[n].M = 0;
					Node[n].n_haloes = 0;
					Node[n].X[0] = Grid.half_grid + i * Grid.grid_size;
					Node[n].X[1] = Grid.half_grid + j * Grid.grid_size;
					Node[n].X[2] = Grid.half_grid + k * Grid.grid_size;
//	fprintf(stdout, "i=%d, j=%d, k=%d; Node[%d], x=%lf, y=%lf, z=%lf\n",
//			i, j, k, n, Node[n].X[0],  Node[n].X[1],  Node[n].X[2]);
					n++;		
				}
			}
		}

//	fprintf(stdout, "\nInitialized %d nodes, grid size=%lf, grid volume=%lf\n",
//		Grid.N3, Grid.grid_size, Grid.grid_size * Grid.grid_size * Grid.grid_size);

}



void assign_halo_to_node(int n)
{
	int N=0, NodeNumber=0;
	int i, j, k;
	double x, y, z;
	
	N = Grid.N;

		x = Haloes[n].Xc;
		y = Haloes[n].Yc;
		z = Haloes[n].Zc;
	
		i = (int) lrint(abs(x - Grid.half_grid)/ Grid.grid_size);
		j = (int) lrint(abs(y - Grid.half_grid)/ Grid.grid_size);
		k = (int) lrint(abs(z - Grid.half_grid)/ Grid.grid_size);
	
//	fprintf(stdout, "x/=%lf, y/=%lf, z/=%lf, rounded to i=%d, j=%d, k=%d\n",
//			x/Grid.grid_size, y/Grid.grid_size, z/Grid.grid_size, 
//			i, j, k);

	NodeNumber = N*N*i + N*j + k; 

	Node[NodeNumber].M += Haloes[n].Mvir;
	Node[NodeNumber].n_haloes++;
	
	

//	fprintf(stdout, "Halo x=%lf, y=%lf, z=%lf, node[%d] x=%lf, y=%lf, z=%lf\n",
//		x, y, z, NodeNumber, Node[NodeNumber].X[0],  Node[NodeNumber].X[1],  Node[NodeNumber].X[2]);
}


void fill_grid()
{
	int i=0, nTot=0;

	nTot = Settings.n_haloes;

	for(i=0; i<nTot; i++)
		assign_halo_to_node(i);

}
