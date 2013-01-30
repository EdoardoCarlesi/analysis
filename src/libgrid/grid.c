#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../general_variables.h"
#include "../general_functions.h"
#include "grid.h"

#ifdef WITH_MPI
#include <mpi.h>
#include "../libparallel/general.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define PTS_CIC 8

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
	Grid.cell_volume = Grid.grid_size * Grid.grid_size * Grid.grid_size;

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
					Node[n].M_CIC = 0;
					Node[n].X_cm[0] = 0;
					Node[n].X_cm[1] = 0;
					Node[n].X_cm[2] = 0;
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



void NGP_assignment(int n)
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

	Node[NodeNumber].X_cm[0] += (Haloes[n].Xc * Haloes[n].Mvir) / Node[NodeNumber].M ; 
	Node[NodeNumber].X_cm[1] += (Haloes[n].Yc * Haloes[n].Mvir) / Node[NodeNumber].M ; 
	Node[NodeNumber].X_cm[2] += (Haloes[n].Zc * Haloes[n].Mvir) / Node[NodeNumber].M ; 

//	fprintf(stdout, "Halo x=%lf, y=%lf, z=%lf, node[%d] x=%lf, y=%lf, z=%lf\n",
//		x, y, z, NodeNumber, Node[NodeNumber].X[0],  Node[NodeNumber].X[1],  Node[NodeNumber].X[2]);
}



void find_nodes_CIC(int *node_numbers, int n)
{
	int N=0, m=0, h=0, j=0, k=0, nearest=0;
	int i[3], nx[3][2];
	double x[3];
	
	N = Grid.N;
	
	INFO_MSG("Finding CIC nodes");

		x[0] = Haloes[n].Xc;
		x[1] = Haloes[n].Yc;
		x[2] = Haloes[n].Zc;
	
		for(m=0; m<3; m++)
			i[m] = (int) lrint(abs(x[m] - Grid.half_grid)/ Grid.grid_size);
	
			nearest = N*N*i[0] + N*i[1] + i[2];

	//		fprintf(stdout, "x=%lf, y=%lf, z=%lf, Nearest grid=%d\n", x[0], x[1], x[2], nearest);		
		
			for(m=0; m<3; m++)
			{
				nx[m][0] = 0;
	
				if(x[m] - Node[nearest].X[m] > 0)
					nx[m][1] = 1;
						else
							nx[m][1] = -1;
			}	
		
			m=0;	
		
			for(h=0; h<2; h++)
			{
				for(j=0; j<2; j++)
				{
					for(k=0; k<2; k++)
					{
						if	(
								i[0] + nx[0][h] < 0 || i[0] + nx[0][h] > N-1 ||
								i[1] + nx[1][j] < 0 || i[1] + nx[1][j] > N-1 ||
								i[2] + nx[2][k] < 0 || i[2] + nx[2][k] > N-1 
							)
						{

							node_numbers[m] = -1;

						}	 
							else 
						{
							node_numbers[m] = N*N*(i[0] + nx[0][h]) +
								N*(i[1] + nx[1][j]) + i[2] + nx[2][k];
						}

				//	fprintf(stdout, "%d) node_number[%d]=%d\n", n, m, node_numbers[m]);
						m++;
					}	
				}
	}
	
	INFO_MSG("Done");
}



void set_fac_CIC(double *r, double *m_fac)
{
	int i=0, cells=0;
	double Len=0., Ltot=0.;

	Len = 10*Grid.grid_size;
	
		for(i=0; i<PTS_CIC; i++)
		{
			if(abs(r[i]) < Len && r[i] != -1)
			{
				Ltot += r[i];	
				cells++;
			}
		}
		//	fprintf(stderr, "Ltot=%lf, Len=%lf, cells=%d\n", Ltot, Len, cells);			

			for(i=0; i<PTS_CIC; i++)
			{
				if(abs(r[i]) < Len && r[i] != -1)
				{
					m_fac[i] = 1 - r[i]/Ltot;

						if(cells > 1)
					m_fac[i] *= 1./( (double) cells - 1.);
				}
			else
		{
			m_fac[i] = 0;
		}
	}
}



void CIC_assignment(int n)
{
	int i=0, index=0;
	double x, y, z;

	int NodeNumber[PTS_CIC];
	double dist[PTS_CIC], mass_fac[PTS_CIC];

	INFO_MSG("CIC assignment");
	
		find_nodes_CIC(NodeNumber, n);

		x = Haloes[n].Xc;
		y = Haloes[n].Yc;
		z = Haloes[n].Zc;

//	fprintf(stdout, "Halo x=%lf, y=%lf, z=%lf, node[%d] x=%lf, y=%lf, z=%lf\n",
//		x, y, z, NodeNumber, Node[NodeNumber].X[0],  Node[NodeNumber].X[1],  Node[NodeNumber].X[2]);
		
		for(i=0; i<PTS_CIC; i++)
		{	
			index = NodeNumber[i];

				if(index != -1)
				{
					dist[i] = sqrt(	(x-Node[index].X[0])*(x-Node[index].X[0]) +
							(y-Node[index].X[1])*(y-Node[index].X[1]) +
							(z-Node[index].X[2])*(z-Node[index].X[2]));		
				}	
						else
					{
						dist[i] = -1;
					}

	//			fprintf(stdout, "%d) dist[%d]=%lf\n", n, i, dist[i]);

				}
	
			set_fac_CIC(dist, mass_fac);
	
		for(i=0; i<PTS_CIC; i++)
		{
			index = NodeNumber[i];
			
		if(index != -1)
		{
			Node[index].M_CIC += mass_fac[i] * Haloes[n].Mvir;
	//		fprintf(stderr, "%d) Node[%d] has m_fac=%lf\n", n, index, mass_fac[i]);
		}
	}	

}	



void fill_grid_NGP()
{
	int i=0, nTot=0;

		nTot = Settings.n_haloes;


#ifdef _OPENMP
		omp_set_num_threads(OMP_THREADS);
#endif

#		pragma omp parallel for 		\
	 	private(i)

			for(i=0; i<nTot; i++)
				NGP_assignment(i);

}



void fill_grid_CIC()
{
	int i=0, nTot=0;

		nTot = Settings.n_haloes;
		//	nTot = 10;

#ifdef _OPENMP
		omp_set_num_threads(OMP_THREADS);
#endif

#		pragma omp parallel for 		\
		private(i)
	
			for(i=0; i<nTot; i++)
				CIC_assignment(i);
}



void free_nodes()
{
	free(Node);
}
