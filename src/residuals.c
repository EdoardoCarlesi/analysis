#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libio/io.h"
#include "libmath/math.h"
#include "libcosmo/cosmo.h"
#include "libhalo/halo.h"

#define CHAR_SIZE 1024
#define fac_max 0.99
#define fac_min 1.0

int SKIP;
int EDGE;
int USE_COLS;


struct vector
{
	int N;
	double *x;
	double *y;

} Vec01, Vec02, VecOUT;


void calloc_vec(struct vector *V, int D)
{
	V->N = D;
	V->x = calloc(D, sizeof(double));
	V->y = calloc(D, sizeof(double));
}


void print_vec(struct vector *V)
{
	int i=0;
	
		for(i=0; i<V->N; i++)
			fprintf(stderr, "%d) %e %e\n", i, V->x[i], V->y[i]);
	
}


void print_vec_to_file(struct vector *V, char *fname)
{
	int i=0;
	FILE *fout = NULL;
	fout = fopen(fname, "w");
	
		if(fout == NULL)
			fprintf(stderr, "Error: could not open file %s to write output\n", fname);	

		for(i=0; i<V->N; i++)
			fprintf(fout, "%e\t%e\n", V->x[i], V->y[i]);
	
	fclose(fout);
}


void load_file_into_vec(FILE *f, struct vector *V, int colX, int colY, int totCol)
{
	char dummy[CHAR_SIZE];
	double *IN;
	int i=0, j=0, count = -1;
	
	//	fprintf(stderr,"x=%d y=%d tot=%d skip=%d\n", colX, colY, totCol, SKIP);
		IN = calloc(totCol, sizeof(double));	

		for(j=0; j<V->N+SKIP; j++)
		{
			
			if(j<SKIP)
				fgets(dummy, CHAR_SIZE, f);
			//	fprintf(stderr, "%d) %s\n", j, dummy);

			else
			{
				for(i=0; i<totCol; i++)
				{
					fscanf(f, "%lf\t", &IN[i]);
					//sscanf(dummy, "%lf\t", &IN[i]);
					//fprintf(stderr, "(%d) %e\t", j, IN[i]);
				}

					V->x[j-SKIP] = IN[colX];
					V->y[j-SKIP] = IN[colY];

				if(USE_COLS == 0)
					fgets(dummy, CHAR_SIZE, f);
			}		
			//fprintf(stderr,"x=%d y=%d skip=%d\n", colX, colY, SKIP);
			//fprintf(stderr,"\n");
		}

		//print_vec(V);

	rewind(f);
}


double get_vec_max(struct vector *V1, struct vector *V2)
{
	double m1 = maximum(V1->x, V1->N);
	double m2 = maximum(V2->x, V2->N);

	// Get the minimum of the two maxima
	if(m2 < m1)
		return m2;
	else
		return m1;
}


double get_vec_min(struct vector *V1, struct vector *V2)
{
	double m1 = minimum(V1->x, V1->N);
	double m2 = minimum(V2->x, V2->N);

	// Get the maximum of the two minima
	if(m2 > m1)
		return m2;
	else
		return m1;
}


/*
 *  USE:
 *  ./residuals file01 file02 fileOUTPUT columnX_1 columX_2 columnY_1 columnY_2 columns_tot n_points_interp bins skip 
  * */
int main(int argc, char **argv)
{
	int i=0, count = 1;

	char *fname01 = argv[count++];
	char *fname02 = argv[count++];
	char *fnameOUT = argv[count++];

	int col_x_01 = atoi(argv[count++])-1;
	int col_x_02 = atoi(argv[count++])-1;
	int col_y_01 = atoi(argv[count++])-1;
	int col_y_02 = atoi(argv[count++])-1;
	int Ntot_col = atoi(argv[count++]);

	int D_out = atoi(argv[count++]);
	int binning = atoi(argv[count++]); // 0 = lin bin, 1 = log bin

	SKIP = atoi(argv[count++]);
	EDGE = atoi(argv[count++]);
	USE_COLS = atoi(argv[count++]);

	int D_01, D_02;
	double min, max, x, y1, y2;	

	//for(i=1; i<count; i++)
	//	fprintf(stderr, "%d) argv=%s\n", i, argv[i]);

	FILE *f01 = fopen(fname01, "r");
	FILE *f02 = fopen(fname02, "r");

	// This also checks that the file has been read in correctly
	D_01 = get_lines(f01, fname01) - SKIP;
	D_02 = get_lines(f02, fname02) - SKIP;

	//fprintf(stderr, "File01 N=%d, x=%d y=%d\n", D_01, col_x_01, col_y_01);
	//fprintf(stderr, "File02 N=%d\n", D_02);

		calloc_vec(&Vec01, D_01);
		calloc_vec(&Vec02, D_02);
		calloc_vec(&VecOUT, D_out);
		
		load_file_into_vec(f01, &Vec01, col_x_01, col_y_01, Ntot_col);
		load_file_into_vec(f02, &Vec02, col_x_02, col_y_02, Ntot_col);

		max = fac_max * get_vec_max (&Vec01, &Vec02);
		min = fac_min * get_vec_min (&Vec01, &Vec02);

			if(binning == 0)
				VecOUT.x = lin_stepper(min, max, VecOUT.N);

			if(binning == 1)
				VecOUT.x = log_stepper(min, max, VecOUT.N);

		//	print_vec(&Vec01);
		//	print_vec(&Vec02);

			for(i=0; i<VecOUT.N; i++)
			{
				x = VecOUT.x[i];
				y1 = get_interpolated_value(Vec01.x, Vec01.y, Vec01.N, x);
				y2 = get_interpolated_value(Vec02.x, Vec02.y, Vec02.N, x);
				//fprintf(stderr, "%d) x=%e, y1=%e, y2=%e\n", i, x, y1, y2);
				VecOUT.y[i] = y1 / y2;

				if(EDGE == 1 && i == VecOUT.N-1)
					VecOUT.y[i] = Vec01.y[Vec01.N-1] / Vec02.y[Vec02.N-1];
		
			}

		//print_vec(&VecOUT);
		print_vec_to_file(&VecOUT, fnameOUT);

	return 0;
}


