#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <malloc.h>
#include <mpi.h>
#include "general.h"

#include "../general_variables.h"


int ThisTask;
int NTask;

int *SizeDispl;
int *SizeHaloes;
int *SizeDisplStructHalo;
int *SizeHaloesStructHalo;



void copy_halo_url(char *url)
{
	pUrls[ThisTask].halo_file = (char*) calloc(strlen(url)-1, sizeof(char));
	strcpy(pUrls[ThisTask].halo_file, url);
}



void init_comm_structures()
{
	pSettings = (struct general_settings *) calloc(NTask, sizeof(struct general_settings));
	pHaloes = (struct halo **) calloc(NTask, sizeof(struct halo *));
	pUrls = (struct internal_urls *) calloc(NTask, sizeof(struct internal_urls));
	pFC = (struct full_catalogue *) calloc(NTask, sizeof(struct full_catalogue));
	SizeDispl = (int*) calloc(NTask, sizeof(int));			
	SizeHaloes = (int*) calloc(NTask, sizeof(int));
	SizeDisplStructHalo = (int*) calloc(NTask, sizeof(int));			
	SizeHaloesStructHalo = (int*) calloc(NTask, sizeof(int));
}



void free_comm_structures()
{
	fprintf(stdout, "\nTask=%d, freeing memory allocated for communication...\n", ThisTask);

	free(pSettings);
	free(pHaloes);
	free(pUrls);
	free(pFC);
	free(SizeDispl);
	free(SizeDisplStructHalo);
	free(SizeHaloes);
	free(SizeHaloesStructHalo);
}



void gather_halo_structures()
{	
	int n, TaskHaloes, TaskOverThreshold, TaskVirialized, TaskSpin, TaskConcentration, TaskAll;
	
	TaskHaloes = pSettings[ThisTask].n_haloes;
	TaskOverThreshold = pSettings[ThisTask].n_threshold;
	TaskVirialized = pSettings[ThisTask].n_virialized;
	TaskSpin = pSettings[ThisTask].n_spin;
	TaskConcentration = pSettings[ThisTask].n_concentration;
	TaskAll = pSettings[ThisTask].n_all;

	MPI_Reduce(&TaskHaloes, &Settings.n_haloes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	MPI_Reduce(&TaskOverThreshold, &Settings.n_threshold, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	MPI_Reduce(&TaskVirialized, &Settings.n_virialized, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	MPI_Reduce(&TaskSpin, &Settings.n_spin, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	MPI_Reduce(&TaskConcentration, &Settings.n_concentration, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	MPI_Reduce(&TaskAll, &Settings.n_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);	

	MPI_Bcast(&Settings, sizeof(struct general_settings), MPI_BYTE, 0, MPI_COMM_WORLD);
		
	pSettings[ThisTask].n_haloes_size = pSettings[ThisTask].n_haloes * sizeof(struct halo);

		if(ThisTask==0)
		{	
			fprintf(stdout, "\nGathering %d haloes...\n", Settings.n_haloes);
			fprintf(stdout, "Gathering %d virialized haloes...\n", Settings.n_virialized);
			fprintf(stdout, "Gathering %d haloes over threshold...\n", Settings.n_threshold);
			fprintf(stdout, "Gathering %d haloes satisfying spin criterion...\n", Settings.n_spin);
			fprintf(stdout, "Gathering %d haloes satisfying concentration criterion...\n", 
				Settings.n_concentration);
			fprintf(stdout, "Gathering %d haloes satisfying all...\n", Settings.n_all);

			SizeDispl[0] = 0;
			SizeDisplStructHalo[0] = 0;
	
			haloes = (struct halo*) calloc(Settings.n_haloes, sizeof(struct halo));
		}

		MPI_Gather(&pSettings[ThisTask].n_haloes, 1, MPI_INT, SizeHaloes, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	
		MPI_Gather(&pSettings[ThisTask].n_haloes_size, 1, MPI_INT, 
			SizeHaloesStructHalo, 1, MPI_INT, 0, MPI_COMM_WORLD);	

				// When using Gatherv the Displ and Haloes size vectors are relevant only at root
			if(ThisTask==0)
			{
				for(n=1; n<NTask; n++)
				{	
					SizeDispl[n] = SizeHaloes[n-1]+SizeDispl[n-1];
					SizeDisplStructHalo[n] = SizeHaloesStructHalo[n-1]+SizeDisplStructHalo[n-1];
				}
			}

	MPI_Gatherv(&pHaloes[ThisTask][0], pSettings[ThisTask].n_haloes*sizeof(struct halo), MPI_BYTE, 
		haloes, SizeHaloesStructHalo, SizeDisplStructHalo, MPI_BYTE, 0, MPI_COMM_WORLD);

}
