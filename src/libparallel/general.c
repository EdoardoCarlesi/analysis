#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <malloc.h>
#include <mpi.h>

#include "../general_def.h"

#include "general.h"


int ThisTask;
int NTask;

int *SizeDisplStructHalo;
int *SizeHaloesStructHalo;

struct Cpu *cpu;


void generate_url_for_tasks()
{
	char halo_list_task[200];
	char profile_list_task[200];
	char subhalo_list_task[200];
	char command[250];	

		sprintf(halo_list_task, "%s.%s",Urls.halo_list, cpu[ThisTask].name);
		sprintf(profile_list_task, "%s.%s",Urls.profile_list, cpu[ThisTask].name);
		sprintf(subhalo_list_task, "%s.%s",Urls.subhalo_list, cpu[ThisTask].name);

		TASK_INFO_MSG(ThisTask, "reading halo, profiles and substructure files lists...");
	
		sprintf(command, "%s s/%s/%s/ <%s >%s", 
			"sed ", "0000", cpu[ThisTask].name, Urls.halo_list, halo_list_task);
				system(command);

		sprintf(command, "%s s/%s/%s/ <%s >%s", 
			"sed ", "0000", cpu[ThisTask].name, Urls.profile_list, profile_list_task);
				system(command);

		sprintf(command, "%s s/%s/%s/ <%s >%s", 
			"sed ", "0000", cpu[ThisTask].name, Urls.subhalo_list, subhalo_list_task);
				system(command);

	pUrls[ThisTask].halo_list = (char*) calloc(strlen(halo_list_task)+1, sizeof(char));
		strcpy(pUrls[ThisTask].halo_list, halo_list_task);

	pUrls[ThisTask].profile_list = (char*) calloc(strlen(profile_list_task)+1, sizeof(char));
		strcpy(pUrls[ThisTask].profile_list, profile_list_task);

	pUrls[ThisTask].subhalo_list = (char*) calloc(strlen(subhalo_list_task)+1, sizeof(char));
		strcpy(pUrls[ThisTask].subhalo_list, subhalo_list_task);

		TASK_INFO_MSG(ThisTask, "has generated halo, profiles and substructure files lists.");
}



void init_comm_structures()
{
	pSettings = (struct general_settings *) calloc(NTask, sizeof(struct general_settings));
	pHaloes = (struct halo **) calloc(NTask, sizeof(struct halo *));
	pUrls = (struct internal_urls *) calloc(NTask, sizeof(struct internal_urls));
	SizeDisplStructHalo = (int*) calloc(NTask, sizeof(int));			
	SizeHaloesStructHalo = (int*) calloc(NTask, sizeof(int));

	// Init the cpu read by every box to default dummy values
	pSettings[ThisTask].box.X[0][0] = 1e+5;
	pSettings[ThisTask].box.X[0][1] = 0;
	pSettings[ThisTask].box.X[1][0] = 1e+5;
	pSettings[ThisTask].box.X[1][1] = 0;
	pSettings[ThisTask].box.X[2][0] = 1e+5;
	pSettings[ThisTask].box.X[2][1] = 0;
}



void init_cpu_struct()
{
		cpu = (struct Cpu*) calloc(NTask, sizeof(struct Cpu));

		sprintf(cpu[ThisTask].name, "%04d", ThisTask);
	//	printf("CPU=%s\n",cpu[ThisTask].name); //, "%04d", ThisTask);
}



void free_comm_structures()
{
	TASK_INFO_MSG(ThisTask, "freeing memory allocated for communication");
	free(pSettings);
	free(pHaloes);
	free(pUrls);
	free(SizeDisplStructHalo);
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

			SizeDisplStructHalo[0] = 0;
	
			Haloes = (struct halo*) calloc(Settings.n_haloes, sizeof(struct halo));
		}

		MPI_Gather(&pSettings[ThisTask].n_haloes_size, 1, MPI_INT, 
			SizeHaloesStructHalo, 1, MPI_INT, 0, MPI_COMM_WORLD);	

				// When using Gatherv the Displ and Haloes size vectors are relevant only at root
			if(ThisTask==0)
			{
				for(n=1; n<NTask; n++)
				{	
					SizeDisplStructHalo[n] = SizeHaloesStructHalo[n-1]+SizeDisplStructHalo[n-1];
				}
			}

	MPI_Gatherv(&pHaloes[ThisTask][0], pSettings[ThisTask].n_haloes*sizeof(struct halo), MPI_BYTE, 
		Haloes, SizeHaloesStructHalo, SizeDisplStructHalo, MPI_BYTE, 0, MPI_COMM_WORLD);
	
if(ThisTask==0)
	INFO_MSG("Done gathering");
}
