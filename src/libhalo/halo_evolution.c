#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "halo_evolution.h"
#include "halo_properties.h"
#include "subhalo_general.h"
#include "subhalo_properties.h"
#include "../libio/halo_io.h"
#include "../general_variables.h"


void compute_halo_and_subhalo_statistics(int j)
{
	fprintf(stdout, "Reading halo url[%d]: %s\n", j, Urls.urls[j]);

	Settings.use_cat = j;
	set_halo_url();

		read_halo_file();

		compute_halo_properties();
		compute_subhalo_properties();

		copy_halo_and_subhalo_properties(j);

		stdout_halo_status(j);

	free_halo_properties();
	free_subhalo_properties();
}



void stdout_halo_status(int j)
{
	fprintf(stdout, "%lf",   HaloProperties[j].z);   
	fprintf(stdout, "\t%lf", HaloProperties[j].c_0);   
	fprintf(stdout, "\t%lf", HaloProperties[j].l_0);   
	fprintf(stdout, "\t%lf", HaloProperties[j].s0);   
	fprintf(stdout, "\t%lf", HaloProperties[j].t0);   
	fprintf(stdout, "\t%lf", HaloProperties[j].avgSub);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].c_0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].l_0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].s0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].t0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].costh0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].cosphi0);   
	fprintf(stdout, "\t%lf", SubHaloProperties[j].vel_0);   
	fprintf(stdout, "\t%e",  SubHaloProperties[j].avgMass);   
	fprintf(stdout, "\n");
}



void initialize_halo_storage()
{
	int k=0, j=0, nTot=0;

	nTot = Urls.nCatalogueFiles;

	HaloProperties = (struct halo_properties *) calloc(nTot, sizeof(struct halo_properties));
	SubHaloProperties = (struct halo_properties *) calloc(nTot, sizeof(struct halo_properties));

		for(j=0; j<nTot; j++)
		{
			k = GrowthFac.npts - j - 1;
			HaloProperties[j].z = GrowthFac.z[k];
			SubHaloProperties[j].z = GrowthFac.z[k];
		}
}



void copy_halo_and_subhalo_properties(int j)
{
	fprintf(stdout, "\nStoring halo average properties at z=%lf.\n", GrowthFac.z[j]);
	HaloProperties[j].avgSub = HaloZ.avgSub;
	HaloProperties[j].l_0 = HaloZ.l_0;
	HaloProperties[j].c_0 = HaloZ.c_0;
	HaloProperties[j].t0 = HaloZ.t0;
	HaloProperties[j].s0 = HaloZ.s0;

	fprintf(stdout, "\nStoring subhalo average properties at z=%lf.\n", GrowthFac.z[j]);
	SubHaloProperties[j].l_0 = SubHaloZ.l_0;
	SubHaloProperties[j].t0 = SubHaloZ.t0;
	SubHaloProperties[j].s0 = SubHaloZ.s0;
	SubHaloProperties[j].cosphi0 = SubHaloZ.cosphi0;
	SubHaloProperties[j].costh0 = SubHaloZ.costh0;
	SubHaloProperties[j].vel_0 = SubHaloZ.vel_0;
	SubHaloProperties[j].avgMass = SubHaloZ.avgMass;
}
