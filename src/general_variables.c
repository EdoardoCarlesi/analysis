#include "general_variables.h"

struct general_settings Settings, *pSettings;
struct internal_urls Urls, *pUrls;
struct cosmology Cosmo;

struct nfw NFW;

struct halo *Haloes, *SubHaloes, **pHaloes;
struct halo_properties SubHaloZ, HaloZ, *HaloProperties, *SubHaloProperties;

struct merger_tree *MergerTree;

struct mass_function MassFunc, ThMassFunc, *MassFuncZ;
struct num_density NumDen;

struct correlation_function Xi, *Xis;
struct power_spectrum *Pks;
struct growth_factor GrowthFac;
