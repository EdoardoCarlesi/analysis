#include "general_def.h"

extern int *cross_correlated_index;

struct general_settings Settings, *pSettings;
struct internal_urls Urls, *pUrls;
struct cosmology Cosmo;

struct halo *Haloes, **pHaloes, *crossHaloes, *tempHaloes;
struct halo_properties *HaloProperties, *SubHaloProperties;
struct sub_structure SubStructure;

struct merger_tree *MergerTree;

struct mass_function *MassFunc, *VelFunc, *ThMassFunc, *GasFunc, *NoGasFunc, *TempFunc, *DarkFunc;
struct num_density NumDen;

struct correlation_function Xi, *Xis;
struct power_spectrum *Pks;
struct growth_factor GrowthFac;
