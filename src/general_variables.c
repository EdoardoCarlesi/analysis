#include "general_variables.h"

struct general_settings Settings;
struct internal_urls Urls_internal;
struct full_catalogue FC;
struct cosmology Cosmo;

struct halo *haloes, *subhaloes;
struct halo_properties SubHaloZ, HaloZ, *HaloProperties, *SubHaloProperties;

struct merger_tree *MTree;

struct mass_function MF, AMF, *mf;
struct num_density ND;

struct correlation_function Xi, *Xis;
struct power_spectrum *Pks;
struct growth_factor GF;
