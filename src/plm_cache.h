#include "floattypes.h"
#include "integration_drude.h"
#include "libcasimir.h"
#include "sfunc.h"


void plm_create_cache(integration_drude_t* int_drude);

void plm_destroy_cache(integration_drude_t* int_drude);


void plm_cache_PlmPlm(struct integ_context* context, float80 x, plm_combination_t* comb,
                      unsigned int index, unsigned int iteration);


void plm_cache_init(struct integ_context* context);

void plm_cache_free(integration_drude_t* int_drude);
