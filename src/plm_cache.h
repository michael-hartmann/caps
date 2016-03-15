#include "floattypes.h"
#include "integration_drude.h"
#include "libcasimir.h"
#include "sfunc.h"


void plm_create_cache(casimir_t* casimir);

void plm_destroy_cache(void);


void plm_cache_PlmPlm(struct integ_context* context, float80 x, plm_combination_t* comb,
                      unsigned int index, unsigned int iteration);


void plm_cache_init(struct integ_context* context, int n);

void plm_cache_free(struct integ_context* context);
