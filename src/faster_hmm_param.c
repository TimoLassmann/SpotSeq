#include "faster_hmm_param.h"


struct faster_param_bag* alloc_faster_param_bag(int num_models, int* K, int L)
{
        struct fast_paramer_bag* b = NULL;
        int i;
        ASSERT(num_models > 0, "No models");

        MMALLOC(b, sizeof(struct fast_param_bag));

        b->fast_params = NULL;
        b->num_models = num_models;
        b->max_last_state = -1;

        MMALLOC(b->fast_params,sizeof(struct faster_hmm_param*)* b->num_models);

        for(i = 0; i < b->num_models;i++){
                b->fast_params[i] = NULL;
                RUNP(b->fast_params[i] = alloc_faster_hmm_param(K[i],L));
        }

        return b;
ERROR:
        return NULL;
}
