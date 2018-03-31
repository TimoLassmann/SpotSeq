#ifndef FAST_HMM_PARAM_H
#define FAST_HMM_PARAM_H


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"
#include "distributions.h"
#include <math.h>
#include <float.h>
#include <stdint.h>


struct fast_t_item{
        float t;
        uint16_t from;
        uint16_t to;
};

struct 
fast_hmm_param{
        struct fast_t_item** list;
        float** emission;
        int8_t* active_states;
        uint16_t last_state; 
        uint32_t alloc_items;
        uint32_t num_items;
        
        uint32_t alloc_num_states;
};

/* Housekeeping function */
extern struct fast_hmm_param* alloc_fast_hmm_param(void);
extern int expand_fast_hmm_param_if_necessary(struct fast_hmm_param* ft, int k);
extern void free_fast_hmm_param(struct fast_hmm_param* ft);

/* key Operations  */
extern int add_state_from_fast_hmm_param(rk_state rndstate,struct fast_hmm_param* ft, float* beta, float alpha, float gamma);
/* Sorting and binary search  */
#endif
