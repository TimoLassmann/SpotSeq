#ifndef FAST_TRANSITION_H
#define FAST_TRANSITION_H


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

struct fast_transition{
        struct fast_t_item** list;
        int8_t* active_states;
        uint16_t last_state; 
        uint32_t alloc_items;
        uint32_t num_items;
        
        uint32_t alloc_num_states;
};

/* Housekeeping function */
extern struct fast_transition* alloc_fast_transition(void);
extern int expand_fast_transition_if_necessary(struct fast_transition* ft, int k);
extern void free_fast_transition(struct fast_transition* ft);

/* key Operations  */
extern int add_state_from_fast_transition(rk_state rndstate,struct fast_transition* ft, float* beta, float alpha, float gamma);
extern int delete_state_from_fast_transition(struct fast_transition* ft, int x);
/* Sorting and binary search  */
#endif
