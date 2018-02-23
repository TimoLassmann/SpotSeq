#ifndef FAST_HMM_PARAM_H
#define FAST_HMM_PARAM_H


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "global.h"
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

struct fast_hmm_param{
        struct fast_t_item** list;
        float** emission;
        uint16_t last_state; 
        uint32_t alloc_items;
        uint32_t num_items;
        
        uint32_t alloc_num_states;
        uint32_t L;
};

/* Housekeeping function */
extern struct fast_hmm_param* alloc_fast_hmm_param(int k,int L);
extern int expand_fast_hmm_param_if_necessary(struct fast_hmm_param* ft, int k);
extern void free_fast_hmm_param(struct fast_hmm_param* ft);


/* Sorting and binary search  */

extern int fast_hmm_param_cmp_by_t_desc(const void *a, const void *b);
extern int fast_hmm_param_cmp_by_to_from_asc(const void *a, const void *b);
extern int fast_hmm_param_cmp_by_from_asc(const void *a, const void *b);
extern int fast_hmm_param_cmp_by_to_asc(const void *a, const void *b);

/* return index of first element < x i.e. we can then do for(i =0; i < return;i++) */

extern int fast_hmm_param_binarySearch_t(struct fast_hmm_param* ft, float x);

/* These functions return the first and last+1 entry in list that has value of x */
extern int fast_hmm_param_binarySearch_to_lower_bound(struct fast_hmm_param* ft, int x);
extern int fast_hmm_param_binarySearch_to_upper_bound(struct fast_hmm_param* ft, int x);      
extern int fast_hmm_param_binarySearch_from_lower_bound(struct fast_hmm_param* ft, int x);
extern int fast_hmm_param_binarySearch_from_upper_bound(struct fast_hmm_param* ft, int x);

#endif
