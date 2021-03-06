#ifndef FAST_HMM_PARAM_H
#define FAST_HMM_PARAM_H



#include "global.h"
#include "tldevel.h"
#include "tlrbtree.h"
#include "distributions.h"
#include <math.h>
#include <float.h>
#include <stdint.h>


struct fast_t_item{
        double t;
        uint16_t from;
        uint16_t to;
};

struct fast_hmm_param{
        struct fast_t_item** list;
        struct fast_t_item** infinity;
        // struct rbtree_root* root;
        double** transition;
        double** emission;
        //double* background_emission;
        int num_trans;
        int alloc_num_trans;
        int last_state;
        int alloc_items;
        int num_items;
        int alloc_num_states;
        int L;
};

/* For when we want to train multiple models  */
struct fast_param_bag{
        struct fast_hmm_param** fast_params;
        int max_last_state;
        int num_models;
};


/* Housekeeping function */
extern struct fast_param_bag* alloc_fast_param_bag(int num_models,  int L);

//extern struct fast_param_bag* alloc_fast_param_bag(int num_models, int* K, int L);
extern void free_fast_param_bag(struct fast_param_bag* b);

extern int expand_num_trans(struct fast_hmm_param* ft);
extern struct fast_hmm_param* alloc_fast_hmm_param(int k,int L);



//extern int expand_fast_hmm_param_if_necessary(struct fast_hmm_param* ft, int new_num_states,int new_items);

extern int expand_ft_if_necessary(struct fast_hmm_param* ft, int new_num_states);
//extern int expand_emission_if_necessary(struct fast_hmm_param* ft, int new_num_states);
//extern int expand_transition_if_necessary(struct fast_hmm_param* ft);

extern void free_fast_hmm_param(struct fast_hmm_param* ft);


/* turn RB tree into a flat indexable structure...  */
extern int make_flat_param_list(struct fast_hmm_param* ft);

/* Sorting and binary search  */

extern int fast_hmm_param_cmp_by_t_desc(const void *a, const void *b);
extern int fast_hmm_param_cmp_by_to_from_asc(const void *a, const void *b);
extern int fast_hmm_param_cmp_by_from_asc(const void *a, const void *b);
extern int fast_hmm_param_cmp_by_to_asc(const void *a, const void *b);

/* return index of first element < x i.e. we can then do for(i =0; i < return;i++) */

extern int fast_hmm_param_binarySearch_t(struct fast_hmm_param* ft, double x);

/* These functions return the first and last+1 entry in list that has value of x */
extern int fast_hmm_param_binarySearch_to_lower_bound(struct fast_hmm_param* ft, int x);
extern int fast_hmm_param_binarySearch_to_upper_bound(struct fast_hmm_param* ft, int x);
extern int fast_hmm_param_binarySearch_from_lower_bound(struct fast_hmm_param* ft, int x);
extern int fast_hmm_param_binarySearch_from_upper_bound(struct fast_hmm_param* ft, int x);

#endif
