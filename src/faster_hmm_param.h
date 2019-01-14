#ifndef FASTER_HMM_PARAM_H
#define FASTER_HMM_PARAM_H

/* The idea here is to use libavl's red-black tree implementations to:

1) insert sort transition parameters as they are sampled
2) manage addition of transitions to new nodes as they are instantiated
3) use the traverse functions to have a one time pass over all sampled u values & set boundaries on transitions
- more explanation is needed:
during dynamic programming, at each position in the sequence(s) transitions with a probablity > u are sampled. Iterating over all transitions and checking if t[i][j] > u is horribly inefficient. Instead we sort transitions based on their probabilities once and simply determine which top "x" transitions are used. In the fast_hmm_param implementation I did this by performing a binary search for the first transition with t[i][j] < u and applied all preceding transition. This is fast but involves N(number of residues in all sequences) binary searches in the inner loop.

Here I will follow the same strategy but will perform a single linear pass over all transitions & u values to pre-compute the boundaries (to be stored in u (maybe... ))

*/

struct faster_t_item{
        double t;
        uint16_t a;
        uint16_t b;;
};


struct faster_hmm_param{
        struct faster_t_item** infinity;
        struct rbtree_root* root;
        double** transition;
        double** emission;
        double* background_emission;
        int last_state;
        int alloc_items;
        int num_items;
        int alloc_num_states;
        int L;
};


struct faster_param_bag{
        struct faster_hmm_param** fast_params;
        int max_last_state;
        int num_models;
};

extern struct faster_param_bag* alloc_fast_param_bag(int num_models, int* K, int L);
extern void free_faster_param_bag(struct faster_param_bag* b);


extern struct faster_hmm_param* alloc_faster_hmm_param(int k,int L);
extern void free_faster_hmm_param(struct faster_hmm_param* ft);

#endif
