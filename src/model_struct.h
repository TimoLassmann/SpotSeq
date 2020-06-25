#ifndef MODEL_STRUCT_H
#define MODEL_STRUCT_H

#include "distributions.h"

struct ihmm_model{
        double** transition_counts;
        double** emission_counts;
        double* beta;
        double* background;
        unsigned int seed;
        rk_state rndstate;
        double gamma;
        double alpha;
        double alpha_a;
        double alpha_b;
        double gamma_a;
        double gamma_b;
        double log_likelihood;
        double alpha_limit;
        double gamma_limit;
        int num_states;         /* this excludes the start and stop states (0,1) */
        int alloc_num_states;
        int L;
        int training_iterations;
};

struct model_bag{
        struct ihmm_model** models;
        struct fhmm** finite_models;
        double* min_u;
        int best_model;
        int max_num_states;
        int num_models;
        unsigned int seed;      /* Starting value */
        rk_state rndstate;      /* main seed */
};


#endif
