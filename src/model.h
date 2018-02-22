#ifndef MODEL_H
#define MODEL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "global.h"
#include "tldevel.h"
#include "distributions.h"
#include <math.h>
#include <float.h>
#include <stdint.h>


#define IHMM_PARAM_PLACEHOLDER -9999.99f

struct seq_buffer;          /* forward declaration  */

struct ihmm_model{
        float** transition_counts;
        float** emission_counts;
        float* beta;
        rk_state rndstate;
        float gamma;
        float alpha;
        float alpha0_a;
        float alpha0_b;
        float gamma_a;
        float gamma_b;
        int num_states;         /* this excludes the start and stop states (0,1) */
        int alloc_num_states;
        int L;
};

/* Housekeeping */

extern struct ihmm_model* alloc_ihmm_model(int K, int L);
extern int resize_ihmm_model(struct ihmm_model* ihmm, int K);
extern void free_ihmm_model(struct ihmm_model* ihmm);




/* Fill counts from sequences  */
extern int fill_counts(struct ihmm_model* ihmm, struct seq_buffer* sb);

extern int remove_unused_states_labels(struct ihmm_model* ihmm, struct seq_buffer* sb);
/* re-estimate hyper parameters */
extern int iHmmHyperSample(struct ihmm_model* model, int iterations);

/* help functions */
extern int print_counts(struct ihmm_model* ihmm);
extern int print_model_parameters(struct ihmm_model* ihmm);


#endif
