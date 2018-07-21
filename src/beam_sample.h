#ifndef BEAM_SAMPLE_H
#define BEAM_SAMPLE_H


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "distributions.h"
#include <math.h>
#include <float.h>
#include <stdint.h>

#include "tldevel.h"
#include "thr_pool.h"
#include "rbtree.h"
#include "fast_hmm_param.h"
#include "ihmm_seq.h"
#include "model.h"
#include "global.h"



/* main runner function  */
int run_beam_sampling(struct ihmm_model* model, struct seq_buffer* sb, struct fast_hmm_param* ft, struct thr_pool* pool, int iterations, int num_threads); 
/* key Operations  */
extern int fill_fast_transitions(struct ihmm_model* model,struct fast_hmm_param* ft);
extern int fill_fast_transitions_only_matrices(struct ihmm_model* model,struct fast_hmm_param* ft);
extern int add_state_from_fast_hmm_param(struct ihmm_model* ihmm,struct fast_hmm_param* ft);

extern int fill_background_emission(struct fast_hmm_param*ft,struct seq_buffer* sb);
extern int fill_background_emission_from_model(struct fast_hmm_param*ft, struct ihmm_model* model);
#endif
