#ifndef BEAM_SAMPLE_H
#define BEAM_SAMPLE_H


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "distributions.h"
#include <math.h>
#include <float.h>
#include <stdint.h>

#include <omp.h>

#include "tldevel.h"
//#include "thr_pool.h"
//#include "rbtree.h"
//#include "fast_hmm_param.h"
#include "ihmm_seq.h"
#include "model.h"
#include "global.h"

#include "hmm_conversion.h"
#include "finite_hmm.h"

#include "thread_data.h"

/* main runner function  */
int run_beam_sampling(struct model_bag* model_bag, struct fast_param_bag*ft_bag, struct seq_buffer* sb,struct wims_thread_data** td, int iterations, int num_threads);


#endif
