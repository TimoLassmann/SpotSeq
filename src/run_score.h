#ifndef RUN_SCORE_H
#define RUN_SCORE_H



#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>


#include "tldevel.h"
#include "thr_pool.h"
#include "ihmm_seq.h"

#include "finite_hmm.h"

#include "thread_data.h"


extern int run_score_sequences(struct fhmm* fhmm, struct seq_buffer* sb, int num_threads);
extern void* do_score_sequences(void* threadarg);

#endif
