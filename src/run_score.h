#ifndef RUN_SCORE_H
#define RUN_SCORE_H


#include <stdlib.h>
#include <stdio.h>

#include <omp.h>

#include "tldevel.h"
//#include "thr_pool.h"

#include "finite_hmm.h"

#include "thread_data.h"

#include "model_struct.h"

extern int run_score_sequences(struct fhmm** fhmm, struct tl_seq_buffer* sb,struct seqer_thread_data** td, int n_model,int mode);

//extern int score_all_vs_all(struct model_bag* mb, struct tl_seq_buffer* sb, struct seqer_thread_data** td);
//extern int run_score_sequences(struct fhmm* fhmm, struct tl_seq_buffer* sb,struct seqer_thread_data** td);
extern void* do_score_sequences(void* threadarg);
//extern void* do_score_sequences_per_model(void* threadarg);

extern int run_label_sequences(struct fhmm* fhmm, struct tl_seq_buffer* sb, int num_threads);
//extern int run_label_sequences(struct fhmm* fhmm, struct seq_buffer* sb, int num_threads);
extern void* do_label_sequences(void* threadarg);

#endif
