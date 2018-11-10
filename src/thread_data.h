#ifndef THREAD_DATA_H
#define THREAD_DATA_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"
#include "randomkit.h"

struct spotseq_thread_data{
        struct fast_param_bag* ft_bag;
        struct fast_hmm_param* ft;
        struct seq_buffer* sb;
        struct fhmm* fhmm;
        float** dyn;
        float** F_matrix;
        float** B_matrix;
        float** t;
        float** e;
        int thread_ID;
        int num_threads;
        unsigned int seed;
        rk_state rndstate;
};

extern  struct spotseq_thread_data** create_spotseq_thread_data(int* num_threads, int max_len, int K,rk_state* random);
extern int resize_spotseq_thread_data(struct spotseq_thread_data** td,int* num_threads, int max_len, int K);
extern void free_spotseq_thread_data(struct spotseq_thread_data** td, int num_threads);

#endif
