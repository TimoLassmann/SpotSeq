#ifndef THREAD_DATA_H
#define THREAD_DATA_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"
#include "randomkit.h"


#define THREAD_DATA_FULL 1
#define THREAD_DATA_BEAM 2


struct seqer_thread_data{
        struct fast_param_bag* ft_bag;
        struct fast_hmm_param* ft;
        struct seq_buffer* sb;
        struct fhmm* fhmm;
        double** dyn;
        double** F_matrix;
        double** B_matrix;
        double** t;             /* for forward backward - to collect estimated */
        double** e;             /* for forward backward - to collect estimated */
        //int num_seq;
        int thread_ID;
        int num_threads;
        unsigned int seed;
        rk_state rndstate;
};

extern struct seqer_thread_data** create_seqer_thread_data(int* num_threads, int max_len, int K,rk_state* random, int mode);
extern int resize_seqer_thread_data(struct seqer_thread_data** td,int* num_threads, int max_len, int K);
extern int compare_seqer_thread_data(struct seqer_thread_data** a , struct seqer_thread_data** b, int num);

extern void free_seqer_thread_data(struct seqer_thread_data** td);

#endif
