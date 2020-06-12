#ifndef THREAD_DATA_H
#define THREAD_DATA_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"
#include "randomkit.h"

struct wims_thread_data{
        struct fast_param_bag* ft_bag;
        struct fast_hmm_param* ft;
        struct seq_buffer* sb;
        struct fhmm* fhmm;
        double** dyn;
        double** F_matrix;
        double** B_matrix;
        double** t;
        double** e;
        //int num_seq;
        int thread_ID;
        int num_threads;
        unsigned int seed;
        rk_state rndstate;
};

extern struct wims_thread_data** create_wims_thread_data(int* num_threads, int max_len, int K,rk_state* random);
extern int resize_wims_thread_data(struct wims_thread_data** td,int* num_threads, int max_len, int K);
extern int compare_wims_data(struct wims_thread_data** a , struct wims_thread_data** b, int num);

extern void free_wims_thread_data(struct wims_thread_data** td);

#endif
