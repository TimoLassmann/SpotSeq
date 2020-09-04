#ifndef THREAD_DATA_H
#define THREAD_DATA_H

#include "randomkit.h"

#ifdef THREAD_DATA_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif



#define THREAD_DATA_FULL 1
#define THREAD_DATA_BEAM 2


struct seqer_thread_data{
        struct fast_param_bag* ft_bag;
        struct fast_hmm_param* ft;
        struct tl_seq_buffer* sb;
        struct fhmm** fhmm;

        struct fhmm_dyn_mat* fmat;
        double** dyn;
        int info;
        //double** F_matrix;
        //double** B_matrix;
        //double** t;             /* for forward backward - to collect estimated */
        //double** e;             /* for forward backward - to collect estimated */
        //int num_seq;
        int thread_ID;
        int model_ID;
        int num_models;
        int num_threads;
        unsigned int seed;
        rk_state rndstate;
};
EXTERN int create_seqer_thread_data(struct seqer_thread_data*** t, int num_threads, int max_len, int K,rk_state* random);
//EXTERN struct seqer_thread_data** create_seqer_thread_data(int* num_threads, int max_len, int K,rk_state* random, int mode);
EXTERN int resize_seqer_thread_data(struct seqer_thread_data** td, int max_len, int K);
EXTERN int compare_seqer_thread_data(struct seqer_thread_data** a , struct seqer_thread_data** b, int num);

EXTERN void free_seqer_thread_data(struct seqer_thread_data** td);

#undef THREAD_DATA_IMPORT
#undef EXTERN

#endif
