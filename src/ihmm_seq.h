#ifndef IHMM_SEQ_H
#define IHMM_SEQ_H



#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"

#include "global.h"
#include "distributions.h"

#include <stdint.h>
#include <ctype.h>
#include "hdf5_glue.h"


#define BLOCK_LEN 70

#define IHMM_SEQ_READ_ALL 0

#define IHMM_SEQ_READ_ONLY_SEQ 1

struct ihmm_sequence{
        uint8_t* seq;
        uint8_t* has_path;
        double** u_arr;
        //double* score_arr;
        int** label_arr;
        int** tmp_label_arr;
        double* u;
        int* label;
        char* name;
        double score;
        double r_score;
        int malloc_len;
        int seq_len;
};

struct seq_buffer{
        struct ihmm_sequence** sequences;
        double* background;

        rk_state rndstate;
        int seed;
        int malloc_num;
        int num_seq;
        int org_num_seq;
        int max_len;
        int L;
};

extern struct ihmm_sequence* alloc_ihmm_seq(void);

extern int translate_internal_to_DNA(struct seq_buffer* sb);


extern int translate_internal_to_PROTEIN(struct seq_buffer* sb);


extern int add_multi_model_label_and_u(struct seq_buffer* sb,int num_models);
extern int check_labels(struct seq_buffer* sb, int num_models);

extern int get_res_counts(struct seq_buffer* sb, double* counts);


extern int realloc_ihmm_seq(struct ihmm_sequence* sequence);

extern struct seq_buffer* get_sequences_from_hdf5_model(char* filename,int mode);
extern int add_sequences_to_hdf5_model(char* filename,struct seq_buffer* sb, int model_index);
extern int random_label_based_on_multiple_models(struct seq_buffer* sb, int K, int model_index, rk_state* random);
extern int random_label_ihmm_sequences(struct seq_buffer* sb, int k,double alpha);

extern int shuffle_sequences_in_buffer(struct seq_buffer* sb);

extern struct seq_buffer* concatenate_sequences(struct seq_buffer* sb);
extern int dirichlet_emission_label_ihmm_sequences(struct seq_buffer* sb, int k, double alpha);
extern int label_ihmm_sequences_based_on_guess_hmm(struct seq_buffer* sb, int k, double alpha);

extern int print_states_per_sequence(struct seq_buffer* sb);
extern struct seq_buffer* create_ihmm_sequences_mem(char** seq, int numseq,rk_state* rndstate);
extern struct seq_buffer* load_sequences(char* in_filename,rk_state* rndstate);

extern int add_reverse_complement_sequences_to_buffer(struct seq_buffer* sb);


extern int print_labelled_ihmm_buffer(struct seq_buffer* sb, rk_state* rndstate);

extern void free_ihmm_sequences(struct seq_buffer* sb);
extern int write_ihmm_sequences_fasta(struct seq_buffer* sb, char* filename, rk_state* rndstate);
extern int write_ihmm_sequences(struct seq_buffer* sb, char* filename, char* comment,rk_state* rndstate);
struct seq_buffer* load_ihmm_sequences(char* in_filename, rk_state* rndstate);
#endif
