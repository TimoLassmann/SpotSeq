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

struct ihmm_sequence{
        uint8_t* seq;
        float** u_arr;
        int** label_arr;
        float* u;
        int* label;
        char* name;
        float score;
        float r_score;
        int malloc_len;
        int seq_len;
};

struct seq_buffer{
        struct ihmm_sequence** sequences;
        float* background;

        rk_state rndstate;
        int seed;
        int malloc_num;
        int num_seq;
        int org_num_seq;
        int max_len;
        int L;
};

extern struct ihmm_sequence* alloc_ihmm_seq(void);

extern int add_multi_model_label_and_u(struct seq_buffer* sb,int num_models);


extern int realloc_ihmm_seq(struct ihmm_sequence* sequence);

extern struct seq_buffer* get_sequences_from_hdf5_model(char* filename);
extern int add_sequences_to_hdf5_model(char* filename,struct seq_buffer* sb);

extern int random_label_ihmm_sequences(struct seq_buffer* sb, int k,float alpha);

extern int shuffle_sequences_in_buffer(struct seq_buffer* sb);

extern struct seq_buffer* concatenate_sequences(struct seq_buffer* sb);
extern int dirichlet_emission_label_ihmm_sequences(struct seq_buffer* sb, int k, float alpha);
extern int label_ihmm_sequences_based_on_guess_hmm(struct seq_buffer* sb, int k, float alpha);

extern int print_states_per_sequence(struct seq_buffer* sb);
extern struct seq_buffer* create_ihmm_sequences_mem(char** seq, int numseq);
extern struct seq_buffer* load_sequences(char* in_filename);

extern int add_reverse_complement_sequences_to_buffer(struct seq_buffer* sb);


extern int print_labelled_ihmm_buffer(struct seq_buffer* sb);

extern void free_ihmm_sequences(struct seq_buffer* sb);
extern int write_ihmm_sequences_fasta(struct seq_buffer* sb, char* filename);
extern int write_ihmm_sequences(struct seq_buffer* sb, char* filename, char* comment);
struct seq_buffer* load_ihmm_sequences(char* in_filename);
#endif
