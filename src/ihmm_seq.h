#ifndef IHMM_SEQ_H
#define IHMM_SEQ_H



#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"


#include "distributions.h"

#include <stdint.h>
#include <ctype.h>


#define ALPHABET_DNA 5
#define ALPHABET_APPROXDNA 16
#define ALPHABET_PROTEIN 20

#define BLOCK_LEN 70

struct ihmm_sequence{
        uint8_t* seq;
        float* u;
        int* label; 
        char* name;
        int malloc_len;
        int seq_len;
};
 
struct seq_buffer{
        struct ihmm_sequence** sequences;
        int malloc_num;
        int num_seq;
        int max_len;
        int L;
};




extern int random_label_ihmm_sequences(struct seq_buffer* sb, int k);

extern int dirichlet_emission_label_ihmm_sequences(struct seq_buffer* sb, int k, float alpha);


extern struct seq_buffer* create_ihmm_sequences_mem(char** seq, int numseq);
extern struct seq_buffer* load_sequences(char* in_filename);

extern int print_labelled_ihmm_buffer(struct seq_buffer* sb);

extern void free_ihmm_sequences(struct seq_buffer* sb);

extern int write_ihmm_sequences(struct seq_buffer* sb, char* filename, char* comment);
struct seq_buffer* load_ihmm_sequences(char* in_filename);
#endif
