#ifndef IHMM_SEQ_H
#define IHMM_SEQ_H



#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"
#include <stdint.h>
#include <ctype.h>

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
};




extern int random_label_ihmm_sequences(struct seq_buffer* sb, int k);

extern struct seq_buffer* create_ihmm_sequences_mem(char** seq, int numseq);
extern struct seq_buffer* load_sequences(char* in_filename);

extern void free_ihmm_sequences(struct seq_buffer* sb);



#endif
