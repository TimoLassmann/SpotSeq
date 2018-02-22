#ifndef IHMM_SEQ_H
#define IHMM_SEQ_H



#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"
#include <stdint.h>

struct ihmm_sequences{
        char** seq;
        float** u;
        int** label;
        int* len;
        int numseq;
        int max_len;
};

extern struct ihmm_sequences* create_ihmm_sequences(char** seq, int numseq);

extern void free_ihmm_sequences(struct ihmm_sequences* iseq);

extern int random_label_ihmm_sequences(struct ihmm_sequences* iseq, int k);

#endif
