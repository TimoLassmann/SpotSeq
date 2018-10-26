#ifndef LABEL_SUFFIX_ARRAY_H
#define LABEL_SUFFIX_ARRAY_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"

#include "ihmm_seq.h"

struct lcs{
        int* str;
        int seq_num;
        int lcp;
        int pos;
};


struct sa{
        struct lcs** lcs;
        int len;
        int alloc_len;
};

extern int search_sa(struct sa* sa, int*p, int len, int* start, int* stop);


extern struct sa* build_sa(struct seq_buffer* sb);
extern void free_sa(struct sa* sa);

#endif
