#ifndef PST_H
#define PST_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "global.h"
#include "tldevel.h"
#include "rbtree.h"

#define MAX_PST_LEN 25

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

struct lcs_count{
        int* counts;
        int unique;
        int len;
        int cur_longest;
        int min_seq;
};


struct motif_struct{
        int* state_list;
        int len;
        int start_in_sa;
        int end_in_sa;
};
extern void free_sa(struct sa* sa);

#endif
