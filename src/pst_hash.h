#ifndef PST_HASH_H
#define PST_HASH_H

#include "khash.h"

#include "tlseqio.h"

#ifdef PST_HASH_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif



KHASH_MAP_INIT_INT64(exact, int)

struct count_hash{
        khash_t(exact) * hash;
        uint64_t* mask;
        int L;
        int len;
};

EXTERN int fill_exact_hash(struct count_hash** hash, struct tl_seq_buffer* sb);
EXTERN void free_exact_hash(struct count_hash* hash);



#undef PST_HASH_IMPORT

#undef EXTERN



#endif
