
#ifndef SEQUENCE_ALLOC_H
#define SEQUENCE_ALLOC_H

#ifdef SEQUENCE_ALLOC_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define MAX_SEQUENCE_NAME_LEN 256


struct seq_ihmm_data;

struct tl_seq;
EXTERN int alloc_ihmm_seq_data(struct tl_seq* s, int num_models, int max_len);
EXTERN int free_ihmm_seq_data(struct seq_ihmm_data** data);

#undef SEQUENCE_ALLOC_IMPORT
#undef EXTERN


#endif
