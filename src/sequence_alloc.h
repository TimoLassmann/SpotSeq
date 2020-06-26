
#ifndef SEQUENCE_ALLOC_H
#define SEQUENCE_ALLOC_H

#ifdef SEQUENCE_ALLOC_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define MAX_SEQUENCE_NAME_LEN 256

struct ihmm_sequence;
struct seq_buffer;

EXTERN int alloc_sequence_buffer(struct seq_buffer** seq_buf, int num_seq);
EXTERN int alloc_ihmm_seq(struct ihmm_sequence** is);
EXTERN int realloc_ihmm_seq(struct ihmm_sequence* sequence, int new_len);
EXTERN void free_ihmm_sequences(struct seq_buffer* sb);

EXTERN int add_multi_model_label_and_u(struct seq_buffer* sb,int num_models);

#undef SEQUENCE_ALLOC_IMPORT
#undef EXTERN


#endif
