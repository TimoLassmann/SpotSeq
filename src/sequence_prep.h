#ifndef SEQUENCE_PREP_H
#define SEQUENCE_PREP_H

#ifdef SEQUENCE_PREP_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct seq_buffer;
struct rng_state;


EXTERN int prep_sequences(struct seq_buffer* sb, struct rng_state* rng, int num_models,int num_states, double sigma);
#undef SEQUENCE_PREP_IMPORT
#undef EXTERN

#endif
