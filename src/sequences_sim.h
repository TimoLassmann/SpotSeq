#ifndef SEQUENCES_SIM_H
#define SEQUENCES_SIM_H

#ifdef SEQUENCES_SIM_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct rng_state;
struct seq_buffer;

EXTERN int sim_sequences(int N, int L,int len,struct seq_buffer** seq_buf,struct rng_state* rng);

#undef SEQUENCES_SIM_IMPORT
#undef EXTERN


#endif
