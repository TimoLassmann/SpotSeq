#ifndef EMIT_RANDOM_H
#define EMIT_RANDOM_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "global.h"
#include "tldevel.h"
#include "ihmm_seq.h"
#include "distributions.h"


#include "finite_hmm.h"
extern struct seq_buffer* emit_sequences_from_fhmm_model(struct fhmm* fhmm, int num, unsigned long seed);
extern struct seq_buffer* emit_sequences_from_random_model(struct seq_buffer* sb_in, int num,unsigned long seed);
#endif
