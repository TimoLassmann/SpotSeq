#ifndef INIT_SEQ_LABEL_H
#define INIT_SEQ_LABEL_H

#include "run_score.h"
#include "model.h"
#include "ihmm_seq.h"
#include "finite_hmm.h"



extern int label_seq_based_on_random_fhmm(struct seq_buffer* sb, int k, double alpha);

#endif
