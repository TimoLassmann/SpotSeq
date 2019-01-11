#ifndef MOTIF_REFINEMENT_H
#define MOTIF_REFINEMENT_H



#include "ihmm_seq.h"
#include "tldevel.h"


int em_algorithm(double** counts,int W, int L, struct seq_buffer* sb);
#endif
