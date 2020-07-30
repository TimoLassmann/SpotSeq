#ifndef MOTIF_REFINEMENT_H
#define MOTIF_REFINEMENT_H



#include "tldevel.h"
#include "sequence_struct.h"

int em_algorithm(double** counts,int W, int L, struct seq_buffer* sb);
#endif
