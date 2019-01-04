#ifndef MOTIF_REFINEMENT_H
#define MOTIF_REFINEMENT_H



#include "ihmm_seq.h"
#include "tldevel.h"

struct motif{
        double** freq_matrix;
        double** count_matrix;
        double* background_freq;
        double* background_counts;
        int W;
        int L;
        double log_likelihood;
};




int em_algorithm(double** counts,int W, int L, struct seq_buffer* sb);
#endif
