#ifndef SEQUENCE_STRUCT_H
#define SEQUENCE_STRUCT_H

#include <inttypes.h>

#include "tlalphabet.h"

struct ihmm_sequence{
        uint8_t* seq;
        uint8_t* has_path;
        double** u_arr;
        double* score_arr;
        uint16_t** label_arr;
        uint16_t** tmp_label_arr;
        double* u;
        int* label;
        char* name;
        double score;
        double r_score;
        int malloc_len;
        int seq_len;
};

struct ihmm_seq_data{
        
};

struct seq_buffer{
        struct ihmm_sequence** sequences;
        struct alphabet* alphabet;
        //int* num_state_arr;
        int malloc_num;
        int num_seq;
        int org_num_seq;
        int max_len;
        int L;
};

#endif
