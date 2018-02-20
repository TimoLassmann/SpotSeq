#ifndef FAST_HMM_PARAM_TEST_FUNCTIONS_H
#define FAST_HMM_PARAM_TEST_FUNCTIONS_H


#include "fast_hmm_param.h"

extern int print_fast_hmm_params(struct fast_hmm_param* ft);
extern int fill_with_random_transitions(struct fast_hmm_param* ft, int k);

#endif
