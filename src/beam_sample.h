#ifndef BEAM_SAMPLE_H
#define BEAM_SAMPLE_H


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "distributions.h"
#include <math.h>
#include <float.h>
#include <stdint.h>

#include "tldevel.h"
#include "fast_hmm_param.h"
#include "ihmm_seq.h"
#include "model.h"
#include "global.h"

struct beam_parameters{
        struct ihmm_model* model;
        char** sequences;
        float alpha_a;
        float alpha_b;
        float gamma_a;
        float gamma_b;
        int numseq;
        int num_thread;
};

/* key Operations  */
extern int fill_fast_transitions(struct ihmm_model* model,struct fast_hmm_param* ft);
extern int add_state_from_fast_hmm_param(struct ihmm_model* ihmm,struct fast_hmm_param* ft);

#endif
