#ifndef HMM_CONVERSION_H
#define HMM_CONVERSION_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"

#include "model.h"
#include "finite_hmm.h"

#include "fast_hmm_param.h"

#include "ihmm_seq.h"

/* Move data from model into fast transition data structure */
extern int fill_fast_transitions(struct ihmm_model* model,struct fast_hmm_param* ft);
extern int fill_fast_transitions_only_matrices(struct ihmm_model* model,struct fast_hmm_param* ft);


/* Move data from h5 file to finite hmm model */
extern int run_build_fhmm_file(char* h5file);
extern struct fhmm* build_finite_hmm_from_infinite_hmm(struct ihmm_model* model);
/* Get background emissions from either sequence buffer or ihmm model */
extern int fill_background_emission(struct fast_hmm_param*ft,struct seq_buffer* sb);
extern int fill_background_emission_from_model(struct fast_hmm_param*ft, struct ihmm_model* model);

#endif
