#ifndef MODEL_H
#define MODEL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "global.h"
#include "tldevel.h"
#include "distributions.h"

#include "finite_hmm.h"
#include <math.h>
#include <stdint.h>


struct seq_buffer;          /* forward declaration  */

//struct hdf5_data;               /* forward declaration`` */

#include "model_struct.h"


/* Initialize number of states.  */
extern int inititalize_model(struct ihmm_model* model, struct seq_buffer* sb, int K);
/* Fill counts from sequences  */
extern int clear_counts(struct ihmm_model* ihmm);
extern int fill_counts(struct ihmm_model* ihmm, struct seq_buffer* sb, int model_index);
extern int add_pseudocounts_emission(struct ihmm_model* model, double alpha);
//extern int remove_unused_states_labels(struct ihmm_model* ihmm, struct seq_buffer* sb);
extern int remove_unused_states_labels(struct ihmm_model* ihmm, struct seq_buffer* sb, int model_index);


/* set hyperparameters  */
extern int set_model_hyper_parameters(struct model_bag* b, double alpha, double gamma);
/* re-estimate hyper parameters */
extern int iHmmHyperSample(struct ihmm_model* model, int iterations);


#endif
