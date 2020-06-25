#ifndef MODEL_CORE_H
#define MODEL_CORE_H

#ifdef MODEL_CORE_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif



struct seq_buffer;          /* forward declaration  */

//struct hdf5_data;               /* forward declaration`` */

#include "model_struct.h"



/* Initialize number of states.  */
EXTERN int inititalize_model(struct seq_buffer* sb, int K);
/* Fill counts from sequences  */
EXTERN int clear_counts(struct ihmm_model* ihmm);
EXTERN int fill_counts(struct ihmm_model* ihmm, struct seq_buffer* sb, int model_index);


EXTERN int add_pseudocounts_emission(struct ihmm_model* model, double alpha);
//extern int remove_unused_states_labels(struct ihmm_model* ihmm, struct seq_buffer* sb);
EXTERN int remove_unused_states_labels(struct ihmm_model* ihmm, struct seq_buffer* sb, int model_index);


/* set hyperparameters  */
EXTERN int set_model_hyper_parameters(struct model_bag* b, double alpha, double gamma);
/* re-estimate hyper parameters */
EXTERN int iHmmHyperSample(struct ihmm_model* model, int iterations);

#undef MODEL_CORE_IMPORT
#undef EXTERN
#endif
