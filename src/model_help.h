#ifndef MODEL_HELP_H
#define MODEL_HELP_H

#ifdef MODEL_HELP_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

/* help functions */
EXTERN int print_counts(struct ihmm_model* ihmm);
EXTERN int print_model_parameters(struct ihmm_model* ihmm);

EXTERN int compare_model_bag(struct model_bag* a, struct model_bag* b);


#undef MODEL_HELP_IMPORT
#undef EXTERN

#endif
