#ifndef MODEL_ALLOC_H
#define MODEL_ALLOC_H


#ifdef MODEL_ALLOC_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct model_bag;

EXTERN struct model_bag* alloc_model_bag(int L, int num_models, int max_states, int seed);

EXTERN void free_model_bag(struct model_bag* b);


EXTERN struct ihmm_model* alloc_ihmm_model(int maxK, int L, unsigned int seed);
EXTERN int resize_ihmm_model(struct ihmm_model* ihmm, int K);
EXTERN void free_ihmm_model(struct ihmm_model* ihmm);

#undef MODEL_ALLOC_IMPORT
#undef EXTERN

#endif
