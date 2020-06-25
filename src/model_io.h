#ifndef MODEL_IO_H
#define MODEL_IO_H


#ifdef MODEL_IO_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


EXTERN struct ihmm_model* read_best_imodel(char* filename, int* best_model);
EXTERN struct fhmm* read_best_fmodel(char* filename, int* best_model);

EXTERN struct model_bag* read_model_bag_hdf5(char* filename);
EXTERN int write_model_bag_hdf5(struct model_bag* bag, char* filename);
EXTERN int write_best_model(char* filename, int best_model);

#undef MODEL_IO_IMPORT
#undef EXTERN

#endif
