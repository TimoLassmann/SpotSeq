#ifndef FINITE_HMM_IO_H
#define FINITE_HMM_IO_H

#include "tlhdf5wrap.h"


#ifdef FINITE_HMM_IO_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


EXTERN struct fhmm*  read_fhmm_parameters(struct hdf5_data* hdf5_data, char* group);
EXTERN struct fhmm* read_best_fmodel(char* filename, int* best_model);
EXTERN int add_fhmm(struct hdf5_data* hdf5_data, struct fhmm* fhmm, char* group);
EXTERN int write_fhmm(char* filename, struct fhmm* fhmm);

EXTERN int write_searchfhmm(char* filename, struct fhmm* fhmm);
EXTERN int read_searchfhmm(char* filename,struct fhmm** ret);


#undef FINITE_HMM_IO_IMPORT
#undef EXTERN
#endif
