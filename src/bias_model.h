#ifndef BIASMODEL_H
#define BIASMODEL_H


#include "finite_hmm_struct.h"

#ifdef BIASMODEL_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


EXTERN int build_bias_model(struct fhmm* model, struct fhmm** bias_model);


EXTERN int write_biashmm(char* filename, struct fhmm* fhmm);
EXTERN int read_biasfhmm(char* filename,struct fhmm** ret);

#undef BIASMODEL_IMPORT
#undef EXTERN

#endif
