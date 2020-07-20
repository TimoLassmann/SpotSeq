#ifndef FINITE_HMM_PLOT_H
#define FINITE_HMM_PLOT_H



#include "finite_hmm_struct.h"

#ifdef FINITE_HMM_PLOT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int plot_finite_hmm_dot(struct fhmm* fhmm,char* filename,float thres);

#undef FINITE_HMM_PLOT_IMPORT
#undef EXTERN

#endif
