#ifndef FINITE_HMM_STATS_H
#define FINITE_HMM_STATS_H

#include "finite_hmm_struct.h"

#ifdef FINITE_HMM_STATS_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif



EXTERN int fhmm_calibrate(struct fhmm* fhmm,struct fhmm_dyn_mat* dm, int seed);

EXTERN double esl_exp_surv(double x, double mu, double lambda);


#undef FINITE_HMM_STATS_IMPORT

#undef EXTERN

#endif
