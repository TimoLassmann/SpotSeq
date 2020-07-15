#ifndef GLOBAL_H
#define GLOBAL_H

#define START_STATE 0
#define END_STATE 1

#define N_STATE 2
#define B_STATE 3
#define E_STATE 4
#define C_STATE 5
#define J_STATE 6


#define EMISSION_H 0.3

#define GB 1073741824LL

#include <float.h>
/* safe logistic function  */
#define LOGISTIC_FLT(x) ((x)) >= log(FLT_MAX) ? 1.0f : (expf((x)) /(1.0f + expf((x))))

#define ALPHABET_DNA 4
#define ALPHABET_APPROXDNA 16
#define ALPHABET_PROTEIN 20

#define MAX_NUM_STATES 1000

#define IHMM_PARAM_PLACEHOLDER -9999.99

extern int double_logsum_init(void);
extern double double_logsum (double a, double b);

#endif
