#ifndef GLOBAL_H
#define GLOBAL_H

#define IHMM_START_STATE 0
#define IHMM_END_STATE 1

#define EMISSION_H 0.3

#define GB 1073741824LL

#include <float.h>
/* safe logistic function  */
#define LOGISTIC_FLT(x) ((x)) >= log(FLT_MAX) ? 1.0f : (expf((x)) /(1.0f + expf((x))))


#endif
