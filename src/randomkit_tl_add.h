#ifndef RANDOMKIT_TL_ADD_H
#define RANDOMKIT_TL_ADD_H

#include "randomkit.h"


#ifdef RANDOMKIT_TL_ADD_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif



EXTERN int copy_rk_state(rk_state* source, rk_state* target);
EXTERN int compare_rk_state(rk_state* a, rk_state* b);

#undef RANDOMKIT_TL_ADD_IMPORT
#undef EXTERN


#endif
