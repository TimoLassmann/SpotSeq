#ifndef RANDOMKIT_IO_H
#define RANDOMKIT_IO_H

#include "randomkit.h"
#include "tlhdf5wrap.h"

#ifdef RANDOMKIT_IO_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


EXTERN int add_RNG_state(struct hdf5_data* hdf5_data, char* group,rk_state* a);

EXTERN int read_RNG_state(struct hdf5_data* hdf5_data, char* group,rk_state* a);

#undef RANDOMKIT_IO_IMPORT
#undef EXTERN


#endif
