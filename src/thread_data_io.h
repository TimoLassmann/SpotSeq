#ifndef THREAD_DATA_IO_H
#define THREAD_DATA_IO_H

#include "thread_data.h"

#include "tlhdf5wrap.h"

#ifdef THREAD_DATA_IO_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN struct seqer_thread_data** read_thread_data_to_hdf5(char* filename);
EXTERN int write_thread_data_to_hdf5(char* filename,struct seqer_thread_data** td,int num_threads,int max_len,int max_K);


#undef THREAD_DATA_IO_IMPORT
#undef EXTERN


#endif
