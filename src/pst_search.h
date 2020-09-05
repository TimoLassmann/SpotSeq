#ifndef PST_SEARCH_H
#define PST_SEARCH_H

#include "tlseqbuffer.h"

#ifdef PST_SEARCH_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct pst;

EXTERN int search_db(struct pst* p, char* filename, double thres,struct tl_seq_buffer** hits, uint64_t* db_size);
//EXTERN int search_db(struct pst* p, char* filename, double thres);
EXTERN int search_db_hdf5(struct pst* p, char* filename, double thres);

#undef PST_SEARCH_IMPORT
#undef EXTERN


#endif
