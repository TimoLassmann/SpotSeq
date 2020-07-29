#ifndef SEARCH_DB_H
#define SEARCH_DB_H

#include <inttypes.h>

#ifdef SEARCH_DB_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct hdf_seq_store{
        uint8_t* seq;
        uint32_t* len;
        char** seq_names;
        int num_seq;
        int alloc_num_seq;
        int seq_buf_size;
        int chunk_number;
};



EXTERN int build_sequence_database(char* filename, char* out,int seed);

EXTERN int read_hdf_seq_store_chunk(struct hdf_seq_store** hs, char* filename);

EXTERN void free_hdf_seq_store(struct hdf_seq_store* h);

#undef SEARCH_DB_IMPORT
#undef EXTERN

#endif
