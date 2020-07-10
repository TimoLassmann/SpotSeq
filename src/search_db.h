#ifndef SEARCH_DB_H
#define SEARCH_DB_H

#ifdef SEARCH_DB_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


EXTERN int build_sequence_database(char* filename, char* out,int seed);

#undef SEARCH_DB_IMPORT
#undef EXTERN

#endif
