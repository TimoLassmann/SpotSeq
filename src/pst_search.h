#ifndef PST_SEARCH_H
#define PST_SEARCH_H


#ifdef PST_SEARCH_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct pst;

EXTERN int search_db(struct pst* p, char* filename, double thres);

#undef PST_SEARCH_IMPORT
#undef EXTERN


#endif
