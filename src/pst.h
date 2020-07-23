#ifndef PST_H
#define PST_H

//#include "tlseqio.h"

//

#include "inttypes.h"
#define MAX_PST_LEN 12

#ifdef PST_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct pst_node;
struct pst;

struct count_hash;

EXTERN int run_build_pst(struct pst** pst,float min_error, float gamma, struct count_hash* h);

//EXTERN int score_pst(const struct pst* pst, const uint8_t* seq,const int len, float* P_M, float* P_R);
EXTERN int score_pst(const struct pst* pst, const uint8_t* seq,const int len, float* score);

EXTERN int z_score_pst(const struct pst* p, int len, float score,double* z_score);

//EXTERN int score_pst(const struct pst* pst, const char* seq,const int len, float* P_M, float* P_R);
EXTERN void free_pst(struct pst* p);

#undef PST_IMPORT
#undef EXTERN

#endif
