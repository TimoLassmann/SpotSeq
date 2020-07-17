#ifndef FINITE_HMM_ALLOC_H
#define FINITE_HMM_ALLOC_H


#ifdef FINITE_HMM_ALLOC_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


struct fhmm;

struct fhmm_dyn_mat;

EXTERN struct fhmm* alloc_fhmm(void);
EXTERN void free_fhmm(struct fhmm* fhmm);


EXTERN int alloc_fhmm_dyn_mat(struct fhmm_dyn_mat** mat,int L,int K);
EXTERN int resize_fhmm_dyn_mat(struct fhmm_dyn_mat* dm,int new_len);
EXTERN int free_fhmm_dyn_mat(struct fhmm_dyn_mat* dm);


#undef FINITE_HMM_ALLOC_IMPORT
#undef EXTERN
#endif
