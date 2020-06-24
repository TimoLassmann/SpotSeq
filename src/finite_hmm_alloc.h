#ifndef FINITE_HMM_ALLOC_H
#define FINITE_HMM_ALLOC_H


#ifdef FINITE_HMM_ALLOC_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


struct fhmm;

EXTERN struct fhmm* alloc_fhmm(void);
EXTERN void free_fhmm(struct fhmm* fhmm);

#undef FINITE_HMM_ALLOC_IMPORT
#undef EXTERN
#endif
