#ifndef PST_CALIBRATE_H
#define PST_CALIBRATE_H

#ifdef PST_CALIBRATE_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct pst;

EXTERN int calibrate_pst(struct pst* p, char* filename,double threshold);

#undef PST_CALIBRATE_IMPORT
#undef EXTERN

#endif
