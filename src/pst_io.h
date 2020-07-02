#ifndef PST_IO_H
#define PST_IO_H



#ifdef PST_IO_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct pst;

EXTERN int write_pst_hdf5(struct pst* p, char* filename);
EXTERN int read_pst_hdf5(struct pst** p, char* filename);


#undef PST_IO_IMPORT
#undef EXTERN


#endif
