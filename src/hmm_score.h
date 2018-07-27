#ifndef HMM_SCORE_H
#define HMM_SCORE_H



#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "tldevel.h"


#include "global.h"
#include "hdf5_glue.h"

struct fhmm{
        float** F_matrix;
        float** B_matrix;
        float** e;
        float** t;
        int** tindex;
        float* background;
        float f_score;
        float b_score;
        int alloc_matrix_len;
        int K;
        int L;
};



extern struct fhmm* init_fhmm(char* filename);
extern void free_fhmm(struct fhmm* fhmm);
#endif

