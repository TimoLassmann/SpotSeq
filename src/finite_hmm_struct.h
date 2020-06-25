#ifndef FINITE_HMM_STRUCT_H
#define FINITE_HMM_STRUCT_H

struct fhmm{
        double** F_matrix;
        double** B_matrix;
        double** e;
        double** t;
        int** tindex;
        double* background;
        double f_score;
        double b_score;
        double r_score;
        int alloc_matrix_len;
        int alloc_K;
        int K;
        int L;
};

#endif
