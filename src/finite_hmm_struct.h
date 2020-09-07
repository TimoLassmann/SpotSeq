#ifndef FINITE_HMM_STRUCT_H
#define FINITE_HMM_STRUCT_H

#define S_STATE 0
#define N_STATE 1
#define B_STATE 2
#define E_STATE 3
#define C_STATE 4
#define J_STATE 5
#define T_STATE 6

struct fhmm{
        float** F_matrix;
        float** B_matrix;

        float** F_NBECJ;
        float** B_NBECJ;
        float** e;
        float** t;
        int** tindex;
        float* m_comp_back; /* Equivalent (hopefully to compo in HMMER - see Biased composition filter.) */
        float* background;
        float f_score;
        float b_score;
        float r_score;
        double H;
        double lambda;
        double tau;
        int alloc_matrix_len;
        int alloc_K;
        int K;
        int L;
        /* int config_len; */

        /* float tSN; */
        /* float tNN; */
        /* float tNB; */

        /* float tBX; */
        /* float tXE; */

        /* float tEC; */
        /* float tCC; */
        /* float tCT; */

        /* float tEJ; */
        /* float tJJ; */
        /* float tJB; */
};

struct fhmm_dyn_mat{
        float** F_matrix;
        float** B_matrix;

        float** F_NBECJ;
        float** B_NBECJ;
        int alloc_matrix_len;
        int alloc_K;
};

#endif
