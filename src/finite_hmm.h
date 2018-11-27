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
        float r_score;
        int alloc_matrix_len;
        int K;
        int L;
};

extern int random_model_score(float* b, float* ret_score, uint8_t* a, int len, int expected_len);
extern int forward(struct fhmm* fhmm,float** matrix,float* ret_score, uint8_t* a, int len);
extern int backward(struct fhmm* fhmm,float** matrix, float* ret_score, uint8_t* a, int len);
extern int posterior_decoding(struct fhmm* fhmm,float** Fmatrix, float** Bmatrix,float score,uint8_t* a, int len,int* path);

//extern int read_fhmm_parameters(struct fhmm* fhmm, char* filename, char* model_name);
struct fhmm*  read_fhmm_parameters(struct hdf5_data* hdf5_data, char* group);

extern struct fhmm* alloc_fhmm(void);
extern int setup_model(struct fhmm* fhmm);
extern int convert_fhmm_scaled_to_prob(struct fhmm* fhmm);
extern struct fhmm* init_fhmm(char* filename);
extern int alloc_dyn_matrices(struct fhmm* fhmm);
extern int realloc_dyn_matrices(struct fhmm* fhmm,int new_len);
extern void free_fhmm(struct fhmm* fhmm);
#endif
