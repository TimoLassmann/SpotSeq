#ifndef HMM_SCORE_H
#define HMM_SCORE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "tldevel.h"
#include "tlhdf5wrap.h"
#include "global.h"


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

extern int remove_state_for_ploting(struct fhmm*fhmm, int state);


extern int calculate_BIC( struct fhmm* fhmm, double ML, double data,double* BIC);
extern int random_model_score(double* b, double* ret_score, uint8_t* a, int len, int expected_len);
extern int forward(struct fhmm* fhmm,double** matrix,double* ret_score, uint8_t* a, int len);
extern int backward(struct fhmm* fhmm,double** matrix, double* ret_score, uint8_t* a, int len);
extern int posterior_decoding(struct fhmm* fhmm,double** Fmatrix, double** Bmatrix,double score,uint8_t* a, int len,int* path);

//extern int read_fhmm_parameters(struct fhmm* fhmm, char* filename, char* model_name);
struct fhmm*  read_fhmm_parameters(struct hdf5_data* hdf5_data, char* group);

extern struct fhmm* alloc_fhmm(void);
extern int setup_model(struct fhmm* fhmm);
extern int convert_fhmm_scaled_to_prob(struct fhmm* fhmm);
extern int convert_fhmm_log_to_prob_for_sampling(struct fhmm* fhmm);
extern struct fhmm* init_fhmm(char* filename);
extern int alloc_dyn_matrices(struct fhmm* fhmm);
extern int realloc_dyn_matrices(struct fhmm* fhmm,int new_len);
extern void free_fhmm(struct fhmm* fhmm);
#endif
