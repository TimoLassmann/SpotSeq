#ifndef HMM_SCORE_H
#define HMM_SCORE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "tldevel.h"
#include "tlhdf5wrap.h"
#include "global.h"


#include "finite_hmm_struct.h"

extern int configure_target_len(struct fhmm* fhmm,int len,  int multihit);

extern int remove_state_for_ploting(struct fhmm*fhmm, int state);


extern int calculate_BIC( struct fhmm* fhmm, double ML, double data,double* BIC);
//extern int random_model_score(double* b, double* ret_score, uint8_t* a, int len, int expected_len);
extern int random_model_score(int len,float* ret_score);
extern int forward(struct fhmm* fhmm , struct fhmm_dyn_mat* m, float* ret_score, uint8_t* a, int len);
extern int backward(struct fhmm* fhmm,struct fhmm_dyn_mat* m , float* ret_score, uint8_t* a, int len);
//extern int backward(struct fhmm* fhmm,float** matrix, float* ret_score, uint8_t* a, int len);
//extern int backward(struct fhmm* fhmm,double** matrix, double* ret_score, uint8_t* a, int len);
extern int posterior_decoding(struct fhmm* fhmm,double** Fmatrix, double** Bmatrix,double score,uint8_t* a, int len,int* path);



extern int setup_model(struct fhmm* fhmm);
extern int convert_fhmm_scaled_to_prob(struct fhmm* fhmm);
extern int convert_fhmm_log_to_prob_for_sampling(struct fhmm* fhmm);
#endif
