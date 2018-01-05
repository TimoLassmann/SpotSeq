
#ifndef ihmm_header


#define ihmm_header

#include "distributions.h"

#define iHMM_START_STATE 0
#define iHMM_STOP_STATE 1

struct iHMM_model{
	
	//model parameters
	rk_state rndstate;
	double* tmp;
	
	//priors
	int soft;
	float alpha0_a;
	float alpha0_b;
	
	float gamma_a;
	float gamma_b;
	int expected_K;
	
	int L; // alphabet length....
	int K; //number of states;
	int infinityghost; //  is the index beyond the allocated states...
	int malloced_states;
	
	
	//model parameters
	float alpha0;
	float gamma;
	
	float collect_alpha;
	float collect_gamma;
	
	float* beta;
	//float* collect_beta;
	
	int* prior_counts;
	int** trans_counts;
	int** emit_counts;
	float* back; 
	
	float* prior;
	
	
	
	// misc structs used internally
	float** transition;
	float** emission;
	
	float* w;
	float* p;
	float* s;
	
	int** M;
	int* sumM;
	int* sumN;
	
	//run info
	
	int numb;
	int nums;
	int numi;
};

extern struct iHMM_model* init_iHMM_model(void);

extern int start_iHMM_model(struct iHMM_model* model, int num_states);

extern int free_iHMM_model(struct iHMM_model* model);

extern int run_make_ihmm(struct iHMM_model* model, char** sequences,int numseq);
extern int particle_gibbs_with_ancestors_controller(struct iHMM_model* model,char** sequences,int numseq);

#endif