#ifndef MODEL_H
#define MODEL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "global.h"
#include "tldevel.h"
#include "distributions.h"

#include "finite_hmm.h"
#include <math.h>
#include <stdint.h>





struct seq_buffer;          /* forward declaration  */

struct hdf5_data;               /* forward declaration`` */


struct ihmm_model{
        double** transition_counts;
        double** emission_counts;
        double* beta;
        unsigned int seed;
        rk_state rndstate;
        double gamma;
        double alpha;
        double alpha_a;
        double alpha_b;
        double gamma_a;
        double gamma_b;
        double log_likelihood;
        double alpha_limit;
        double gamma_limit;
        int num_states;         /* this excludes the start and stop states (0,1) */
        int alloc_num_states;
        int L;
        int training_iterations;
};

struct model_bag{
        struct ihmm_model** models;
        struct fhmm** finite_models;
        double* min_u;
        int best_model;
        int max_num_states;
        int num_models;
        unsigned int seed;      /* Starting value */
        rk_state rndstate;      /* main seed */
};

/* Housekeeping */

extern struct model_bag* alloc_model_bag(int* num_state_array, int L, int num_models, rk_state* rndstate);
extern void free_model_bag(struct model_bag* b);


extern struct ihmm_model* alloc_ihmm_model(int K, int L, unsigned int seed);
extern int resize_ihmm_model(struct ihmm_model* ihmm, int K);
extern void free_ihmm_model(struct ihmm_model* ihmm);

/* Model IO */
extern struct ihmm_model* read_best_imodel(char* filename, int* best_model);
extern struct fhmm* read_best_fmodel(char* filename, int* best_model);

extern struct model_bag* read_model_bag_hdf5(char* filename);
extern int write_model_bag_hdf5(struct model_bag* bag, char* filename);
//extern int add_fhmm(char* filename,double** transition,double** emission, int N, int L);
extern int add_fhmm(struct hdf5_data* hdf5_data, struct fhmm* fhmm, char* group);
extern int add_background_emission(char* filename,double* background,int L);
extern int add_annotation( char* filename, char* name, char* value);

//extern int write_model_hdf5(struct ihmm_model* model, char* filename);
//struct ihmm_model* read_model_hdf5(char* filename);
extern int write_model(struct ihmm_model* model, char* filename);
extern struct ihmm_model* read_model( char* filename);

extern int write_best_model(char* filename, int best_model);

/* Write RNG states in threads to file to ensure reproducibility.... */
struct wims_thread_data** read_thread_data_to_hdf5(char* filename);
int write_thread_data_to_hdf5(char* filename,struct wims_thread_data** td,int num_threads,int max_len,int max_K);

/* Initialize number of states.  */
extern int inititalize_model(struct ihmm_model* model, struct seq_buffer* sb, int K);
/* Fill counts from sequences  */
extern int clear_counts(struct ihmm_model* ihmm);
extern int fill_counts(struct ihmm_model* ihmm, struct seq_buffer* sb, int model_index);
extern int add_pseudocounts_emission(struct ihmm_model* model, double* background, double alpha);
//extern int remove_unused_states_labels(struct ihmm_model* ihmm, struct seq_buffer* sb);
extern int remove_unused_states_labels(struct ihmm_model* ihmm, struct seq_buffer* sb, int model_index);


/* set hyperparameters  */
extern int set_model_hyper_parameters(struct model_bag* b, double alpha, double gamma);
/* re-estimate hyper parameters */
extern int iHmmHyperSample(struct ihmm_model* model, int iterations);

/* help functions */
extern int print_counts(struct ihmm_model* ihmm);
extern int print_model_parameters(struct ihmm_model* ihmm);


#endif
