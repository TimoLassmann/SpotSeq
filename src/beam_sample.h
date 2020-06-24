#ifndef BEAM_SAMPLE_H
#define BEAM_SAMPLE_H

#ifdef BEAM_SAMPLE_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct model_bag;
struct seqer_thread_data;
struct fast_param_bag;

EXTERN int run_beam_sampling(struct model_bag* model_bag, struct fast_param_bag*ft_bag, struct seq_buffer* sb,struct seqer_thread_data** td, int iterations, int num_threads);


#undef BEAM_SAMPLE_IMPORT
#undef EXTERN




/* main runner function  */


#endif
