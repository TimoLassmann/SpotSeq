#ifndef PST_BUILD_H
#define PST_BUILD_H



#ifdef PST_BUILD
#define EXTERN
#else
#define EXTERN extern
#endif
EXTERN int create_pst_model(struct rng_state* rng,struct tl_seq_buffer* sb, char* train_seq,char* seq_db,char* out_model, double p_min, double gamma, double z_thres);

//EXTERN int create_pst_model(struct rng_state* rng,char* ihmm_model_file, char* train_seq,char* seq_db,char* out_model, double p_min, double gamma, double z_thres);

#undef PST_BUILD
#undef EXTERN
#endif
