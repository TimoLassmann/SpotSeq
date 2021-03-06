#include "tllogsum.h"
#include "tlseqbuffer.h"
#include "run_score.h"

#include "sequence_struct.h"

#include "finite_hmm_alloc.h"
#include "finite_hmm_score.h"

#include "thread_data.h"

#include "model_struct.h"

int run_score_sequences(struct fhmm** fhmm, struct tl_seq_buffer* sb,struct seqer_thread_data** td, int n_model,int mode)
{
        //struct fhmm** container = NULL;
        int i;
        int num_threads;
        ASSERT(fhmm != NULL,"no model");
        ASSERT(sb != NULL, "no parameters");
        /* just to be 100% safe... */
        init_logsum();
        /* always a good idea */
        /* allocate dyn programming matrices.  */
        //RUN(realloc_dyn_matrices(fhmm, sb->max_len+1));
        //LOG_MSG("new len: %d states:%d", sb->max_len,fhmm->K);
        num_threads = td[0]->num_threads;

        //MMALLOC(container, sizeof(struct fhmm*) * 1);
        //container[0] = fhmm;
        for(i = 0; i < num_threads;i++){
                td[i]->thread_ID = i;
                td[i]->sb = sb;
                td[i]->fhmm = fhmm;
                td[i]->num_models = n_model;
                td[i]->info = mode;
        }

#ifdef HAVE_OPENMP
        omp_set_num_threads(num_threads);
#pragma omp parallel shared(td) private(i)
        {
#pragma omp for schedule(dynamic) nowait
#endif
                for(i = 0; i < num_threads;i++){
                        do_score_sequences(td[i]);
                }
#ifdef HAVE_OPENMP
        }
#endif
        return OK;
ERROR:
        return FAIL;
}

void* do_score_sequences(void* threadarg)
{
        struct seqer_thread_data *data;
        struct tl_seq* seq = NULL;
        struct seq_ihmm_data* d = NULL;
        struct fhmm_dyn_mat* m;
        double* s = NULL;
        int mode;
        int i;
        int j;
        int num_threads;
        int thread_id;
        int num_models;
        double fwd;
        double bias;
        double null;
        /* function pointer  */
        /* remember:
           return_type (*function_name)(arguments)
        */

        int (*func_ptr)(struct fhmm*, struct fhmm_dyn_mat*, uint8_t* , int,int , double*);

        data = (struct seqer_thread_data *) threadarg;

        num_threads = data->num_threads;
        thread_id = data->thread_ID;
        num_models = data->num_models;
        mode = data->info;

        if(mode == FHMM_SCORE_P_LODD){
                func_ptr = fhmm_score_p_lodd;
        }else if(mode == FHMM_SCORE_LODD){
                func_ptr = fhmm_score_lodd;
        }else if(mode == FHMM_SCORE_FULL){
                /* Run in full  more */
        }else{
                ERROR_MSG("Unknown function type!");
        }
//model_id = data->model_ID;
        m = data->fmat;
        //fhmm = data->fhmm;
        j = 0;
        for(i = 0; i < data->num_models;i++){

                if(j < data->fhmm[i]->K){
                        j = data->fhmm[i]->K;
                }
        }
        /* make sure we have enough memory  */

        if(m->alloc_matrix_len < data->sb->max_len || m->alloc_K < j){
                //LOG_MSG("have: %d %d", m->alloc_matrix_len, m->alloc_K);
                //LOG_MSG("want: %d %d", data->sb->max_len,  j);
                resize_fhmm_dyn_mat(m,
                                    MACRO_MAX(m->alloc_matrix_len, data->sb->max_len),
                                    MACRO_MAX(m->alloc_K,j)
                        );
        }

        if(mode == FHMM_SCORE_FULL){
                /* run forward on model, bias model and random  */
                for(i =0; i < data->sb->num_seq;i++){
                        if( i% num_threads == thread_id){
                                seq = data->sb->sequences[i];
                                s = seq->data;
                                fhmm_score_fwd(data->fhmm[0],m,seq->seq, seq->len,1, &fwd);
                                fhmm_score_fwd(data->fhmm[1],m,seq->seq, seq->len,1, &bias);
                                fhmm_score_null(data->fhmm[1],m,seq->seq, seq->len,1, &null);
                                s[0] =(fwd - null) / 0.69314718055994529;
                                s[1] =(fwd - bias) / 0.69314718055994529;
                                s[2] = esl_exp_surv(s[0], data->fhmm[0]->tau,data->fhmm[0]->lambda);
                                s[3] = esl_exp_surv(s[1], data->fhmm[0]->tau,data->fhmm[0]->lambda);
                        }
                }
                return NULL;
        }
        //LOG_MSG("Average sequence length: %d",expected_len);
        for(i =0; i < data->sb->num_seq;i++){
                if( i% num_threads == thread_id){
                        seq = data->sb->sequences[i];
                        s = seq->data;
                        for(j = 0; j < num_models;j++){
                                func_ptr(data->fhmm[j],m,seq->seq, seq->len,1, &s[j]);
                                //LOG_MSG("Thread: %d;model:%d Seq: %s %f %f",thread_id,model_id, data->sb->sequences[i]->name, f_score, logP);
                                //RUN(forward(fhmm, m, &f_score, seq->seq, seq->seq_len ));
                                //seq->score_arr[thread_id] = logP
                                //seq->score = f_score;
                        }
                }
        }
        return NULL;
ERROR:
        WARNING_MSG("Something really went wrong....");
        return NULL;
}




void* do_score_sequences_per_model(void* threadarg)
{
        struct seqer_thread_data *data;

        struct fhmm* fhmm = NULL;
        struct fhmm_dyn_mat* m = NULL;
        struct tl_seq* seq = NULL;
        struct seq_ihmm_data* d = NULL;
        int i;
        //int num_threads;
        int thread_id;          /* this is actually the model number */

        double f_score;
        double logP;
        //double r_score;
        data = (struct seqer_thread_data *) threadarg;

        //num_threads = data->num_threads;
        thread_id = data->thread_ID;
        fhmm = data->fhmm;
        m = data->fmat;

        //LOG_MSG("Average sequence length: %d",expected_len);

        for(i =0; i < data->sb->num_seq;i++){
                seq = data->sb->sequences[i];
                d = seq->data;
                //score_seq_fwd(fhmm,m,seq->seq, seq->len,1, &f_score, &logP);
                LOG_MSG("Seq: %s %f %f", data->sb->sequences[i]->name, f_score, logP);
                //RUN(forward(fhmm, m, &f_score, seq->seq, seq->seq_len ));
                d->score_arr[thread_id] = logP;
        }
        return NULL;
ERROR:
        return NULL;
}

void* do_label_sequences(void* threadarg)
{
        struct seqer_thread_data *data;
        struct fhmm* fhmm = NULL;
        struct tl_seq* seq  = NULL;
        //struct ihmm_sequence* seq = NULL;
        int i;
        int num_threads;
        int thread_id;
        struct fhmm_dyn_mat* m;
        float f_score;
        float b_score;
        data = (struct seqer_thread_data *) threadarg;

        num_threads = data->num_threads;
        thread_id = data->thread_ID;
        fhmm = data->fhmm;
        //LOG_MSG("Average sequence length: %d",expected_len);

        for(i =0; i < data->sb->num_seq;i++){
                if( i% num_threads == thread_id){
                        seq = data->sb->sequences[i];
                        RUN( forward(fhmm, m, &f_score, seq->seq, seq->len,1));
                        RUN(backward(fhmm, m, &b_score, seq->seq, seq->len,1));
                        //RUN(posterior_decoding(fhmm,data->F_matrix,data->B_matrix,f_score,seq->seq, seq->seq_len, seq->label));
                }
        }
        return NULL;
ERROR:
        return NULL;
}


int run_label_sequences(struct fhmm* fhmm, struct tl_seq_buffer* sb, int num_threads)
{

        //struct thr_pool* pool = NULL;
        struct seqer_thread_data** td = NULL;
        int i;
        ASSERT(fhmm != NULL,"no model");
        ASSERT(sb != NULL, "no parameters");

        init_logsum();
        /* start threadpool  */
        //if((pool = thr_pool_create(num_threads ,num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");


        /* allocate data for threads; */
        RUN(create_seqer_thread_data(&td,num_threads,(sb->max_len+2)  , fhmm->K ,NULL));// & sb->rndstate));

        /* score sequences  */

        /*for(i = 0; i < num_threads;i++){
                td[i]->sb = sb;
                td[i]->fhmm = fhmm;
                if(thr_pool_queue(pool, do_label_sequences,td[i]) == -1){
                        fprintf(stderr,"Adding job to queue failed.");
                }
        }
        thr_pool_wait(pool);
        */
        for(i = 0; i < num_threads;i++){
                td[i]->sb = sb;
                td[i]->fhmm = fhmm;
        }


#ifdef HAVE_OPENMP
        omp_set_num_threads(num_threads);
#pragma omp parallel shared(td) private(i)
        {
#pragma omp for schedule(dynamic) nowait
#endif
                for(i = 0; i < num_threads;i++){
                        do_label_sequences(td[i]);
                }
#ifdef HAVE_OPENMP
        }
#endif


        free_seqer_thread_data(td);
        //thr_pool_destroy(pool);

        return OK;
ERROR:
        free_seqer_thread_data(td);
        //thr_pool_destroy(pool);
        return FAIL;
}
