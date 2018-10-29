
#include "run_score.h"


int run_score_sequences(struct fhmm* fhmm, struct seq_buffer* sb, int num_threads)
{
        struct thr_pool* pool = NULL;
        struct spotseq_thread_data** td = NULL;
        int i;
        ASSERT(fhmm != NULL,"no model");
        ASSERT(sb != NULL, "no parameters");

        /* start threadpool  */
        if((pool = thr_pool_create(num_threads ,num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");


        /* allocate dyn programming matrices.  */
        RUN(realloc_dyn_matrices(fhmm, sb->max_len+1));

        /* allocate data for threads; */
        RUNP(td = create_spotseq_thread_data(&num_threads,(sb->max_len+2)  , fhmm->K ));

        /* score sequences  */

        for(i = 0; i < num_threads;i++){
                td[i]->sb = sb;
                td[i]->fhmm = fhmm;
                if(thr_pool_queue(pool, do_score_sequences,td[i]) == -1){
                        fprintf(stderr,"Adding job to queue failed.");
                }
        }
        thr_pool_wait(pool);

        free_spotseq_thread_data(td,num_threads);
        thr_pool_destroy(pool);

        return OK;
ERROR:
        free_spotseq_thread_data(td,num_threads);
        thr_pool_destroy(pool);
        return FAIL;
}

int run_label_sequences(struct fhmm* fhmm, struct seq_buffer* sb, int num_threads)
{

        struct thr_pool* pool = NULL;
        struct spotseq_thread_data** td = NULL;
        int i;
        ASSERT(fhmm != NULL,"no model");
        ASSERT(sb != NULL, "no parameters");

        /* start threadpool  */
        if((pool = thr_pool_create(num_threads ,num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");


        /* allocate data for threads; */
        RUNP(td = create_spotseq_thread_data(&num_threads,(sb->max_len+2)  , fhmm->K ));

        /* score sequences  */

        for(i = 0; i < num_threads;i++){
                td[i]->sb = sb;
                td[i]->fhmm = fhmm;
                if(thr_pool_queue(pool, do_label_sequences,td[i]) == -1){
                        fprintf(stderr,"Adding job to queue failed.");
                }
        }
        thr_pool_wait(pool);

        free_spotseq_thread_data(td,num_threads);
        thr_pool_destroy(pool);

        return OK;
ERROR:
        free_spotseq_thread_data(td,num_threads);
        thr_pool_destroy(pool);
        return FAIL;
}



void* do_score_sequences(void* threadarg)
{
        struct spotseq_thread_data *data;
        struct fhmm* fhmm = NULL;
        struct ihmm_sequence* seq = NULL;
        int i;
        int num_threads;
        int thread_id;
        int expected_len;
        float f_score;
        float r_score;
        data = (struct spotseq_thread_data *) threadarg;

        num_threads = data->num_threads;
        thread_id = data->thread_ID;
        fhmm = data->fhmm;

        expected_len = 0;
        for(i = 0; i < data->sb->num_seq;i++){
                expected_len += data->sb->sequences[i]->seq_len;
        }
        expected_len = expected_len / data->sb->num_seq;
        //LOG_MSG("Average sequence length: %d",expected_len);

        for(i =0; i < data->sb->num_seq;i++){
                if( i% num_threads == thread_id){
                        seq = data->sb->sequences[i];
                        RUN(forward(fhmm, data->F_matrix, &f_score, seq->seq, seq->seq_len ));
                        RUN(random_model_score(fhmm->background, &r_score, seq->seq, seq->seq_len,expected_len));
                        //fprintf(stdout,"seq:%d %f %f log-odds: %f  p:%f\n",i, f_score,r_score,f_score - r_score, LOGISTIC_FLT(f_score - r_score));
                        seq->score = (f_score - r_score) / logf(2.0);
                }
        }
        return NULL;
ERROR:
        return NULL;
}

void* do_label_sequences(void* threadarg)
{
        struct spotseq_thread_data *data;
        struct fhmm* fhmm = NULL;
        struct ihmm_sequence* seq = NULL;
        int i;
        int num_threads;
        int thread_id;
        float f_score;
        float b_score;
        data = (struct spotseq_thread_data *) threadarg;

        num_threads = data->num_threads;
        thread_id = data->thread_ID;
        fhmm = data->fhmm;
        //LOG_MSG("Average sequence length: %d",expected_len);

        for(i =0; i < data->sb->num_seq;i++){
                if( i% num_threads == thread_id){
                        seq = data->sb->sequences[i];
                        RUN( forward(fhmm, data->F_matrix, &f_score, seq->seq, seq->seq_len));
                        RUN(backward(fhmm, data->B_matrix, &b_score, seq->seq, seq->seq_len));
                        RUN(posterior_decoding(fhmm,data->F_matrix,data->B_matrix,f_score,seq->seq, seq->seq_len, seq->label));
                }
        }
        return NULL;
ERROR:
        return NULL;
}
