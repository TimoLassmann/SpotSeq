
#include "beam_sample.h"

#include "fast_hmm_param_test_functions.h"

struct beam_thread_data{
        struct fast_hmm_param* ft;
        struct seq_buffer* sb;
        struct fhmm* fhmm;
        float** dyn;
        float** F_matrix;
        float** B_matrix;
        float** t;
        float** e;
        int thread_ID;
        int num_threads;
};

static void* do_score_sequences(void* threadarg);
static void* do_sample_path_and_posterior(void* threadarg);
static void* do_dynamic_programming(void *threadarg);
static void* do_forward_backward(void *threadarg);


int approximatelyEqual(float a, float b, float epsilon);

static struct beam_thread_data** create_beam_thread_data(int* num_threads, int max_len, int K);
static int sum_counts_from_multiple_threads(struct beam_thread_data** td,int* num_threads,int K);
static int resize_beam_thread_data(struct beam_thread_data** td,int* num_threads, int max_len, int K);
static void free_beam_thread_data(struct beam_thread_data** td, int num_threads);

static int transfer_counts(struct ihmm_model* ihmm, float** t, float** e);

static int assign_posterior_probabilities_to_sampled_path(float** F,float** B,float** E, struct ihmm_sequence* ihmm_seq );

//static int set_u(struct seq_buffer* sb, struct ihmm_model* model, float* min_u);
static int set_u(struct seq_buffer* sb, struct ihmm_model* model, struct fast_hmm_param* ft, float* min_u);
static int unset_u(struct seq_buffer* sb);

static int get_max_to_last_state_transition(struct fast_hmm_param*ft,float* max);
//static int check_if_ft_is_indexable(struct fast_hmm_param* ft, int num_states);

static int dynamic_programming(float** matrix,struct fast_hmm_param* ft, struct ihmm_sequence* ihmm_seq);

static int forward_slice(float** matrix,struct fast_hmm_param* ft, struct ihmm_sequence* ihmm_seq, float* score);
static int backward_slice(float** matrix,struct fast_hmm_param* ft, struct ihmm_sequence* ihmm_seq, float* score);
static int collect_slice(struct beam_thread_data* data,struct ihmm_sequence* ihmm_seq, float total);

int run_beam_sampling(struct ihmm_model* model, struct seq_buffer* sb, struct fast_hmm_param* ft,struct thr_pool* pool, int iterations, int num_threads)
{
        int i;
        int iter;
        float min_u;
        float max;
        //float** matrix = NULL;
        struct thr_pool* local_pool = NULL;
        struct beam_thread_data** td = NULL;
        int need_local_pool;

        int no_path;

        ASSERT(model != NULL, "nop model.");
        ASSERT(sb,"no sequence buffer");
        ASSERT(sb->num_seq > 0, "No sequences");
        ASSERT(ft != NULL, "No transition struct");
        ASSERT(iterations > 1, "No iterations");
        ASSERT(num_threads > 0, "No threads");

        init_logsum();
        //for(i = 0;i < 10;i++){
        //        RUN(iHmmHyperSample(model, 20));
        //}

        /* sample transitions / emission */
        RUN(fill_fast_transitions(model,ft));
        //RUN(print_fast_hmm_params(ft));

        //LOG_MSG("LASTSTATE in FT: %d", ft->last_state);
        //exit(0);

        /* Threading setup...  */
        need_local_pool = 0;
        if(pool){
                local_pool = pool;
        }else{
                if((local_pool = thr_pool_create(num_threads,num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");
                need_local_pool =1;
        }

        RUNP(td = create_beam_thread_data(&num_threads,(sb->max_len+2)  , ft->last_state));
        LOG_MSG("Will use %d threads.", num_threads);
        no_path = 0;
        for(iter = 0;iter < iterations;iter++){//}iterations;iter++){
                /* Set U */
                RUN(set_u(sb,model,ft, &min_u));

                /* I only want to add states if the last iteration was successful */
                if(!no_path){
                        RUN(get_max_to_last_state_transition(ft, &max));
                        //fprintf(stdout,"MAX:%f\n", max);
                        while(max > min_u && model->num_states < 300){//}sb->max_len){
                                //fprintf(stdout,"ITER: %d Add state! MAX:%f min_U:%f max_len: %d \n",iter , max, min_u,sb->max_len);
                                RUN(add_state_from_fast_hmm_param(model,ft));
                                RUN(get_max_to_last_state_transition(ft, &max));
                                //fprintf(stdout,"MAX:%f min_U:%f\n", max, min_u);
                                //exit(0);
                        }
                }
                 RUN(make_flat_param_list(ft));
                LOG_MSG("Iteration %d (%d states)  alpha = %f, gamma = %f", iter, model->num_states, model->alpha ,model->gamma);

                RUN(resize_beam_thread_data(td, &num_threads,(sb->max_len+2)  ,model->num_states));

                //dyn prog + labelling
                for(i = 0; i < num_threads;i++){
                        td[i]->ft = ft;
                        td[i]->sb = sb;
                        if(thr_pool_queue(local_pool,do_dynamic_programming,td[i]) == -1){
                                fprintf(stderr,"Adding job to queue failed.");
                        }
                        /* if(thr_pool_queue(local_pool,do_forward_backward,td[i]) == -1){ */
                        /*         fprintf(stderr,"Adding job to queue failed."); */
                        /* } */
                        /* if(thr_pool_queue(local_pool, do_sample_path_and_posterior,td[i]) == -1){ */
                        /*         fprintf(stderr,"Adding job to queue failed."); */
                        /* } */
                }
                thr_pool_wait(local_pool);

                /*struct fhmm* fhmm = NULL;
                RUNP(fhmm = build_finite_hmm_from_infinite_hmm(model));
                for(i = 0; i < num_threads;i++){
                        td[i]->ft = ft;
                        td[i]->sb = sb;
                        td[i]->fhmm = fhmm;
                        if(thr_pool_queue(local_pool, do_score_sequences,td[i]) == -1){
                                fprintf(stderr,"Adding job to queue failed.");
                        }
                }
                thr_pool_wait(local_pool);
                free_fhmm(fhmm);*/
                no_path =0;
                for(i = 0; i < sb->num_seq;i++){
                        if(sb->sequences[i]->u[0] == -1){
                                no_path = no_path + 1;
                                LOG_MSG("weird split must have happened in seq %d",i);
                        }
                        //if(i < 5){
                        //        LOG_MSG("seq:%d had score of %f", i, sb->sequences[i]->score);
                        //}
                }
                /* if more than 1% of sequences don't have a path redo */
                //if((double) no_path / (double) sb->num_seq >= 0.01){
                //        no_path = 1;
                //}else{
                //        no_path = 0;
                //}

                if(no_path){
                        LOG_MSG("weird split must have happened. %d",iter);
                        RUN(fill_fast_transitions(model,ft));
                        iterations++;
                }else{
                        /* I am doing this as a pre-caution. I don't want the inital model
                         * contain states that are not visited.. */
                        RUN(remove_unused_states_labels(model, sb));
                        RUN(fill_counts(model,sb));

                        RUN(iHmmHyperSample(model, 1));
                        model->gamma = 0.1;
                        model->alpha = 0.5;
                        RUN(fill_fast_transitions(model,ft));
                }
                /* print out model - used for plotting  */
                /*if((iter+1) % 10 == 0){
                //LOG_MSG("print %d\n",iter);
                char tmp_buffer[BUFFER_LEN];
                        snprintf(tmp_buffer,BUFFER_LEN,"model_at_%07d.h5",iter+1);
                        RUN(write_model_hdf5(model,tmp_buffer));
                        RUN(add_background_emission(tmp_buffer,ft->background_emission,ft->L));
                        RUN(run_build_fhmm_file(tmp_buffer));
                        }*/
        }

        if(need_local_pool){
                thr_pool_destroy(local_pool);
        }
        free_beam_thread_data(td, num_threads);
        return OK;
ERROR:
        free_beam_thread_data(td, num_threads);
        return FAIL;
}

void* do_forward_backward(void *threadarg)
{
        struct beam_thread_data *data;
        int i,j;
        int num_threads;
        int thread_id;

        float f_score;
        float b_score;
        data = (struct beam_thread_data *) threadarg;

        num_threads = data->num_threads;
        thread_id = data->thread_ID;

        /* clear e and t count tables.  */
        for(i = 0; i < data->ft->last_state;i++){
                for(j =0; j < data->ft->last_state;j++){
                        data->t[i][j] = -INFINITY;
                }

        }
        for(i = 0; i < ALPHABET_PROTEIN;i++){
                for(j =0; j < data->ft->last_state;j++){
                        data->e[i][j] = -INFINITY;
                }
        }
        for(i =0; i < data->sb->num_seq;i++){
                if( i% num_threads == thread_id){
                        // LOG_MSG("Thread %d running sequence %d",thread_id, i);
                        RUN(forward_slice(data->F_matrix,data->ft, data->sb->sequences[i],&f_score));
                        if(f_score == -INFINITY){
                                data->sb->sequences[i]->u[0] = -1;
                        }else{
                                RUN(backward_slice(data->B_matrix,data->ft, data->sb->sequences[i],&b_score));
                                //fprintf(stdout,"%f (f)\n%f (b)\n", f_score,b_score);
                                RUN(collect_slice(data, data->sb->sequences[i], f_score));
                        }
                }
        }
        return NULL;
ERROR:
        return NULL;
}

void* do_score_sequences(void* threadarg)
{
        struct beam_thread_data *data;
        struct fhmm* fhmm = NULL;
        struct ihmm_sequence* seq = NULL;
        int i;
        int num_threads;
        int thread_id;
        int expected_len;
        float f_score;
        float r_score;
        data = (struct beam_thread_data *) threadarg;

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
                        fprintf(stdout,"seq:%d %f %f log-odds: %f  p:%f\n",i, f_score,r_score,f_score - r_score, LOGISTIC_FLT(f_score - r_score));
                }
        }
        return NULL;
ERROR:
        return NULL;
}

void* do_sample_path_and_posterior(void* threadarg)
{
        struct beam_thread_data *data;
        struct ihmm_sequence* seq = NULL;
        int i;
        int num_threads;
        int thread_id;
        float f_score;
        float b_score;
        float r_score;

        data = (struct beam_thread_data *) threadarg;

        num_threads = data->num_threads;
        thread_id = data->thread_ID;

        for(i =0; i < data->sb->num_seq;i++){
                if( i% num_threads == thread_id){
                        seq = data->sb->sequences[i];
                        //               LOG_MSG("Thread %d running sequence %d",thread_id, i);
                        RUN(dynamic_programming(data->dyn,data->ft, seq));

                        if(seq->u[0] != -1){
                                RUN(forward_slice(data->F_matrix, data->ft, seq, &f_score));
                                RUN(backward_slice(data->B_matrix, data->ft, seq, &b_score));
                                RUN(random_model_score(data->ft->background_emission, &r_score, seq->seq, seq->seq_len,seq->seq_len));
                                if(!approximatelyEqual(f_score, b_score, 10e-5)){
                                        fprintf(stdout,"%f %f %d (%0.8f)\n", f_score,b_score, approximatelyEqual(f_score, b_score, 10e-5), 10e-5);
                                }
                                fprintf(stdout,"seq: %d\tp:%f f:%f r:%f diff:%f  %f\t%f \n",i,seq->score, f_score,r_score, seq->score - f_score,f_score-r_score, LOGISTIC_FLT(f_score-r_score));
                                seq->score = f_score;

                                RUN(assign_posterior_probabilities_to_sampled_path(data->F_matrix,data->B_matrix,data->ft->emission, seq));

                        }
                }
        }
        return NULL;
ERROR:
        return NULL;
}

void* do_dynamic_programming(void *threadarg)
{
        struct beam_thread_data *data;
        int i;
        int num_threads;
        int thread_id;

        data = (struct beam_thread_data *) threadarg;

        num_threads = data->num_threads;
        thread_id = data->thread_ID;

        for(i =0; i < data->sb->num_seq;i++){
                if( i% num_threads == thread_id){
                        //               LOG_MSG("Thread %d running sequence %d",thread_id, i);
                        RUN(dynamic_programming(data->dyn,data->ft, data->sb->sequences[i]));

                        /*while(data->sb->sequences[i]->score == -INFINITY){
                                RUN(dynamic_programming(data->dyn,data->ft, data->sb->sequences[i]));
                                }*/
                }
        }
        return NULL;
ERROR:
        return NULL;
}


/* This function assumes (oh no!) that beta has space for an additional
   g * element */
int add_state_from_fast_hmm_param(struct ihmm_model* ihmm,struct fast_hmm_param* ft)
{
        struct fast_t_item** infinity = NULL;
        struct fast_t_item* tmp = NULL;
        float* tmp_prob = NULL;

        float* beta;
        float alpha;
        float gamma;
        rk_state rndstate;

        float sum,be,bg,pe,pg, a,b;
        int i,new_k;//,list_index;
        //intl,r;

        //int pg_hack;            /* I don't want add states that are not reachable. */
        //float* tmp_pg = NULL;


        ASSERT(ihmm != NULL, "No model");
        ASSERT(ft != NULL, "No ft.");
        /* Sorting is only strictly necessary if this is called after another function re-sorted it */
        //qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_from_asc);


        rndstate = ihmm->rndstate;


        //list_index = ft->num_items;

        /* First add empty space to host the newstate -> old state transitions. */
        //if(list_index + ft->last_state + ft->last_state + 1 >= ft->alloc_num_states){
        //        LOG_MSpG("requesting more memory in add state...");
        //RUN(expand_fast_hmm_param_if_necessary(ft, list_index + ft->last_state + ft->last_state + 1));
        //}
        /* Check if model needs to be extended (mainly beta of course) */

        RUN(resize_ihmm_model(ihmm, ihmm->num_states + 1));

        ihmm->num_states = ihmm->num_states + 1;
        RUN(expand_ft_if_necessary(ft, ihmm->num_states));

        MMALLOC(tmp_prob, sizeof(float) *(ihmm->num_states));

        beta = ihmm->beta;
        alpha = ihmm->alpha;
        gamma = ihmm->gamma;

        new_k = ft->last_state;
        infinity = ft->infinity;
        //fprintf(stdout,"LAST: %d\n",new_k);
        /* fill out transition FROM new state  */
        sum = 0.0f;
        for(i = 0;i <= new_k;i++){
                tmp_prob[i] =  rk_gamma(&rndstate, beta[i] * alpha, 1.0);
                if(i == IHMM_START_STATE){
                        tmp_prob[i] = 0.0f;
                }
                sum += tmp_prob[i];
        }
        for(i = 0;i < new_k;i++){
                tmp = NULL;
                MMALLOC(tmp, sizeof(struct fast_t_item));
                tmp->from = new_k;
                tmp->to = i;
                tmp->t = tmp_prob[i] / sum;
                ft->root->tree_insert(ft->root,tmp);
                ft->transition[new_k][i] = tmp->t;
        }
        infinity[new_k]->from = new_k;
        infinity[new_k]->to = new_k;
        infinity[new_k]->t = tmp_prob[new_k] / sum;
        ft->transition[new_k][new_k] = infinity[new_k]->t;

        /*list = ft->list;
          list_index = ft->num_items;
          sum = 0.0;
          for(i = 0;i <= ft->last_state;i++){
          list[list_index]->from = new_k;
          list[list_index]->to = i;
          if(i!= IHMM_START_STATE){
          list[list_index]->t = rk_gamma(&rndstate, beta[i] * alpha, 1.0);
          }else{
          list[list_index]->t = 0.0;
          }
          sum += list[list_index]->t;
          list_index++;
          if(list_index == ft->alloc_items){
          RUN(expand_transition_if_necessary(ft));
          list = ft->list;
          }

          }
          for(i = ft->num_items;i < list_index;i++){
          list[i]->t /= sum;
          }
          ft->num_items = list_index;*/
        //first get beta for new column
        be = beta[new_k];
        bg = rk_beta(&rndstate, 1.0,gamma );

        beta[new_k] = bg*be;
        beta[new_k+1] = (1.0 - bg) *be;

        ihmm->beta = beta;
        //now split prob in last columns...
        a = alpha * beta[new_k];
        b = 0.0;
        for(i = 0; i <= new_k;i++){

                b += beta[i];
        }
        b = alpha * (1.0 - b);

        /*
          MMALLOC(tmp_pg, sizeof(float)* (ft->last_state+1));

          pg_hack = -1;
          while(pg_hack == -1){
          for(i = 0; i < ft->last_state+1;i++){

          if(a < 1e-2 || b < 1e-2){     // % This is an approximation when a or b are really small.
          pg = rk_binomial(&rndstate, 1.0, a / (a+b));
          }else{
          pg = rk_beta(&rndstate, a, b);
          }
          tmp_pg[i] = pg;

          }
          for(i = 0; i < ft->last_state;i++){
          if(i != IHMM_END_STATE){
          if(tmp_pg[i] != 1){
          pg_hack = 1;
          }
          }
          }

          }
          for(i = 0; i < ft->last_state+1;i++){
          fprintf(stdout,"from:%d pg:%f\n",i,tmp_pg[i]);
          }
        */

        // split last column - i.e. play with infinity.

        for(i = 0 ; i <= new_k;i++){
                if(a < 1e-2 || b < 1e-2){     // % This is an approximation when a or b are really small.
                        pg = rk_binomial(&rndstate, 1.0, a / (a+b));
                }else{
                        pg = rk_beta(&rndstate, a, b);
                }
                pe = infinity[i]->t;

                //transition to state just instantiated will go into the RB tree.
                tmp = NULL;
                MMALLOC(tmp, sizeof(struct fast_t_item));
                tmp->from = i;
                tmp->to = new_k;
                tmp->t = pg * pe;
                ft->root->tree_insert(ft->root,tmp);
                ft->transition[i][new_k] = tmp->t;

                //transition into infinity will remain in the infinity array...
                infinity[i]->from = i;
                infinity[i]->to = new_k+1;
                infinity[i]->t = (1.0-pg) * pe;
                ft->transition[i][new_k+1] = infinity[i]->t;



        }

        /*qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_asc);


          l = fast_hmm_param_binarySearch_to_lower_bound(ft,ft->last_state);
          r = fast_hmm_param_binarySearch_to_upper_bound(ft,ft->last_state);

          for(i = l;i < r;i++){
          if(a < 1e-2 || b < 1e-2){     // % This is an approximation when a or b are really small.
          pg = rk_binomial(&rndstate, 1.0, a / (a+b));
          }else{
          pg = rk_beta(&rndstate, a, b);
          }
          pe = list[i]->t;
          //fprintf(stdout,"Filling in %d -> %d : %f to %f   PG:%f\n",list[i]->from,list[i]->to,pe,pg*pe ,pg );
          list[i]->t = pg * pe;

          list[list_index]->from = list[i]->from;
          list[list_index]->to = new_k+1;
          list[list_index]->t = (1.0-pg) * pe;

          //fprintf(stdout,"Filling in %d -> %d : %f to %f\n",list[i]->from,list[i]->to,pe,(1.0-pg) * pe);

          list_index++;
          if(list_index == ft->alloc_items){
          RUN(expand_transition_if_necessary(ft));
          list = ft->list;
          }
          }*/


        /* add emission  */
        sum = 0.0;
        for(i = 0; i < ihmm->L;i++){
                ft->emission[i][new_k] =rk_gamma(&rndstate, EMISSION_H, 1.0);
                sum += ft->emission[i][new_k];
        }
        for(i = 0; i < ihmm->L;i++){
                ft->emission[i][new_k] /= sum;
        }


        //MFREE(tmp_pg);
        //ft->num_items = list_index;
        ft->last_state = new_k+1;
        ihmm->rndstate = rndstate;
        MFREE(tmp_prob);
        return OK;
ERROR:
        //if(tmp_pg){
        //        MFREE(tmp_pg);
        // }
        if(tmp_prob){
                MFREE(tmp_prob);
        }
        return FAIL;
}


struct beam_thread_data** create_beam_thread_data(int* num_threads, int max_len, int K)
{
        struct beam_thread_data** td = NULL;
        int i;
        int local_num_treads;
        size_t mem_needed;
        ASSERT(*num_threads> 0, "no threads");
        local_num_treads = *num_threads;


        mem_needed = sizeof(float) * local_num_treads * max_len * K;

        while(mem_needed  > GB){
                local_num_treads--;
                ASSERT(local_num_treads != 0, "No space! %d asked for but the limit is %d",mem_needed,GB);
                mem_needed = sizeof(float) * local_num_treads * max_len * K;

        }

        MMALLOC(td, sizeof(struct beam_thread_data*) * local_num_treads);
        for(i = 0; i < local_num_treads;i++){
                td[i] = NULL;
                MMALLOC(td[i], sizeof(struct beam_thread_data));
                td[i]->dyn = NULL;
                td[i]->F_matrix = NULL;
                td[i]->B_matrix = NULL;
                td[i]->t = NULL;
                td[i]->e = NULL;
                td[i]->fhmm = NULL;
                //  RUNP(matrix = malloc_2d_float(matrix,sb->max_len+1, ft->last_state, 0.0f));
                RUNP(td[i]->dyn = malloc_2d_float(td[i]->dyn, max_len, K, 0.0f));
                RUNP(td[i]->F_matrix = malloc_2d_float(td[i]->F_matrix, max_len, K, 0.0f));
                RUNP(td[i]->B_matrix = malloc_2d_float(td[i]->B_matrix, max_len, K, 0.0f));

                RUNP(td[i]->t = malloc_2d_float(td[i]->t, K,K, -INFINITY));
                RUNP(td[i]->e = malloc_2d_float(td[i]->e,ALPHABET_PROTEIN,K,-INFINITY));


                td[i]->ft = NULL;
                td[i]->sb = NULL;
                td[i]->thread_ID = i;
                td[i]->num_threads = local_num_treads;


        }

        *num_threads = local_num_treads;

        return td;
ERROR:
        free_beam_thread_data(td, *num_threads);
        return NULL;
}

int transfer_counts(struct ihmm_model* ihmm, float** t, float** e)
{
        float* used_states = NULL;
        float sum;
        int K = ihmm->num_states;
        int new_K;
        int i,j,a,b;

        MMALLOC(used_states, sizeof(float) * K);
        for(i = 0; i < K;i++){
                used_states[i] = 0.0;
        }
        used_states[IHMM_END_STATE] = 100;
        used_states[IHMM_START_STATE] = 100;

        for(i = 0; i <K; i++){
                for(j = 0; j < K; j++){
                        ihmm->transition_counts[i][j] = 0.0;
                }

        }

        for(i = 0; i < ihmm->L; i++){
                for(j = 0; j < K; j++){
                        used_states[j] += scaledprob2prob(e[i][j]);
                        ihmm->emission_counts[i][j] = 0.0;
                }
        }

        new_K = 0;
        sum = 0;
        for(i = 0; i < K;i++){
                fprintf(stdout,"%d : %0.10f  beta: %f \n",i , used_states[i], ihmm->beta[i]);
                if(used_states[i]){
                        ihmm->beta[new_K] = ihmm->beta[i];
                        used_states[i] = new_K;
                        new_K++;
                }else{
                        used_states[i] = -1;
                        sum += ihmm->beta[i];
                }
        }

        ihmm->beta[new_K] = sum;

        ihmm->num_states = new_K+1;

        RUN(resize_ihmm_model(ihmm, new_K+1));



        sum = 0;
        fprintf(stdout,"\n");
        for(i = 0; i < K;i++){
                if(i <= new_K){
                        sum += ihmm->beta[i];
                }
                fprintf(stdout,"%d : %f  beta: %f \n",i , used_states[i],ihmm->beta[i]);
        }
        fprintf(stdout,"SUM:%f \n", sum);
        for(i = 0; i < K; i++){
                if(used_states[i] != -1){
                        a = used_states[i];

                        for(j = 0; j < K; j++){
                                if(used_states[j] != -1){
                                        b = used_states[j];

                                        ihmm->transition_counts[a][b] = scaledprob2prob(t[i][j]);
                                }
                        }
                }

        }
        for(i = 0; i < ihmm->L; i++){
                for(j = 0; j < K; j++){
                        if(used_states[j] != -1){
                                b = used_states[j];
                                ihmm->emission_counts[i][b] = scaledprob2prob(e[i][j]);
                        }
                }
        }

        //ihmm->num_states = new_K;
        MFREE(used_states);
        return OK;
ERROR:
        return FAIL;
}

int sum_counts_from_multiple_threads(struct beam_thread_data** td,int* num_threads,int K)
{
        int i,j,c;
        int local_num_treads;

        local_num_treads = *num_threads;
        for(c = 1; c <  local_num_treads;c++){
                for(i = 0; i < K; i++){
                        for(j = 0; j < K; j++){
                                td[0]->t[i][j] = logsum(td[0]->t[i][j],  td[c]->t[i][j]);
                        }

                }
                for(i = 0; i < ALPHABET_PROTEIN; i++){
                        for(j = 0; j < K; j++){
                                td[0]->e[i][j] = logsum(td[0]->e[i][j],  td[c]->e[i][j]);

                        }

                }


        }

        return OK;
}

int resize_beam_thread_data(struct beam_thread_data** td,int* num_threads, int max_len, int K)
{
        int i;
        int local_num_treads;
        int cur_threads;
        size_t mem_needed;
        ASSERT(*num_threads> 0, "no threads");
        local_num_treads = *num_threads;
        cur_threads =  *num_threads;

        mem_needed = sizeof(float) * local_num_treads * max_len * K;

        while(mem_needed  > GB){
                local_num_treads--;
                ASSERT(local_num_treads != 0, "No space! %d asked for but the limit is %d",mem_needed,GB);
                mem_needed = sizeof(float) * local_num_treads * max_len * K;
        }


        for(i = local_num_treads; i < cur_threads;i++){

                free_2d((void**) td[i]->dyn);
                MFREE(td[i]);
        }

        //LOG_MSG("mallocing auxiliary datastructures to %d %d", max_len,K);
        for(i = 0; i < local_num_treads;i++){
                RUNP(td[i]->dyn = malloc_2d_float(td[i]->dyn, max_len, K, 0.0f));

                RUNP(td[i]->F_matrix = malloc_2d_float(td[i]->F_matrix, max_len, K, 0.0f));
                RUNP(td[i]->B_matrix = malloc_2d_float(td[i]->B_matrix, max_len, K, 0.0f));

                RUNP(td[i]->t = malloc_2d_float(td[i]->t, K,K, -INFINITY));
                RUNP(td[i]->e = malloc_2d_float(td[i]->e,ALPHABET_PROTEIN,K,-INFINITY));

                td[i]->num_threads = local_num_treads;
        }
        *num_threads = local_num_treads;
        return OK;
ERROR:
        return FAIL;
}


void free_beam_thread_data(struct beam_thread_data** td, int num_threads)
{
        int i;
        if(td){
                for(i = 0; i < num_threads;i++){
                        free_2d((void**) td[i]->dyn);
                        free_2d((void**) td[i]->F_matrix);
                        free_2d((void**) td[i]->B_matrix);
                        free_2d((void**) td[i]->t);
                        free_2d((void**) td[i]->e);

                        MFREE(td[i]);
                }

                MFREE(td);
        }
}





int approximatelyEqual(float a, float b, float epsilon)
{
        return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


int assign_posterior_probabilities_to_sampled_path(float** F,float** B,float** E, struct ihmm_sequence* ihmm_seq )
{
        /* got through path assign emission probs, store in u... */
        float* emission = NULL;
        int* label = NULL;
        float* u = NULL;
        uint8_t* seq = NULL;
        float total = 0.0;
        int i,l,len;

        u = ihmm_seq->u;
        len = ihmm_seq->seq_len;
        seq = ihmm_seq->seq;
        label = ihmm_seq->label;
        total = ihmm_seq->score;

        for (i = 0; i < len; i++) {
                l = label[i];
                emission = E[seq[i]];

                u[i] = F[i][l] + (B[i][l] - prob2scaledprob(emission[l])) - total;
                //fprintf(stdout,"%d letter:%d label:%d F:%f B:%f total=%f  result:%f \n",i,seq[i],l,F[i][l],B[i][l],total,scaledprob2prob(u[i]));
        }
        //exit(0);

        return OK;
}




int collect_slice(struct beam_thread_data* data,struct ihmm_sequence* ihmm_seq, float total)
{
        float** e = data->e;
        float** t = data->t;
        float** F = data->F_matrix;
        float** B = data->B_matrix;
        float* emission = NULL;

        struct fast_hmm_param* ft = data->ft;
        struct fast_t_item** list = NULL;
        float* u = NULL;
        uint8_t* seq = NULL;
        int i,j,a,b,l,len,boundary;

        u = ihmm_seq->u;
        len = ihmm_seq->seq_len;
        seq = ihmm_seq->seq;

        list = ft->list;

        l = ft->last_state;

        boundary = fast_hmm_param_binarySearch_t(ft, u[0]);
        //fill first row.
        for(j = 0; j < boundary;j++){
                if(list[j]->from == IHMM_START_STATE){
                        t[IHMM_START_STATE][list[j]->to] = logsum(t[IHMM_START_STATE][list[j]->to], prob2scaledprob(list[j]->t) + B[0][list[j]->to] - total);

                }
        }
        emission = ft->emission[seq[0]];
        //fprintf(stdout,"L:%d\n",seq[0]);
        for(i = 0; i < l;i++){
                e[seq[0]][i] = logsum(e[seq[0]][i], (F[0][i] + (B[0][i] - prob2scaledprob(emission[i]) ))  - total);
        }

        for(i = 1; i < len;i++){
                boundary = fast_hmm_param_binarySearch_t(ft, u[i]);

                for(j = 0; j < boundary;j++){
                        a = list[j]->from;
                        b = list[j]->to;
                        t[a][b]  = logsum( t[a][b], F[i-1][a] + prob2scaledprob(list[j]->t) + B[i][b] - total);
                }
                emission = ft->emission[seq[i]];
                //fprintf(stdout,"L:%d\n",seq[i]);
                for(j = 0; j < l;j++){
                        e[seq[i]][j] = logsum(e[seq[i]][j],  (F[i][j] + (B[i][j] - prob2scaledprob(emission[j] )))  - total);
                }
        }


        /* First let's check if there is a path! i.e. end is reachable.  */



        boundary = fast_hmm_param_binarySearch_t(ft, u[len]);
        for(j = 0; j < boundary;j++){
                a = list[j]->from;
                b = list[j]->to;
                if(b == IHMM_END_STATE){
                        t[a][b]  = logsum( t[a][b], F[len-1][a] + prob2scaledprob(list[j]->t)  - total);

                }
        }
        return OK;
}


int forward_slice(float** matrix,struct fast_hmm_param* ft, struct ihmm_sequence* ihmm_seq, float* score)
{
        struct fast_t_item** list = NULL;
        float* u = NULL;
        uint8_t* seq = NULL;
        float* emission = NULL;
        float* cur = NULL;
        float* last = NULL;



        int i,j,len,boundary;
        int l;
        int a,b;

        ASSERT(ft!= NULL, "no parameters");
        ASSERT(matrix != NULL,"No dyn matrix");

        u = ihmm_seq->u;
        len = ihmm_seq->seq_len;
        seq = ihmm_seq->seq;

        list = ft->list;

        cur = matrix[0];
        l = ft->last_state;

        for(i = 0; i < l;i++){
                cur[i]  = -INFINITY;
        }

        boundary = fast_hmm_param_binarySearch_t(ft, u[0]);
        //fill first row.
        for(j = 0; j < boundary;j++){
                if(list[j]->from == IHMM_START_STATE){
                        cur[list[j]->to] = logsum(cur[list[j]->to], prob2scaledprob( list[j]->t));
                }
        }
        emission = ft->emission[seq[0]];

        for(i = 0; i < l;i++){
                cur[i] += prob2scaledprob(emission[i]);
        }

        for(i = 1; i < len;i++){
                last = cur;
                cur = matrix[i];
                for(j = 0; j < l;j++){
                        cur[j] = -INFINITY;
                }

                boundary = fast_hmm_param_binarySearch_t(ft, u[i]);

                for(j = 0; j < boundary;j++){
                        a = list[j]->from;
                        b = list[j]->to;
                        //fprintf(stdout,"lastA: %f  prob:  %f from: %d to: %d\n",last[a] , prob2scaledprob(list[j]->t) ,a,b);
                        cur[b] = logsum(cur[b],last[a] + prob2scaledprob(list[j]->t));

                }
                emission = ft->emission[seq[i]];
                for(j = 0; j < l;j++){
                        cur[j] += prob2scaledprob(emission[j]);
                }
        }


        /* First let's check if there is a path! i.e. end is reachable.  */

        *score = prob2scaledprob(0.0f);

        boundary = fast_hmm_param_binarySearch_t(ft, u[len]);
        for(j = 0; j < boundary;j++){
                a = list[j]->from;
                b = list[j]->to;
                if(b == IHMM_END_STATE){
                        *score = logsum(*score, matrix[len-1][a] + prob2scaledprob(list[j]->t));
                }
        }
        //fprintf(stdout,"SCORE: %f\t(%f)\n",*score,u[len]);


        return OK;
ERROR:
        return FAIL;
}

int backward_slice(float** matrix,struct fast_hmm_param* ft, struct ihmm_sequence* ihmm_seq, float* score)
{
        struct fast_t_item** list = NULL;
        float* u = NULL;
        uint8_t* seq = NULL;
        float* emission = NULL;
        float* cur = NULL;
        float* next = NULL;

        int i,j,len,boundary;
        int l;
        int a,b;


        ASSERT(ft!= NULL, "no parameters");
        ASSERT(matrix != NULL,"No dyn matrix");

        u = ihmm_seq->u;
        len = ihmm_seq->seq_len;
        seq = ihmm_seq->seq;

        list = ft->list;

        cur = matrix[len-1];
        l = ft->last_state;

        for(i = 0; i < l;i++){
                cur[i]  = -INFINITY;
        }

        boundary = fast_hmm_param_binarySearch_t(ft, u[len]);
        //fill first row.
        for(j = 0; j < boundary;j++){
                if(list[j]->to  == IHMM_END_STATE){
                        cur[list[j]->from] = logsum(cur[list[j]->from], prob2scaledprob( list[j]->t));
                }
        }
        emission = ft->emission[seq[len-1]];
        for(i = 0; i < l;i++){
                cur[i] += prob2scaledprob(emission[i]);
        }
        for(i = len -2; i >= 0;i--){
                next = cur;
                cur = matrix[i];
                for(j = 0; j < l;j++){
                        cur[j] = -INFINITY;
                }

                boundary = fast_hmm_param_binarySearch_t(ft, u[i+1]);

                for(j = 0; j < boundary;j++){
                        a = list[j]->from;
                        b = list[j]->to;
                        cur[a] = logsum(cur[a], next[b] + prob2scaledprob(list[j]->t));
                }
                emission = ft->emission[seq[i]];


                for(j = 0; j < l;j++){
                        cur[j] += prob2scaledprob(emission[j]);
                }
        }


        /* First let's check if there is a path! i.e. end is reachable.  */

        *score = prob2scaledprob(0.0f);

        boundary = fast_hmm_param_binarySearch_t(ft, u[0]);
        for(j = 0; j < boundary;j++){
                a = list[j]->from;
                b = list[j]->to;
                if(a == IHMM_START_STATE){
                        *score = logsum(*score, matrix[0][b]+ prob2scaledprob(list[j]->t));
                }
        }
//        fprintf(stdout,"SCORE: %f\t(%f)\n",*score,u[0]);

        return OK;
ERROR:
        return FAIL;
}



int dynamic_programming(float** matrix,struct fast_hmm_param* ft, struct ihmm_sequence* ihmm_seq)
{
        int i,j,len,boundary;
        float* u = NULL;
        uint8_t* seq = NULL;
        int* label = NULL;
        int a,b;
        float score;
        float sum;
        float* emission;
        float* tmp_row;
        float r;
        int l;
        struct fast_t_item** list = NULL;
        ASSERT(ft!= NULL, "no parameters");
        ASSERT(matrix != NULL,"No dyn matrix");

        u = ihmm_seq->u;
        len = ihmm_seq->seq_len;
        seq = ihmm_seq->seq;
        label = ihmm_seq->label;

        list = ft->list;

        tmp_row = matrix[len];

        l = ft->last_state;

        boundary = fast_hmm_param_binarySearch_t(ft, u[0]);
        for(i = 0; i < l;i++){
                matrix[0][i] = 0.0f;
        }
        //fill first row.
        for(j = 0; j < boundary;j++){
                if(list[j]->from == IHMM_START_STATE){
                        matrix[0][list[j]->to] += list[j]->t;
                }
        }
        sum = 0;
        emission = ft->emission[seq[0]];
        for(i = 0; i < l;i++){

                matrix[0][i] *=  emission[i];

                sum += matrix[0][i];
        }
        for(i = 0; i < l;i++){
                matrix[0][i] /= sum;
        }
        /*sum = 0.0;

          fprintf(stdout,"%d %f s:%d",i, u[0],seq[0]);
          for(i = 0; i < ft->last_state;i++){
          fprintf(stdout," %f",matrix[0][i]);
          sum += matrix[0][i];
          }
          fprintf(stdout," sum: %f\n",sum);*/
        for(i = 1; i < len;i++){
                emission = ft->emission[seq[i]];
                for(j = 0; j < ft->last_state;j++){
                        matrix[i][j] = 0.0f;
                }

                boundary = fast_hmm_param_binarySearch_t(ft, u[i]);

                for(j = 0; j < boundary;j++){
                        a = list[j]->from;
                        b = list[j]->to;
                        matrix[i][b] += matrix[i-1][a];

                }
                sum = 0.0;

                for(j = 0; j < l;j++){
                        matrix[i][j] *=  emission[j];
                        sum += matrix[i][j];

                }
                for(j = 0; j < l;j++){
                        matrix[i][j] /= sum;
                }
                /*sum = 0.0;
                  fprintf(stdout,"%d %f s:%d",i, u[i],seq[i]);
                  for(j = 0; j < ft->last_state;j++){
                  fprintf(stdout," %f",matrix[i][j]);
                  sum += matrix[i][j];
                  }
                  fprintf(stdout," sum: %f\n",sum);*/
        }
        /* Backtracking...  */
        /* Pick last label based on probabilities in last row. Then look for
         * transitions to that label with a prob > min_u and select ancestor
         * based on probs in previous row */
        l = IHMM_END_STATE;
        /* First let's check if there is a path! i.e. end is reachable.  */

        sum = 0.0f;

        score = prob2scaledprob(1.0f);

        boundary = fast_hmm_param_binarySearch_t(ft, u[len]);
        for(j = 0; j < boundary;j++){
                a = list[j]->from;
                b = list[j]->to;
                if(b == l){
                        sum += matrix[len-1][a];
                }
        }

        if(sum != 0.0 && !isnan(sum)){
                l = IHMM_END_STATE;
                for(i = len-1; i >= 0; i--){
                        //fprintf(stdout,"pick: %d %d\n",i,l);
                        for(j = 0; j < ft->last_state;j++){
                                tmp_row[j] = -1.0;
                        }
                        sum = 0.0f;
                        boundary = fast_hmm_param_binarySearch_t(ft, u[i+1]);
                        for(j = 0; j < boundary;j++){
                                a = list[j]->from;
                                b = list[j]->to;
                                if(b == l){
                                        tmp_row[a] = matrix[i][a];
                                        sum += matrix[i][a];
                                }
                        }
                        r =  random_float_zero_to_x(sum);
                        for(j = 0; j < boundary;j++){

                                a = list[j]->from;
                                b = list[j]->to;
                                if(list[j]->to == l){

                                        r -= tmp_row[a];
                                        if(r <= 0.0f){
                                                l = a;
                                                score = score + prob2scaledprob(list[j]->t);
                                                break;
                                        }
                                }
                        }
                        score = score + prob2scaledprob( ft->emission[seq[i]][l]);
                        label[i] = l;
                }
                /* go from start to first state used in path... */
                score = score + prob2scaledprob(  ft->transition[IHMM_START_STATE][l]);
                ihmm_seq->score = score;

        }else{
                //u[0] = -1.0f;
                ihmm_seq->score = -INFINITY;
                //LOG_MSG("No PATH!: %f",sum);
        }
        //ihmm_seq->score = score;
        return OK;
ERROR:
        return FAIL;
}


int set_u(struct seq_buffer* sb, struct ihmm_model* model, struct fast_hmm_param* ft, float* min_u)
{
        int i,j;
        float* u = 0;
        int* label =0;
        float x;
        int len;

        float local_min_u = 1.0;
        ASSERT(sb != NULL, "No sequences.");
        ASSERT(model != NULL, "No model.");
        //qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_from_asc);
        //last_state = ft->last_state;

        for(i = 0; i < sb->num_seq;i++){
                label = sb->sequences[i]->label;
                u = sb->sequences[i]->u;
                len = sb->sequences[i]->seq_len;

                x = ft->transition[IHMM_START_STATE][label[0]];
                //c = IHMM_START_STATE * last_state + label[0];
                //c = a* (num_states-1) + b;
                u[0] =  rk_double(&model->rndstate) *x;
                //ASSERT(ft->list[c]->t != 0.0f,"BAD %d -> %d %f",ft->list[c]->from,ft->list[c]->to,ft->list[c]->t);

                local_min_u = MACRO_MIN(local_min_u, u[0]);
                for (j = 1; j < len;j++){
                        //c = label[j-1] * last_state + label[j];
                        x = ft->transition[label[j-1]][label[j]];
                        u[j] =  rk_double(&model->rndstate) * x;//rk_double(&model->rndstate) *
                        //if(!i && j < 5){
                        //       fprintf(stdout,"%d->%d %f\n",label[j-1],label[j],ft->list[c]->t );
                        //}
                        //fprintf(stdout,"%d %d  ;; %d %d\n",label[j-1],label[j],ft->list[c]->from ,ft->list[c]->to);
                        local_min_u = MACRO_MIN(local_min_u, u[j]);
                        //ASSERT(ft->list[c]->t != 0.0f,"BAD %d -> %d %f",ft->list[c]->from,ft->list[c]->to,ft->list[c]->t);


                }
                x = ft->transition[label[len-1]][IHMM_END_STATE];

                u[len] =  rk_double(&model->rndstate) * x;//(ft->list[c]->t);
                //ASSERT(ft->list[c]->t != 0.0f,"BAD %d -> %d %f",ft->list[c]->from,ft->list[c]->to,ft->list[c]->t);
                //fprintf(stdout,"%d %d -> %d: %f  \n",label[len-1],ft->list[c]->from ,ft->list[c]->to, ft->list[c]->t );
                local_min_u = MACRO_MIN(local_min_u, u[len]);

        }

        *min_u = local_min_u;
        return OK;
ERROR:
        return FAIL;
}

int unset_u(struct seq_buffer* sb)
{
        int i,j;
        float* u = 0;
        int len;

        ASSERT(sb != NULL, "No sequences.");
        //qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_from_asc);
        //last_state = ft->last_state;

        for(i = 0; i < sb->num_seq;i++){
                u = sb->sequences[i]->u;
                len = sb->sequences[i]->seq_len;
                for (j = 0; j < len;j++){

                        u[j] = 0.0f;
                }

        }
        return OK;
ERROR:
        return FAIL;
}

int get_max_to_last_state_transition(struct fast_hmm_param*ft,float* max)
{

        int i;
        float local_max;

        ASSERT(ft != NULL, "No fast hmm parameters.");

        local_max = -1.0f;


        for(i = 0; i< ft->last_state;i++){
                if(ft->infinity[i]->t > local_max){
                        local_max = ft->infinity[i]->t;
                }
                //fprintf(stdout,"%d->%d %f\n", ft->infinity[i]->from, ft->infinity[i]->to, ft->infinity[i]->t);
        }
        *max = local_max;
        return OK;
ERROR:
        return FAIL;
}



#ifdef ITESTBEAM
/* These are test functions. */

static int add_state_integration_test(void);
static int shrink_grow_integration_test(void);
static int full_run_test(void);
static int full_run_test_protein(void);


int main(const int argc,const char * argv[])
{
        RUN(print_program_header((char * const*)argv,"Integration Test"));
        LOG_MSG("Start add state test");
        RUN(add_state_integration_test());

        LOG_MSG("DONE add state test");

        LOG_MSG("Start shrink / grow test");
        RUN(shrink_grow_integration_test());

        LOG_MSG("DONE shrink / grow test");

        LOG_MSG("START run full test");
        RUN(full_run_test());
        LOG_MSG("DONE run full test");

        LOG_MSG("START run full test (protein)");
        RUN(full_run_test_protein());
        LOG_MSG("DONE run full test");
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int full_run_test(void)
{

        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;
        struct seq_buffer* sb = NULL;
        char *tmp_seq[8] = {
                "ACAGGCTAAAGGAGGGGGCAGTCCCCA",
                "AGGCTAAAGGAGGGGGCAGTCCCCACC",
                "AGGCTAAAGGAGGGGGCAGTCCCCACC",
                "AGTCCCCACCATATTTGAGTCTTTCTC",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGGAGTGGATATCACAGGCTAAAGGAGGGG"};

        int numseq = 8;
        int initial_states = 8;
        int i;


        /* First initialize beam param struct.  */

        RUNP(sb = create_ihmm_sequences_mem(tmp_seq ,numseq));

        RUNP(model = alloc_ihmm_model(initial_states, sb->L));
        rk_randomseed(&model->rndstate);

        model->alpha0_a = 6.0f;
        model->alpha0_b = 15.0f;
        model->gamma_a = 16.0f;
        model->gamma_b = 4.0f;
        model->alpha = IHMM_PARAM_PLACEHOLDER;
        model->gamma = IHMM_PARAM_PLACEHOLDER;

        RUN(inititalize_model(model, sb,initial_states));// initial_states) );
        for(i = 0;i < 10;i++){
                RUN(iHmmHyperSample(model, 20));
        }

        RUNP(ft = alloc_fast_hmm_param(initial_states,sb->L));
        RUN(fill_background_emission(ft, sb));

        RUN(run_beam_sampling( model, sb, ft,NULL, 10, 2));


        //sb, num thread, guess for aplha and gamma.. iterations.


        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        free_ihmm_sequences(sb);
        return OK;
ERROR:
        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        free_ihmm_sequences(sb);
        return FAIL;
}

int full_run_test_protein(void)
{

        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;
        struct seq_buffer* sb = NULL;
        char *tmp_seq[18] = {
                "RRRAHTQAEQKRRDAIKRGYDDLQTIVPTCQQQDFSIGSQKLSKAIVLQKTIDYIQFLH",
                "RREAHTQAEQKRRDAIKKGYDSLQELVPRCQPNDSSGYKLSKALILQKSIEYIGYL",
                "RRITHISAEQKRRFNIKLGFDTLHGLVSTLSAQPSLKVSKATTLQKTAEYILMLQ",
                "RRAGHIHAEQKRRYNIKNGFDTLHALIPQLQQNPNAKLSKAAMLQKGADHIKQLR",
                "KRILHLHAEQNRRSALKDGFDQLMDIIPDLYSGGVKPTNAVVLAKSADHIRRLQ",
                "KKATHLRCERQRREAINSGYSDLKDLIPQTTTSLGCKTTNAAILFRACDFMSQLK",
                "LRTSHKLAERKRRKEIKELFDDLKDALPLDKSTKSSKWGLLTRAIQYIEQLK",
                "YRRTHTANERRRRGEMRDLFEKLKITLGLLHSSKVSKSLILTRAFSEIQGLT",
                "TRKSVSERKRRDEINELLENLKTIVQNPSDSNEKISHETILFRVFERVSGVD",
                "GHRSETEKQRRDDTNDLLNEFKKIVQKSESEKLSKEEVLFRIVKLLSGIQ",
                "KRAHHNALERKRRDHIKDSFHSLRDSVPSLQGEKASRAQILDKATEYIQYMR",
                "RRAHHNELERRRRDHIKDHFTILKDAIPLLDGEKSSRALILKRAVEFIHVMQ",
                "KRAHHNALERRRRDHIKESFTNLREAVPTLKGEKASRAQILKKTTECIQTMR",
                "GRHVHNELEKRRRAQLKRCLEQLRQQMPLGVDHTRYTTLSLLRGARMHIQKLE",
                "NRSSHNELEKHRRAKLRLYLEQLKQLVPLGPDSTRHTTLSLLKRAKVHIKKLE",
                "SRSTHNEMEKNRRAHLRLCLEKLKGLVPLGPESSRHTTLSLLTKAKLHIKKLE",
                "NRTSHNELEKNRRAHLRNCLDGLKAIVPLNQDATRHTTLGLLTQARALIENLK",
                "NRSTHNELEKNRRAHLRLCLERLKVLIPLGPDCTRHTTLGLLNKAKAHIKKLE"};

        int numseq = 18;
        int initial_states = 8;
        int i;

        RUNP(sb = create_ihmm_sequences_mem(tmp_seq ,numseq));

        RUNP(model = alloc_ihmm_model(initial_states, sb->L));
        rk_randomseed(&model->rndstate);

        model->alpha0_a = 6.0f;
        model->alpha0_b = 15.0f;
        model->gamma_a = 16.0f;
        model->gamma_b = 4.0f;
        model->alpha = IHMM_PARAM_PLACEHOLDER;
        model->gamma = IHMM_PARAM_PLACEHOLDER;

        RUN(inititalize_model(model, sb,initial_states));// initial_states) );
        for(i = 0;i < 10;i++){
                RUN(iHmmHyperSample(model, 20));
        }

        RUNP(ft = alloc_fast_hmm_param(initial_states,sb->L));
        RUN(fill_background_emission(ft, sb));

        RUN(run_beam_sampling( model, sb, ft,NULL,  10,2));
        //sb, num thread, guess for aplha and gamma.. iterations.


        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        free_ihmm_sequences(sb);

        return OK;
ERROR:
        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        free_ihmm_sequences(sb);
        return FAIL;
}


static int shrink_grow_integration_test(void)
{
        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;
        struct seq_buffer* iseq = NULL;

        char *tmp_seq[4] = {
                "ACGT",
                "ACGT",
                "ACGT",
                "ACGT"};

        int initial_states = 4;
        int i;

        RUNP(ft = alloc_fast_hmm_param(initial_states,initial_states));

        RUNP(model = alloc_ihmm_model(initial_states, 4));
        RUNP(iseq = create_ihmm_sequences_mem(tmp_seq ,4));
        RUN(random_label_ihmm_sequences(iseq, initial_states * 10,0.3));
        RUN(fill_counts(model,iseq));
        /* I am doing this as a pre-caution. I don't want the inital model
         * contain states that are not visited.. */
        RUN(remove_unused_states_labels(model, iseq));;
        RUN(fill_counts(model,iseq));
        RUN(print_counts(model));

        model->alpha0_a = 4.0f;
        model->alpha0_b = 2.0f;
        model->gamma_a = 3.0f;
        model->gamma_b = 6.0f;
        model->alpha = IHMM_PARAM_PLACEHOLDER;
        model->gamma = IHMM_PARAM_PLACEHOLDER;
        RUN(iHmmHyperSample(model, 10));
        RUN(print_model_parameters(model));

        RUN(random_label_ihmm_sequences(iseq, 2,0.3));
        RUN(fill_counts(model,iseq));
        RUN(print_counts(model));

        RUN(iHmmHyperSample(model, 10));

        RUN(print_model_parameters(model));
        LOG_MSG("Fill transitions 1");
        RUN(fill_fast_transitions(model,ft));
        RUN(print_fast_hmm_params(ft));
        LOG_MSG("Fill transitions 2");
        RUN(fill_fast_transitions(model,ft));
        RUN(print_fast_hmm_params(ft));
        LOG_MSG("Fill transitions 3");

        RUN(fill_fast_transitions(model,ft));
        RUN(print_fast_hmm_params(ft));
        LOG_MSG("Add 10 states");
        for(i = 0;i < 10;i++){
                RUN(add_state_from_fast_hmm_param(model,ft));
        }
        RUN(print_model_parameters(model));
        RUN(print_fast_hmm_params(ft));

        LOG_MSG("Fill transitions 4");
        RUN(random_label_ihmm_sequences(iseq, 2,0.3));
        RUN(fill_counts(model,iseq));
        RUN(print_counts(model));
        RUN(fill_fast_transitions(model,ft));
        RUN(print_fast_hmm_params(ft));


        RUN(make_flat_param_list(ft));
        //for(i = 0; i< 10;i++){
        //       fprintf(stdout,"%d->%d: %f\n", ft->list[i]->from,ft->list[i]->to,ft->list[i]->t);
        //

        int test = fast_hmm_param_binarySearch_t(ft, ft->list[3]->t);
        fprintf(stdout,"Selected for > %f\n",ft->list[3]->t);
        for(i = 0; i< test;i++){
                fprintf(stdout,"%d->%d: %f\n", ft->list[i]->from,ft->list[i]->to,ft->list[i]->t);
        }

        RUN(random_label_ihmm_sequences(iseq, 6,0.3));
        RUN(fill_counts(model,iseq));
        RUN(print_counts(model));
        RUN(iHmmHyperSample(model, 10));

        RUN(print_model_parameters(model));


        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        free_ihmm_sequences(iseq);

        return OK;
ERROR:
        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        free_ihmm_sequences(iseq);
        return FAIL;
}

static int add_state_integration_test(void)
{
        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;
        struct seq_buffer* iseq = NULL;

        char *tmp_seq[4] = {
                "ACGT",
                "ACGT",
                "ACGT",
                "ACGT"};


        int initial_states = 4;



        /* Let's start allocating all structures  */
        RUNP(ft = alloc_fast_hmm_param(initial_states,initial_states));

        RUNP(model = alloc_ihmm_model(initial_states, 4));
        RUNP(iseq = create_ihmm_sequences_mem(tmp_seq ,4));
        RUN(random_label_ihmm_sequences(iseq, initial_states,0.3));


        RUN(fill_counts(model,iseq));
        RUN(print_counts(model));

        /* I am doing this as a pre-caution. I don't want the inital model
         * contain states that are not visited.. */
        RUN(remove_unused_states_labels(model, iseq));
        RUN(fill_counts(model,iseq));
        RUN(print_counts(model));
        /* Now there are counts but no model parameters. */
        model->alpha0_a = 4.0f;
        model->alpha0_b = 2.0f;
        model->gamma_a = 3.0f;
        model->gamma_b = 6.0f;
        model->alpha = IHMM_PARAM_PLACEHOLDER;
        model->gamma = IHMM_PARAM_PLACEHOLDER;
        RUN(iHmmHyperSample(model, 10));

        /* Now I should have everything ready to go.  */
        RUN(print_model_parameters(model));
        /* Just to verify if everything works..  */
        RUN(remove_unused_states_labels(model, iseq));
        RUN(fill_counts(model,iseq));
        RUN(print_counts(model));
        RUN(print_model_parameters(model));

        LOG_MSG("Fill transitions from counts.");
        RUN(fill_fast_transitions(model,ft));
        RUN(print_fast_hmm_params(ft));
        RUN(print_model_parameters(model));
        LOG_MSG("Add a state.");
        RUN(add_state_from_fast_hmm_param(model,ft));
        RUN(print_model_parameters(model));
        RUN(print_fast_hmm_params(ft));

        LOG_MSG("Add a state.");
        RUN(add_state_from_fast_hmm_param(model,ft));
        RUN(print_model_parameters(model));
        RUN(print_fast_hmm_params(ft));

        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        free_ihmm_sequences(iseq);


        return OK;
ERROR:
        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        free_ihmm_sequences(iseq);
        return FAIL;
}

#endif

