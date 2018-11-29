
#include "beam_sample.h"

#include "fast_hmm_param_test_functions.h"

static void* do_sample_path_and_posterior(void* threadarg);
static void* do_dynamic_programming(void *threadarg);
static void* do_forward_backward(void *threadarg);

int approximatelyEqual(float a, float b, float epsilon);

static int sum_counts_from_multiple_threads(struct spotseq_thread_data** td,int* num_threads,int K);

static int transfer_counts(struct ihmm_model* ihmm, float** t, float** e);

static int assign_posterior_probabilities_to_sampled_path(float** F,float** B,float** E, struct ihmm_sequence* ihmm_seq );

//static int set_u(struct seq_buffer* sb, struct ihmm_model* model, float* min_u);
static int set_u_multi(struct model_bag* model_bag, struct fast_param_bag*  ft_bag, struct seq_buffer* sb);
static int set_u(struct seq_buffer* sb, struct ihmm_model* model, struct fast_hmm_param* ft, float* min_u,int model_index);
static int reset_u_if_no_path(struct fast_hmm_param* ft, float* u,int * label, int len, rk_state* rndstate);
static int unset_u(struct seq_buffer* sb);


static int detect_valid_path(struct seq_buffer* sb,int num_models, int* no_path);
static int reset_valid_path(struct seq_buffer* sb,int num_models);

static int expand_ihmms(struct model_bag* model_bag, struct fast_param_bag* ft_bag);
static int add_state_from_fast_hmm_param(struct ihmm_model* ihmm,struct fast_hmm_param* ft);

static int get_max_to_last_state_transition(struct fast_hmm_param*ft,float* max);
//static int check_if_ft_is_indexable(struct fast_hmm_param* ft, int num_states);

static int dynamic_programming(struct spotseq_thread_data* data, int target);
static int dynamic_programming_clean(struct fast_hmm_param* ft,  float** matrix,uint8_t* seq,int* label,float* u,int len,uint8_t* has_path ,rk_state* random);
static int forward_slice(float** matrix,struct fast_hmm_param* ft, struct ihmm_sequence* ihmm_seq, float* score);
static int backward_slice(float** matrix,struct fast_hmm_param* ft, struct ihmm_sequence* ihmm_seq, float* score);
static int collect_slice(struct spotseq_thread_data* data,struct ihmm_sequence* ihmm_seq, float total);

int run_beam_sampling(struct model_bag* model_bag, struct fast_param_bag* ft_bag, struct seq_buffer* sb,struct spotseq_thread_data** td, struct thr_pool* pool, int iterations, int num_threads)
{
        int i;
        int iter;
        //float min_u;
        //float max;
        //float** matrix = NULL;
        //struct thr_pool* local_pool = NULL;
        //struct spotseq_thread_data** td = NULL;
        //int need_local_pool;

        int no_path;

        //struct ihmm_model* model = NULL;
        struct fast_hmm_param* ft = NULL;

        ASSERT(model_bag != NULL, "no model.");
        ASSERT(sb,"no sequence buffer");
        ASSERT(sb->num_seq > 0, "No sequences");
        ASSERT(ft_bag != NULL, "No transition struct");
        ASSERT(iterations > 1, "No iterations");
        ASSERT(num_threads > 0, "No threads");

        init_logsum();

        /* Threading setup...  */
        /*need_local_pool = 0;
        if(pool){
                local_pool = pool;
        }else{
                if((local_pool = thr_pool_create(num_threads,num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");
                need_local_pool =1;
        }

        RUNP(td = create_spotseq_thread_data(&num_threads,(sb->max_len+2)  ,ft_bag->max_last_state, &model_bag->rndstate));
        LOG_MSG("Will use %d threads.", num_threads);*/


        //sb->org_num_seq = sb->num_seq;
        no_path = 0;                            /* Assume that we don't have a path in the first iteration */
        for(iter = 0;iter < iterations;iter++){//}iterations;iter++){
                /* shuffle and sub-sample sequences (or not...) */
                //RUN(shuffle_sequences_in_buffer(sb));
                /* sample transitions / emission */
                ft_bag->max_last_state = -1;
                model_bag->max_num_states = -1;

                //LOG_MSG("Check labelling at start..(%d)", iter);
                //RUN(check_labels(sb,model_bag->num_models ));
                //LOG_MSG("Done");
                if(!no_path){
                        for(i = 0; i < model_bag->num_models;i++){
                                RUN(remove_unused_states_labels(model_bag->models[i], sb,i ));
                                RUN(fill_counts(model_bag->models[i], sb,i));
                                RUN(iHmmHyperSample(model_bag->models[i], 20));
                                model_bag->max_num_states  = MACRO_MAX(model_bag->max_num_states ,model_bag->models[i]->num_states);
                        }
                }

                no_path = 1;
                while(no_path){
                        no_path = 0;
                        ft_bag->max_last_state = -1;
                        for(i = 0; i < model_bag->num_models;i++){
                                RUN(fill_fast_transitions(model_bag->models[i], ft_bag->fast_params[i]));
                                ft_bag->max_last_state = MACRO_MAX(ft_bag->max_last_state,ft_bag->fast_params[i]->last_state);
                        }
                        /* Set U */ //for(i = 0; i < model_bag->num_models;i++){
                        //       RUN(fill_fast_transitions(model_bag->models[i], ft_bag->fast_params[i]));

                        //      ft_bag->max_last_state = MACRO_MAX(ft_bag->max_last_state,ft_bag->fast_params[i]->last_state);
                        //}
                        RUN(reset_valid_path(sb,model_bag->num_models));
                        RUN(set_u_multi(model_bag, ft_bag, sb));

                        //RUN(set_u(sb,model,ft, &min_u));

                        //exit(0);
                        /* I only want to add states if the last iteration was successful */
                        //if(!no_path){
                        RUN(expand_ihmms(model_bag, ft_bag));
                                //}

                        RUN(resize_spotseq_thread_data(td, &num_threads,(sb->max_len+2)  , ft_bag->max_last_state));

                        /*for(i = 0; i < model_bag->num_models;i++){
                          LOG_MSG("Iteration %d Model %d (%d states)  alpha = %f, gamma = %f", iter,i, model_bag->models[i]->num_states, model_bag->models[i]->alpha ,model_bag->models[i]->gamma);
                          }*/

                        //LOG_MSG("Iteration %d (%d states) sampling %d ", iter, model->num_states,sb->num_seq);
                        //exit(0);
                        //dyn prog + labelling
                        for(i = 0; i < num_threads;i++){
                                td[i]->ft_bag = ft_bag;
                                td[i]->ft = ft;
                                td[i]->sb = sb;
                                if(thr_pool_queue(pool,do_dynamic_programming,td[i]) == -1){
                                        fprintf(stderr,"Adding job to queue failed.");
                                }
                                /* if(thr_pool_queue(local_pool,do_forward_backward,td[i]) == -1){ */
                                /*         fprintf(stderr,"Adding job to queue failed."); */
                                /* } */
                                /* if(thr_pool_queue(local_pool, do_sample_path_and_posterior,td[i]) == -1){ */
                                /*         fprintf(stderr,"Adding job to queue failed."); */
                                /* } */
                        }
                        thr_pool_wait(pool);

                        //LOG_MSG("Check labelling after dyn (%d).",iter);
                        //RUN(check_labels(sb,model_bag->num_models ));
                        //LOG_MSG("Done");
                        no_path = 0;


                        RUN(detect_valid_path(sb,model_bag->num_models, &no_path));
                        if(no_path){
                                LOG_MSG("weird split must have happened. %d",iter);
                                //exit(0);
                                iterations++;
                        }
                }

                /* swap tmp label with label */
                int** tmp = NULL;

                for(i = 0; i < sb->num_seq;i++){
                        tmp= sb->sequences[i]->label_arr;
                        sb->sequences[i]->label_arr = sb->sequences[i]->tmp_label_arr;
                        sb->sequences[i]->tmp_label_arr = tmp;
                }

                /* if more than 1% of sequences don't have a path redo */

                ///if(no_path){
                //        LOG_MSG("weird split must have happened. %d",iter);
                        //exit(0);
                        //for(i = 0; i < model_bag->num_models;i++){
                        //       RUN(fill_fast_transitions(model_bag->models[i], ft_bag->fast_params[i]));

                        //      ft_bag->max_last_state = MACRO_MAX(ft_bag->max_last_state,ft_bag->fast_params[i]->last_state);
                        //}


                        //RUN(fill_fast_transitions(model,ft));
                        //iterations++;
                //}else{
                        /* I am doing this as a pre-caution. I don't want the inital model
                         * contain states that are not visited.. */
                        //ft_bag->max_last_state = -1;
                        //model_bag->max_num_states = -1;
                        //for(i = 0; i < model_bag->num_models;i++){
                        //        RUN(remove_unused_states_labels(model_bag->models[i], sb,i ));
                        //        RUN(fill_counts(model_bag->models[i], sb,i));

                                //RUN(fill_fast_transitions(model_bag->models[i], ft_bag->fast_params[i]));

                                //ft_bag->max_last_state = MACRO_MAX(ft_bag->max_last_state,ft_bag->fast_params[i]->last_state);
                        //        model_bag->max_num_states = MACRO_MAX(model_bag->max_num_states, model_bag->models[i]->num_states);
                                //print_model_parameters(model_bag->models[i]);
                        //}

                //}
                //print_fast_hmm_params(ft_bag->fast_params[0]);
                /* print out model - used for plotting  */
                /*if((iter+1) % 10 == 0){
                //LOG_MSG("print %d\n",iter);
                char tmp_buffer[BUFFER_LEN];
                        snprintf(tmp_buffer,BUFFER_LEN,"model_at_%07d.h5",iter+1);
                        RUN(write_model_hdf5(model,tmp_buffer));
                        RUN(add_background_emission(tmp_buffer,ft->background_emission,ft->L));
                        RUN(run_build_fhmm_file(tmp_buffer));
                        }*/
                for(i = 0; i < model_bag->num_models;i++){
                        //LOG_MSG("Iteration %d Model %d (%d states)  alpha = %f, gamma = %f", iter,i, model_bag->models[i]->num_states, model_bag->models[i]->alpha ,model_bag->models[i]->gamma);
                        model_bag->models[i]->training_iterations++;
                }
        }
        //sb->num_seq = sb->org_num_seq;
        //if(need_local_pool){
        //        thr_pool_destroy(local_pool);
        //}
        //free_spotseq_thread_data(td, num_threads);
        return OK;
ERROR:
        //free_spotseq_thread_data(td, num_threads);
        return FAIL;
}


int detect_valid_path(struct seq_buffer* sb,int num_models, int* no_path)
{
        int i,j;

        *no_path = 0;
        for(i = 0; i < sb->num_seq;i++){
                for(j = 0; j < num_models;j++){
                        if(sb->sequences[i]->has_path[j] == 0){
                                //LOG_MSG("weird split must have happened in seq %d m%d",i,j);
                                *no_path = 1;
                                return OK;
                        }
                }
        }

        return OK;

}

int reset_valid_path(struct seq_buffer* sb,int num_models)
{
        int i,j;
        for(i = 0; i < sb->num_seq;i++){
                for(j = 0; j < num_models;j++){
                        sb->sequences[i]->has_path[j] = 0;
                }
        }
        return OK;
}

void* do_forward_backward(void *threadarg)
{
        struct spotseq_thread_data *data;
        int i,j;
        int num_threads;
        int thread_id;

        float f_score;
        float b_score;
        data = (struct spotseq_thread_data *) threadarg;

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
                                if(i  < 5){
                                        fprintf(stdout,"%d %f (f)\n%d %f (b)\n",i, f_score,i,b_score);
                                }
                                RUN(collect_slice(data, data->sb->sequences[i], f_score));
                        }
                }
        }
        return NULL;
ERROR:
        return NULL;
}


void* do_sample_path_and_posterior(void* threadarg)
{
        struct spotseq_thread_data *data;
        struct ihmm_sequence* seq = NULL;
        int i;
        int num_threads;
        int thread_id;
        float f_score;
        float b_score;
        float r_score;

        data = (struct spotseq_thread_data *) threadarg;

        num_threads = data->num_threads;
        thread_id = data->thread_ID;

        for(i =0; i < data->sb->num_seq;i++){
                if( i% num_threads == thread_id){
                        seq = data->sb->sequences[i];
                        //               LOG_MSG("Thread %d running sequence %d",thread_id, i);
                        //RUN(dynamic_programming(data->dyn,data->ft, seq, data->seed));

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
        struct spotseq_thread_data *data;
        struct ihmm_sequence* s = NULL;
        int i;
        int j;
        int num_threads;
        int thread_id;
        //int safety = 10;
        data = (struct spotseq_thread_data *) threadarg;

        num_threads = data->num_threads;
        thread_id = data->thread_ID;

        for(i =0; i < data->sb->num_seq;i++){

                if( i% num_threads == thread_id){

                        s = data->sb->sequences[i];
                        for(j = 0; j < data->ft_bag->num_models; j++){
                                //LOG_MSG("Run seq: %d M:%d (thread%d)",i,j, data->thread_ID);
                                //s->has_path[j] = 0;
                                //safety = 10;
                                //while(!s->has_path[j]){


                                //if(!s->has_path[j]){
                                        RUN(dynamic_programming_clean(data->ft_bag->fast_params[j],
                                                                      data->dyn,
                                                                      s->seq,
                                                                      s->tmp_label_arr[j],
                                                                      s->u_arr[j],
                                                                      s->seq_len,
                                                                      &s->has_path[j],
                                                                      &data->rndstate));

                                        //}
                        }

                        //LOG_MSG("Thread %d running sequence %d   %f %d",thread_id, i,data->sb->sequences[i]->score,data->seed);
                        //RUN(dynamic_programming(data,i));
                        //

                        /*while(data->sb->sequences[i]->score == -INFINITY){
                          RUN(dynamic_programming(data->dyn,data->ft, data->sb->sequences[i]));
                                }*/
                }
        }
        return NULL;
ERROR:
        return NULL;
}








int expand_ihmms(struct model_bag* model_bag, struct fast_param_bag* ft_bag)
{
        struct ihmm_model* model = NULL;
        struct fast_hmm_param* ft = NULL;
        int i;
        float max;
        float min_u;
        for(i = 0; i < model_bag->num_models;i++){

                min_u = model_bag->min_u[i];
                model = model_bag->models[i];
                ft = ft_bag->fast_params[i];


                RUN(get_max_to_last_state_transition(ft, &max));
                while(max >= min_u && model->num_states < MAX_NUM_STATES && max > 0.0 ){//}sb->max_len){
                        //fprintf(stdout,"ITER: %d Add state! MAX:%f min_U:%f max_len: %d \n",iter , max, min_u,sb->max_len);
                        RUN(add_state_from_fast_hmm_param(model,ft));
                        RUN(get_max_to_last_state_transition(ft, &max));
                        //fprintf(stdout,"MAX:%f min_U:%f\n", max, min_u);
                        //exit(0);
                        //      break;
                }
                RUN(make_flat_param_list(ft));
                //print_fast_hmm_params(ft);
                ft_bag->max_last_state = MACRO_MAX(ft_bag->max_last_state, ft->last_state);

        }
        return OK;
ERROR:
        return FAIL;
}

/* This function assumes (oh no!) that beta has space for an additional
p   g * element */
int add_state_from_fast_hmm_param(struct ihmm_model* ihmm,struct fast_hmm_param* ft)
{
        struct fast_t_item** infinity = NULL;
        struct fast_t_item* tmp = NULL;
        float* tmp_prob = NULL;

        float* beta;
        float alpha;
        float gamma;
        //rk_state rndstate;

        float sum,be,bg,pe,pg, a,b;
        int i,new_k;//,list_index;
        //intl,r;

        //int pg_hack;            /* I don't want add states that are not reachable. */
        //float* tmp_pg = NULL;

        ASSERT(ihmm != NULL, "No model");
        ASSERT(ft != NULL, "No ft.");
        /* Sorting is only strictly necessary if this is called after another function re-sorted it */
        //qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_from_asc);

        //rndstate = ihmm->rndstate;

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
                tmp_prob[i] =  rk_gamma(&ihmm->rndstate, beta[i] * alpha, 1.0);
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
        bg = rk_beta(&ihmm->rndstate, 1.0,gamma );

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
                        pg = rk_binomial(&ihmm->rndstate, 1.0, a / (a+b));
                }else{
                        pg = rk_beta(&ihmm->rndstate, a, b);
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
                ft->emission[i][new_k] = rk_gamma(&ihmm->rndstate, ft->background_emission[i], 1.0);
                sum += ft->emission[i][new_k];
        }
        for(i = 0; i < ihmm->L;i++){
                ft->emission[i][new_k] /= sum;
        }


        //MFREE(tmp_pg);
        //ft->num_items = list_index;
        ft->last_state = new_k+1;
        //ihmm->rndstate = rndstate;
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

int sum_counts_from_multiple_threads(struct spotseq_thread_data** td,int* num_threads,int K)
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




int collect_slice(struct spotseq_thread_data * data,struct ihmm_sequence* ihmm_seq, float total)
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

int dynamic_programming_clean(struct fast_hmm_param* ft,  float** matrix,uint8_t* seq,int* label,float* u,int len,uint8_t* has_path,rk_state* random)
{
        struct fast_t_item** list = NULL;
        int i,j,boundary;
        int state;
        int a,b;
        float sum;
        float* emission;
        float* tmp_row;
        float r, tmp_r;
        int K;

        K = ft->last_state;

        list = ft->list;
        tmp_row = matrix[len];

        boundary = fast_hmm_param_binarySearch_t(ft, u[0]);
        for(i = 0; i < K;i++){
                matrix[0][i] = 0.0f;
        }
        //fill first row.
        for(j = 0; j < boundary;j++){
                if(list[j]->from == IHMM_START_STATE){
                        matrix[0][list[j]->to] = list[j]->t;
                }
        }
        sum = 0;
        emission = ft->emission[seq[0]];
        for(i = 0; i < K;i++){

                matrix[0][i] *=  emission[i];

                sum += matrix[0][i];
        }
        for(i = 0; i < K;i++){
                matrix[0][i] /= sum;
        }
        for(i = 1; i < len;i++){
                emission = ft->emission[seq[i]];
                for(j = 0; j < K;j++){
                        matrix[i][j] = 0.0f;
                }

                boundary = fast_hmm_param_binarySearch_t(ft, u[i]);

                for(j = 0; j < boundary;j++){
                        a = list[j]->from;
                        b = list[j]->to;
                        matrix[i][b] += matrix[i-1][a];

                }
                sum = 0.0;

                for(j = 0; j < K;j++){
                        matrix[i][j] *=  emission[j];
                        sum += matrix[i][j];

                }
                for(j = 0; j < K;j++){
                        matrix[i][j] /= sum;
                }

        }


        sum = 0.0f;
        //float tmp_r;

        boundary = fast_hmm_param_binarySearch_t(ft, u[len]);
        for(j = 0; j < boundary;j++){
                a = list[j]->from;
                b = list[j]->to;
                if(b == IHMM_END_STATE){
                        sum += matrix[len-1][a];
                }
        }

        if(sum != 0.0 && !isnan(sum)){
                state = IHMM_END_STATE;
                for(i = len-1; i >= 0; i--){
                        //fprintf(stdout,"pick: %d %d\n", i,state);
                        for(j = 0; j < K;j++){
                                tmp_row[j] = -1.0;
                        }
                        sum = 0.0f;
                        boundary = fast_hmm_param_binarySearch_t(ft, u[i+1]);
                        for(j = 0; j < boundary;j++){
                                a = list[j]->from;
                                b = list[j]->to;
                                if(b == state && a != IHMM_START_STATE){
                                        tmp_row[a] = matrix[i][a];
                                        sum += matrix[i][a];
                                }
                        }
                        //r =  random_float_zero_to_x(sum);
                        //r = rand_r(&seed) / (float) RAND_MAX *sum;
                        //tmp_r = rk_double(random);
                        r = rk_double(random)*sum;
                        tmp_r = r;
                        //r = random_float_zero_to_x_thread(sum, &data->seed);
                        for(j = 0; j < boundary;j++){
                                //if(j == 0 && i == len-1){
                                //        fprintf(stdout,"%f thread: %f  %f \n",random_float_zero_to_x(sum), random_float_zero_to_x_thread(sum, &seed) , rand_r(&seed) / (float) RAND_MAX);
                                //}
                                a = list[j]->from;
                                b = list[j]->to;
                                if(list[j]->to == state && a != IHMM_START_STATE){
                                        r -= tmp_row[a];

                                        label[i] = a;
                                        if(r <= FLT_EPSILON){
                                                state = a;
                                                break;
                                        }
                                }
                        }

                }
                /* sanitycheck!  */
                *has_path = 1;
        }else{
                *has_path = 0;
                //u[0] = -1.0f;
        }

        return OK;

}

int dynamic_programming(struct spotseq_thread_data* data, int target)
{
        float** matrix = NULL;
        struct fast_hmm_param* ft = NULL;
        struct ihmm_sequence* ihmm_seq = NULL;

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

        ASSERT(data != NULL, "no thread data");

        matrix = data->dyn;
        ft = data->ft;
        ihmm_seq =  data->sb->sequences[target];

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
                        matrix[0][list[j]->to] = list[j]->t;
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
                        //r =  random_float_zero_to_x(sum);

                        //r = rand_r(&seed) / (float) RAND_MAX *sum;
                        r = random_float_zero_to_x_thread(sum, &data->seed);

                        for(j = 0; j < boundary;j++){
                                //if(j == 0 && i == len-1){
                                //        fprintf(stdout,"%f thread: %f  %f \n",random_float_zero_to_x(sum), random_float_zero_to_x_thread(sum, &seed) , rand_r(&seed) / (float) RAND_MAX);
                                //}
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

        }
        return OK;
ERROR:
        return FAIL;
}


int set_u_multi(struct model_bag* model_bag, struct fast_param_bag*  ft_bag, struct seq_buffer* sb)
{
        int i;

        for(i = 0; i < model_bag->num_models;i++){
                RUN(set_u(sb, model_bag->models[i], ft_bag->fast_params[i], &model_bag->min_u[i],i));
        }
        return OK;
ERROR:
        return FAIL;
}

int set_u(struct seq_buffer* sb, struct ihmm_model* model, struct fast_hmm_param* ft, float* min_u, int model_index)
{
        int i,j;
        float* u = 0;
        int* label =0;
        float x;
        //float r;
        int len;

        float local_min_u = 1.0;
        ASSERT(sb != NULL, "No sequences.");
        ASSERT(model != NULL, "No model.");
        //qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_from_asc);
        //last_state = ft->last_state;

        for(i = 0; i < sb->num_seq;i++){
                label = sb->sequences[i]->label_arr[model_index];
                u = sb->sequences[i]->u_arr[model_index];
                len = sb->sequences[i]->seq_len;

                x = ft->transition[IHMM_START_STATE][label[0]];
                //c = IHMM_START_STATE * last_state + label[0];
                //c = a* (num_states-1) + b;
                //u[0] = rk_beta(&model->rndstate, 1.0, 11) * x;
                //r = rk_beta(&model->rndstate, 1.0, 1.0) * x;
                //while(fabs(r-0.0) < FLT_EPSILON ){
                //        r = rk_beta(&model->rndstate, 1.0, 1.1) * x;
                //}
                //u[0] = r;


                u[0] = rk_double(&model->rndstate) *x;
                //ASSERT(ft->list[c]->t != 0.0f,"BAD %d -> %d %f",ft->list[c]->from,ft->list[c]->to,ft->list[c]->t);

                local_min_u = MACRO_MIN(local_min_u, u[0]);
                for (j = 1; j < len;j++){
                        //c = label[j-1] * last_state + label[j];
                        x = ft->transition[label[j-1]][label[j]];
                        //r = rk_beta(&model->rndstate, 1.0, 1.0) * x;
                        //while(fabs(r-0.0) < FLT_EPSILON ){
                        //        r = rk_beta(&model->rndstate, 1.0, 1.1) * x;
                        //}
                        //u[j] = r;
                        //u[j] = rk_beta(&model->rndstate, 1.0, 11) * x;
                        u[j] =  rk_double(&model->rndstate) * x;//rk_double(&model->rndstate) *
                        //if(!i && j < 5){
                        //       fprintf(stdout,"%d->%d %f\n",label[j-1],label[j],ft->list[c]->t );
                        //}
                        //fprintf(stdout,"%d %d  ;; %d %d\n",label[j-1],label[j],ft->list[c]->from ,ft->list[c]->to);
                        local_min_u = MACRO_MIN(local_min_u, u[j]);
                        //ASSERT(ft->list[c]->t != 0.0f,"BAD %d -> %d %f",ft->list[c]->from,ft->list[c]->to,ft->list[c]->t);


                }
                x = ft->transition[label[len-1]][IHMM_END_STATE];

                //r = rk_beta(&model->rndstate, 1.0, 1.0) * x;
                //while(fabs(r-0.0) < FLT_EPSILON ){
                //        r = rk_beta(&model->rndstate, 1.0, 1.1) * x;
                // }
                //u[len] = r;
                //u[len] = rk_beta(&model->rndstate, 1.0, 11) * x;
                u[len] = rk_double(&model->rndstate) * x;//(ft->list[c]->t);
                //ASSERT(ft->list[c]->t != 0.0f,"BAD %d -> %d %f",ft->list[c]->from,ft->list[c]->to,ft->list[c]->t);
                //fprintf(stdout,"%d %d -> %d: %f  \n",label[len-1],ft->list[c]->from ,ft->list[c]->to, ft->list[c]->t );
                local_min_u = MACRO_MIN(local_min_u, u[len]);

        }

        *min_u = local_min_u;
        return OK;
ERROR:
        return FAIL;
}

int reset_u_if_no_path(struct fast_hmm_param* ft, float* u,int * label, int len, rk_state* rndstate)
{
        float x;
        int j;
        x = ft->transition[IHMM_START_STATE][label[0]];
        u[0] = rk_double(rndstate) *x;
        for (j = 1; j < len;j++){
                x = ft->transition[label[j-1]][label[j]];
                u[j] = rk_double(rndstate) * x;
        }
        x = ft->transition[label[len-1]][IHMM_END_STATE];
        u[len] = rk_double(rndstate) * x;
        return OK;
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

        model->alpha_a = 6.0f;
        model->alpha_b = 15.0f;
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

        model->alpha_a = 6.0f;
        model->alpha_b = 15.0f;
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

        RUNP(model = alloc_ihmm_model(initial_states, 4,42));

        RUNP(iseq = create_ihmm_sequences_mem(tmp_seq ,4));
        add_multi_model_label_and_u(iseq, 1);
        RUN(random_label_ihmm_sequences(iseq, initial_states * 10,0.3));
        RUN(fill_counts(model,iseq));
        /* I am doing this as a pre-caution. I don't want the inital model
         * contain states that are not visited.. */
        RUN(remove_unused_states_labels(model, iseq));;
        RUN(fill_counts(model,iseq));
        RUN(print_counts(model));

        model->alpha_a = 4.0f;
        model->alpha_b = 2.0f;
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
        model->alpha_a = 4.0f;
        model->alpha_b = 2.0f;
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

