#include "hmm_conversion.h"

int fill_fast_transitions_only_matrices(struct ihmm_model* model,struct fast_hmm_param* ft)
{
        float* tmp_prob = NULL;
        int i,j;
        //int list_index;
        //int last_index;
        int last_state;
        float sum;

        ASSERT(model != NULL, "No model");
        ASSERT(ft != NULL,"No fast_hmm_param structure");

        // delete old tree...

        if(ft->root){

                free_rbtree(ft->root->node,ft->root->fp_free);
                ft->root->node = NULL;
                ft->root->num_entries = 0;
                ft->root->cur_data_nodes = 0;
                if(ft->root->data_nodes){
                        MFREE(ft->root->data_nodes);
                }
        }

        MMALLOC(tmp_prob, sizeof(float) *(model->num_states));

        last_state = model->num_states -1;
        //fprintf(stdout,"%d last state\n",last_state);
        /* check if there is enough space to hold new transitions... */
        /* This is slightly to generous as I am allocating memory for the
         * infinity state as well */

        RUN(expand_ft_if_necessary(ft, model->num_states));
        //RUN(expand_fast_hmm_param_if_necessary(ft, model->num_states *model->num_states  ));
        /* Empty previous transitions by setting the index to zero  */
        /* Perhaps clear RBtree...  */

        /* Let's start with the transitions from the start state!  */
        //last_index = list_index;
        sum = 0.0;

        /* Disallow Start to start transitions */
        ft->transition[IHMM_START_STATE][IHMM_START_STATE] = 0.0;
        /* Disallow Start to end transitions i.e. zero length sequences are not allowed*/
        ft->transition[IHMM_START_STATE][IHMM_END_STATE] = 0.0f;
        /* Now to the remaining existing transitions... */
        for(i = 2; i < last_state;i++){
                tmp_prob[i] = rk_gamma(&model->rndstate, model->transition_counts[IHMM_START_STATE][i] + model->beta[i] * model->alpha,1.0);
                sum += tmp_prob[i];
        }
        /* the last to infinity transition (this is just used in stick breaking
         * when adding states ). here there should be no counts as this
         * possibility was not observed in the last transition. */
        tmp_prob[last_state] = rk_gamma(&model->rndstate, model->beta[last_state] * model->alpha,1.0);
        sum += tmp_prob[last_state];

        /* Normalize!  */
        for(i = 2; i < last_state;i++){
                ft->transition[IHMM_START_STATE][i] =  tmp_prob[i] / sum;
        }
        ft->transition[IHMM_START_STATE][last_state] = tmp_prob[last_state] / sum;

        /* And now the stop state...  */
        /* There is no possibility escape the end state - all transitions from
         * end are zero. I am not sure if this matters in my dyn prig. code but
         * why not! */
        for(i = 0; i < last_state;i++){
                ft->transition[IHMM_END_STATE][i] = 0.0f;
        }
        ft->transition[IHMM_END_STATE][last_state] = 0.0f;

        for(i = 2; i < last_state;i++){
                /* Remeber where I started filling...  */
                //last_index = list_index;
                sum = 0.0;

                for(j = 1; j < last_state;j++){
                        tmp_prob[j] = rk_gamma(&model->rndstate, model->transition_counts[i][j] + model->beta[j] * model->alpha,1.0);
                        sum += tmp_prob[j];
                }
                tmp_prob[last_state] = rk_gamma(&model->rndstate, model->beta[last_state] * model->alpha,1.0);
                sum += tmp_prob[last_state];

                for(j = 1; j < last_state;j++){
                        ft->transition[i][j] = tmp_prob[j] / sum;
                }
                ft->transition[i][last_state] = tmp_prob[last_state] / sum;
        }
        /* init emission probabilities.  */
        /*  model->emission[i][j]  = rk_gamma(&model->rndstate, model->emission[i][j]+ EMISSION_H, 1.0);
        //model->emission[i][j]  =  model->emission[i][j]+ EMISSION_H;
        sum +=model->emission[i][j];
        }
        for(j = 0;j <model->L;j++){
        model->emission[i][j] /= sum;*/
        for(i = 0; i < model->L;i++){
                for(j = 0; j < last_state;j++){
                        ft->emission[i][j] = rk_gamma(&model->rndstate, model->emission_counts[i][j] + EMISSION_H, 1.0);
                }
        }

        for(j = 0; j < last_state;j++){
                sum = 0.0f;
                for(i = 0; i < model->L;i++){
                        sum += ft->emission[i][j];
                }
                for(i = 0; i < model->L;i++){
                        ft->emission[i][j] /= sum;
                }
        }

        /* kind of important... */
        ft->last_state = last_state;
        MFREE(tmp_prob);
        return OK;
ERROR:
        return FAIL;
}

int fill_fast_transitions(struct ihmm_model* model,struct fast_hmm_param* ft)
{
        struct fast_t_item* tmp = NULL;
        struct fast_t_item** infinity = NULL;
        float* tmp_prob = NULL;
        int i,j;
        //int list_index;
        //int last_index;
        int last_state;
        float sum;

        ASSERT(model != NULL, "No model");
        ASSERT(ft != NULL,"No fast_hmm_param structure");

        // delete old tree...

        if(ft->root){
                free_rbtree(ft->root->node,ft->root->fp_free);
                ft->root->node = NULL;
                ft->root->num_entries = 0;
                ft->root->cur_data_nodes = 0;
                if(ft->root->data_nodes){
                        MFREE(ft->root->data_nodes);
                }
        }

        MMALLOC(tmp_prob, sizeof(float) *(model->num_states));
        last_state = model->num_states -1;

        /* check if there is enough space to hold new transitions... */
        /* This is slightly to generous as I am allocating memory for the
         * infinity state as well */

        RUN(expand_ft_if_necessary(ft, model->num_states));
        //RUN(expand_fast_hmm_param_if_necessary(ft, model->num_states *model->num_states  ));
        /* Empty previous transitions by setting the index to zero  */
        /* Perhaps clear RBtree...  */
        infinity = ft->infinity;

        /* Let's start with the transitions from the start state!  */
        //last_index = list_index;
        sum = 0.0;

        /* Disallow Start to start transitions */

        tmp = NULL;
        MMALLOC(tmp, sizeof(struct fast_t_item));
        tmp->from = IHMM_START_STATE;
        tmp->to = IHMM_START_STATE;
        tmp->t =  0.0f;
        ft->root->tree_insert(ft->root,tmp);
        // insert into transition matrix.

        ft->transition[IHMM_START_STATE][IHMM_START_STATE] = 0.0;


        /* Disallow Start to end transitions i.e. zero length sequences are not allowed*/
        tmp = NULL;
        MMALLOC(tmp, sizeof(struct fast_t_item));
        tmp->from = IHMM_START_STATE;
        tmp->to = IHMM_END_STATE;
        tmp->t =  0.0f;
        ft->root->tree_insert(ft->root,tmp);

        ft->transition[IHMM_START_STATE][IHMM_END_STATE] = 0.0f;
        /* Now to the remaining existing transitions... */
        for(i = 2; i < last_state;i++){
                tmp_prob[i] = rk_gamma(&model->rndstate, model->transition_counts[IHMM_START_STATE][i] + model->beta[i] * model->alpha,1.0);
                sum += tmp_prob[i];
        }
        /* the last to infinity transition (this is just used in stick breaking
         * when adding states ). here there should be no counts as this
         * possibility was not observed in the last transition. */
        tmp_prob[last_state] = rk_gamma(&model->rndstate, model->beta[last_state] * model->alpha,1.0);
        sum += tmp_prob[last_state];

        /* Normalize!  */
        for(i = 2; i < last_state;i++){
                tmp = NULL;
                MMALLOC(tmp, sizeof(struct fast_t_item));
                tmp->from = IHMM_START_STATE;
                tmp->to = i;
                tmp->t =  tmp_prob[i] / sum;
                ft->root->tree_insert(ft->root,tmp);
                ft->transition[IHMM_START_STATE][i] = tmp->t;
        }
        infinity[IHMM_START_STATE]->from = IHMM_START_STATE;
        infinity[IHMM_START_STATE]->to = last_state;
        infinity[IHMM_START_STATE]->t = tmp_prob[last_state] / sum;;
        ft->transition[IHMM_START_STATE][last_state] = infinity[IHMM_START_STATE]->t;


        /* And now the stop state...  */
        /* There is no possibility escape the end state - all transitions from
         * end are zero. I am not sure if this matters in my dyn prig. code but
         * why not! */
        for(i = 0; i < last_state;i++){
                tmp = NULL;
                MMALLOC(tmp, sizeof(struct fast_t_item));
                tmp->from = IHMM_END_STATE;
                tmp->to = i;
                tmp->t = 0.0f;
                ft->root->tree_insert(ft->root,tmp);
                ft->transition[IHMM_END_STATE][i] = 0.0f;
        }
        infinity[IHMM_END_STATE]->from = IHMM_END_STATE;
        infinity[IHMM_END_STATE]->to = last_state;
        infinity[IHMM_END_STATE]->t = 0.0f;
        ft->transition[IHMM_END_STATE][last_state] = 0.0f;


        for(i = 2; i < last_state;i++){
                /* Remeber where I started filling...  */
                //last_index = list_index;
                sum = 0.0;

                for(j = 1; j < last_state;j++){
                        tmp_prob[j] = rk_gamma(&model->rndstate, model->transition_counts[i][j] + model->beta[j] * model->alpha,1.0);
                        sum += tmp_prob[j];
                }
                tmp_prob[last_state] = rk_gamma(&model->rndstate, model->beta[last_state] * model->alpha,1.0);
                sum += tmp_prob[last_state];

                for(j = 1; j < last_state;j++){
                        tmp = NULL;
                        MMALLOC(tmp, sizeof(struct fast_t_item));
                        tmp->from = i;
                        tmp->to = j;
                        tmp->t = tmp_prob[j] / sum;
                        ft->root->tree_insert(ft->root,tmp);
                        ft->transition[i][j] = tmp->t;
                }
                infinity[i]->from = i;
                infinity[i]->to = last_state;
                infinity[i]->t = tmp_prob[last_state] / sum;
                ft->transition[i][last_state] = infinity[i]->t;
        }
        /* init emission probabilities.  */
        /*  model->emission[i][j]  = rk_gamma(&model->rndstate, model->emission[i][j]+ EMISSION_H, 1.0);
        //model->emission[i][j]  =  model->emission[i][j]+ EMISSION_H;
        sum +=model->emission[i][j];
        }
        for(j = 0;j <model->L;j++){
        model->emission[i][j] /= sum;*/
        for(i = 0; i < model->L;i++){
                for(j = 0; j < last_state;j++){
                        ft->emission[i][j] = rk_gamma(&model->rndstate, model->emission_counts[i][j] + EMISSION_H, 1.0);
                }
        }


        for(j = 0; j < last_state;j++){
                sum = 0.0f;
                for(i = 0; i < model->L;i++){
                        sum += ft->emission[i][j];
                }
                for(i = 0; i < model->L;i++){
                        ft->emission[i][j] /= sum;
                }
        }

        /* kind of important... */
        ft->last_state = last_state;
        MFREE(tmp_prob);
        return OK;
ERROR:
        return FAIL;
}

struct fhmm* build_finite_hmm_from_infinite_hmm(struct ihmm_model* model)
{
        struct fhmm* fhmm = NULL;
        struct fast_hmm_param* ft = NULL;

        float** s1_e = NULL;
        float** s2_e = NULL;

        float** s1_t = NULL;
        float** s2_t = NULL;

        int initial_states = 10;
        float sum;
        int iterations = 1000;
        int iter,c,i,j;

        ASSERT(model != NULL, "No model");

        RUNP(fhmm = alloc_fhmm());



        fhmm->K = model->num_states;
        fhmm->L = model->L;

        RUNP(ft = alloc_fast_hmm_param(initial_states,model->L));

        RUN(fill_background_emission_from_model(ft,model));

        /* copy background probabilitys into fhmm */

        MMALLOC(fhmm->background, sizeof(float) * fhmm->L);
        for (i = 0; i < fhmm->L; i++){
                fhmm->background[i] = ft->background_emission[i];
        }

        /* Set alphabet and number of states. */

        s1_e = malloc_2d_float(s1_e, model->num_states, model->L, 0.0);
        s2_e = malloc_2d_float(s2_e, model->num_states, model->L, 0.0);

        s1_t = malloc_2d_float(s1_t, model->num_states, model->num_states, 0.0);
        s2_t = malloc_2d_float(s2_t, model->num_states, model->num_states, 0.0);

        for(iter=  0;iter < iterations;iter++){
                RUN(fill_fast_transitions_only_matrices(model,ft));
                for(i = 0;i < model->num_states;i++){
                        for(c = 0; c < model->L;c++){
                                s1_e[i][c] += ft->emission[c][i];
                                s2_e[i][c] += (ft->emission[c][i] * ft->emission[c][i]);
                        }
                }
                for(i = 0;i < model->num_states;i++){
                        for(c = 0;c < model->num_states;c++){
                                s1_t[i][c] += ft->transition[i][c];
                                s2_t[i][c] += (ft->transition[i][c] * ft->transition[i][c]);
                        }
                }

        }
        for(i = 0; i < model->num_states;i++){
                sum = 0.0;
                for(j = 0; j < model->L;j++){
                        s2_e[i][j] = sqrt(  ((double) iterations * s2_e[i][j] - s1_e[i][j] * s1_e[i][j])/ ((double) iterations * ((double) iterations -1.0)));
                        s1_e[i][j] = s1_e[i][j] / (double) iterations;
                        sum+= s1_e[i][j];
                        //fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_e[i][j], s2_e[i][j]);
                }
                if(sum){
                        //fprintf(stdout,"Emission:%d\n",i);
                        for(j = 0; j < model->L;j++){
                                s1_e[i][j] /= sum;
                                //fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_e[i][j], s2_e[i][j]);
                        }
                }

                sum = 0;
                for(j = 0;j < model->num_states;j++){
                        s2_t[i][j] = sqrt(  ((double) iterations * s2_t[i][j] - s1_t[i][j] * s1_t[i][j])/ ((double) iterations * ((double) iterations -1.0)));
                        s1_t[i][j] = s1_t[i][j] / (double) iterations;
                        sum+= s1_t[i][j];
                        //fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_t[i][j], s2_t[i][j]);
                }
                if(sum){
                        //fprintf(stdout,"transition:%d\n",i);
                        for(j = 0;j < model->num_states;j++){
                                s1_t[i][j] /= sum;
                                //fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_t[i][j], s2_t[i][j]);
                        }
                }
        }

        fhmm->e = s1_e;
        fhmm->t = s1_t;
        /* alloc dyn matrices (now that I know how many states there are) */
        RUN(alloc_dyn_matrices(fhmm));

        /* convert probs into log space/ set tindex to allow for fast-ish dyn
         * programming in case there is a sparse transition matrix */
        RUN(setup_model(fhmm));

        free_2d((void**) s2_e);

        free_2d((void**) s2_t);
        free_fast_hmm_param(ft);
        return fhmm;
ERROR:
        free_2d((void**) s1_e);
        free_2d((void**) s2_e);

        free_2d((void**) s1_t);
        free_2d((void**) s2_t);
        free_fast_hmm_param(ft);

        return NULL;
}

int run_build_fhmm_file(char* h5file)
{
        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;

        int initial_states = 10;
        int iter;
        int i,j,c;
        float** s1_e = NULL;
        float** s2_e = NULL;

        float** s1_t = NULL;
        float** s2_t = NULL;
        float sum;
        int iterations = 1000;

        ASSERT(h5file != NULL, "No parameters found.");


        RUNP(model = read_model_hdf5(h5file));
        RUNP(ft = alloc_fast_hmm_param(initial_states,model->L));


        /* first index is state * letter ; second is sample (max = 100) */

        s1_e = malloc_2d_float(s1_e, model->num_states, model->L, 0.0);
        s2_e = malloc_2d_float(s2_e, model->num_states, model->L, 0.0);

        s1_t = malloc_2d_float(s1_t, model->num_states, model->num_states, 0.0);
        s2_t = malloc_2d_float(s2_t, model->num_states, model->num_states, 0.0);


        for( iter=  0;iter < iterations;iter++){
                RUN(fill_fast_transitions_only_matrices(model,ft));
                for(i = 0;i < model->num_states;i++){
                        for(c = 0; c < model->L;c++){
                                s1_e[i][c] += ft->emission[c][i];
                                s2_e[i][c] += (ft->emission[c][i] * ft->emission[c][i]);
                        }
                }
                for(i = 0;i < model->num_states;i++){
                        for(c = 0;c < model->num_states;c++){
                                s1_t[i][c] += ft->transition[i][c];
                                s2_t[i][c] += (ft->transition[i][c] * ft->transition[i][c]);
                        }
                }

        }

        for(i = 0; i < model->num_states;i++){
                sum = 0.0;
                for(j = 0; j < model->L;j++){
                        s2_e[i][j] = sqrt(  ((double) iterations * s2_e[i][j] - s1_e[i][j] * s1_e[i][j])/ ((double) iterations * ((double) iterations -1.0)));
                        s1_e[i][j] = s1_e[i][j] / (double) iterations;
                        sum+= s1_e[i][j];
                        //fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_e[i][j], s2_e[i][j]);
                }
                if(sum){
                        //fprintf(stdout,"Emission:%d\n",i);
                        for(j = 0; j < model->L;j++){
                                s1_e[i][j] /= sum;
                                //fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_e[i][j], s2_e[i][j]);
                        }
                }

                sum = 0;
                for(j = 0;j < model->num_states;j++){
                        s2_t[i][j] = sqrt(  ((double) iterations * s2_t[i][j] - s1_t[i][j] * s1_t[i][j])/ ((double) iterations * ((double) iterations -1.0)));
                        s1_t[i][j] = s1_t[i][j] / (double) iterations;
                        sum+= s1_t[i][j];
                        //fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_t[i][j], s2_t[i][j]);
                }
                if(sum){
                        //fprintf(stdout,"transition:%d\n",i);
                        for(j = 0;j < model->num_states;j++){
                                s1_t[i][j] /= sum;
                                //fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_t[i][j], s2_t[i][j]);
                        }
                }
        }

        RUN(add_fhmm(h5file,s1_t,s1_e, model->num_states, model->L  ));

        free_2d((void**) s1_e);
        free_2d((void**) s2_e);

        free_2d((void**) s1_t);
        free_2d((void**) s2_t);

        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        return OK;
ERROR:

        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        return FAIL;
}

int fill_background_emission_from_model(struct fast_hmm_param*ft, struct ihmm_model* model)
{


        int i,j;
        float sum = 0.0f;
        ASSERT(ft != NULL, "No parameters");
        ASSERT(model != NULL, "No model ");
        for(i = 0; i < ft->L;i++){
                ft->background_emission[i] = 0.0f;
        }

        for(i = 0; i < model->L;i++){
                for(j = 0 ; j < model->num_states;j++){
                        ft->background_emission[i] += model->emission_counts[i][j];
                        sum +=  model->emission_counts[i][j];
                }
        }

        ASSERT(sum != 0.0f,"No sequence counts found");
        for(i = 0; i < ft->L;i++){
                ft->background_emission[i] /= sum;
        }

        return OK;
ERROR:
        return FAIL;
}

int fill_background_emission(struct fast_hmm_param*ft,struct seq_buffer* sb)
{

        int i,j;
        float sum = 0.0f;

        ASSERT(ft != NULL, "No parameters");
        ASSERT(sb != NULL, "No sequences");

        for(i = 0; i < ft->L;i++){
                ft->background_emission[i] = 0.0f;
        }

        for(i = 0; i < sb->num_seq;i++){
                for(j = 0;j < sb->sequences[i]->seq_len;j++){
                        ft->background_emission[sb->sequences[i]->seq[j]]++;
                        sum++;
                }
        }
        ASSERT(sum != 0.0f,"No sequence counts found");
        for(i = 0; i < ft->L;i++){
                ft->background_emission[i] /= sum;
        }

        return OK;
ERROR:
        return FAIL;
}
