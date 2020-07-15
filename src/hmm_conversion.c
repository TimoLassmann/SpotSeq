#include <omp.h>

#include "hmm_conversion.h"


#include "randomkit_tl_add.h"
int convert_ihmm_to_fhmm(struct ihmm_model* model,struct fhmm* fhmm, int allow_zero_counts );

int fill_fast_transitions_only_matrices(struct ihmm_model* model,struct fast_hmm_param* ft)
{
        double* tmp_prob = NULL;
        int i,j;
        //int list_index;
        //int last_index;
        int last_state;
        double sum;
        rk_state local_rk_copy;
        ASSERT(model != NULL, "No model");
        ASSERT(ft != NULL,"No fast_hmm_param structure");

        RUN(copy_rk_state(&model->rndstate,&local_rk_copy));
        // delete old tree...
        /*if(ft->root){
                //ft->root->free_tree(ft->root);

                free_rbtree(ft->root->node,ft->root->fp_free);
                ft->root->node = NULL;
                ft->root->num_entries = 0;
                ft->root->cur_data_nodes = 0;
                if(ft->root->data_nodes){
                        MFREE(ft->root->data_nodes);
                }
                }*/

        MMALLOC(tmp_prob, sizeof(double) *(model->num_states));

        last_state = model->num_states -1;
        //fprintf(stdout,"%d last state\n",last_state);
        /* check if there is enough space to hold new transitions... */
        /* This is slightly to generous as I am allocating memory for the
         * infinity state as well */
        RUN(expand_ft_if_necessary(ft, model->num_states));

        //LOG_MSG("size: %d %d",DIM1(ft->transition),DIM2(ft->transition));
        //RUN(expand_fast_hmm_param_if_necessary(ft, model->num_states *model->num_states  ));
        /* Empty previous transitions by setting the index to zero  */
        /* Perhaps clear RBtree...  */

        /* Let's start with the transitions from the start state!  */
        //last_index = list_index;

        //memcpy((void *)local_rk_copy, (void *)&model->rndstate, sizeof (struct rk_state_));

        sum = 0.0;

        /* Disallow Start to start transitions */
        ft->transition[START_STATE][START_STATE] = 0.0;
        /* Disallow Start to end transitions i.e. zero length sequences are not allowed*/
        ft->transition[START_STATE][END_STATE] = 0.0;
        /* Now to the remaining existing transitions... */
        for(i = 2; i < last_state;i++){
                tmp_prob[i] = rk_gamma(&local_rk_copy, model->transition_counts[START_STATE][i] + model->beta[i] * model->alpha,1.0);
                sum += tmp_prob[i];
        }
        /* the last to infinity transition (this is just used in stick breaking
         * when adding states ). here there should be no counts as this
         * possibility was not observed in the last transition. */
        tmp_prob[last_state] = rk_gamma(&local_rk_copy, model->beta[last_state] * model->alpha,1.0);
        sum += tmp_prob[last_state];

        /* Normalize!  */
        for(i = 2; i < last_state;i++){
                ft->transition[START_STATE][i] =  tmp_prob[i] / sum;
        }
        ft->transition[START_STATE][last_state] = tmp_prob[last_state] / sum;

        /* And now the stop state...  */
        /* There is no possibility escape the end state - all transitions from
         * end are zero. I am not sure if this matters in my dyn prig. code but
         * why not! */
        for(i = 0; i < last_state;i++){
                ft->transition[END_STATE][i] = 0.0;
        }
        ft->transition[END_STATE][last_state] = 0.0;

        for(i = 2; i < last_state;i++){
                /* Remeber where I started filling...  */
                //last_index = list_index;
                sum = 0.0;

                for(j = 1; j < last_state;j++){
                        tmp_prob[j] = rk_gamma(&local_rk_copy, model->transition_counts[i][j] + model->beta[j] * model->alpha,1.0);
                        sum += tmp_prob[j];
                }
                tmp_prob[last_state] = rk_gamma(&local_rk_copy, model->beta[last_state] * model->alpha,1.0);
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
                        ft->emission[i][j] = rk_gamma(&local_rk_copy, model->emission_counts[i][j] + EMISSION_H, 1.0);
                }
        }

        for(j = 0; j < last_state;j++){
                sum = 0.0;
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
        double* tmp_prob = NULL;
        int i,j;
        //int list_index;
        //int last_index;
        int last_state;
        double sum;

        ASSERT(model != NULL, "No model");
        ASSERT(ft != NULL,"No fast_hmm_param structure");

        // delete old tree...
        if(!ft->list){
                ft->alloc_num_trans = 65536;
                MMALLOC(ft->list, sizeof(struct fast_t_item*) * ft->alloc_num_trans);
                for(i = 0; i < ft->alloc_num_trans;i++){
                        ft->list[i] = NULL;
                        MMALLOC(ft->list[i], sizeof(struct fast_t_item));
                }
        }
        ft->num_trans = 0;
        /*
        if(ft->root){
                //ft->root->free_tree(ft->root);
                free_rbtree(ft->root->node,ft->root->fp_free);
                ft->root->node = NULL;
                ft->root->num_entries = 0;
                ft->root->cur_data_nodes = 0;
                if(ft->root->data_nodes){
                        MFREE(ft->root->data_nodes);
                }
                }*/

        MMALLOC(tmp_prob, sizeof(double) *(model->num_states));
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

        tmp = ft->list[ft->num_trans];
                //MMALLOC(tmp, sizeof(struct fast_t_item));
        tmp->from = START_STATE;
        tmp->to = START_STATE;
        tmp->t =  0.0;

        ft->num_trans++;
        if(ft->num_trans == ft->alloc_num_trans){
                RUN(expand_num_trans(ft));
        }
        //ft->root->tree_insert(ft->root,tmp);
        // insert into transition matrix.

        ft->transition[START_STATE][START_STATE] = 0.0;


        /* Disallow Start to end transitions i.e. zero length sequences are not allowed*/
        tmp = ft->list[ft->num_trans];
                //MMALLOC(tmp, sizeof(struct fast_t_item));
        tmp->from = START_STATE;
        tmp->to = END_STATE;
        tmp->t =  0.0;
        ft->num_trans++;
        if(ft->num_trans == ft->alloc_num_trans){
                RUN(expand_num_trans(ft));
        }
        //ft->root->tree_insert(ft->root,tmp);

        ft->transition[START_STATE][END_STATE] = 0.0;
        /* Now to the remaining existing transitions... */
        for(i = 2; i < last_state;i++){
                tmp_prob[i] = rk_gamma(&model->rndstate, model->transition_counts[START_STATE][i] + model->beta[i] * model->alpha,1.0);
                sum += tmp_prob[i];
        }
        /* the last to infinity transition (this is just used in stick breaking
         * when adding states ). here there should be no counts as this
         * possibility was not observed in the last transition. */
        tmp_prob[last_state] = rk_gamma(&model->rndstate, model->beta[last_state] * model->alpha,1.0);
        sum += tmp_prob[last_state];

        /* Normalize!  */
        for(i = 2; i < last_state;i++){
                tmp = ft->list[ft->num_trans];
                //MMALLOC(tmp, sizeof(struct fast_t_item));
                tmp->from = START_STATE;
                tmp->to = i;
                tmp->t =  tmp_prob[i] / sum;
                ft->num_trans++;
                if(ft->num_trans == ft->alloc_num_trans){
                        RUN(expand_num_trans(ft));
                }
                //ft->root->tree_insert(ft->root,tmp);
                ft->transition[START_STATE][i] = tmp->t;
        }
        infinity[START_STATE]->from = START_STATE;
        infinity[START_STATE]->to = last_state;
        infinity[START_STATE]->t = tmp_prob[last_state] / sum;;
        ft->transition[START_STATE][last_state] = infinity[START_STATE]->t;


        /* And now the stop state...  */
        /* There is no possibility escape the end state - all transitions from
         * end are zero. I am not sure if this matters in my dyn prig. code but
         * why not! */
        for(i = 0; i < last_state;i++){
                tmp = ft->list[ft->num_trans];
                //MMALLOC(tmp, sizeof(struct fast_t_item));
                tmp->from = END_STATE;
                tmp->to = i;
                tmp->t = 0.0;
                ft->num_trans++;
                if(ft->num_trans == ft->alloc_num_trans){
                        RUN(expand_num_trans(ft));
                }
                //ft->root->tree_insert(ft->root,tmp);
                ft->transition[END_STATE][i] = 0.0;
        }
        infinity[END_STATE]->from = END_STATE;
        infinity[END_STATE]->to = last_state;
        infinity[END_STATE]->t = 0.0;
        ft->transition[END_STATE][last_state] = 0.0;


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
                if(sum == 0.0){
                        sum = 1.0;
                }

                for(j = 1; j < last_state;j++){
                        tmp = ft->list[ft->num_trans];
                        //MMALLOC(tmp, sizeof(struct fast_t_item));
                        tmp->from = i;
                        tmp->to = j;
                        tmp->t = tmp_prob[j] / sum;
                        ft->num_trans++;
                        if(ft->num_trans == ft->alloc_num_trans){
                                RUN(expand_num_trans(ft));
                        }
                        //ft->root->tree_insert(ft->root,tmp);
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
                sum = 0.0;
                for(i = 0; i < model->L;i++){
                        sum += ft->emission[i][j];
                }
                if(sum == 0.0){
                        sum = 1.0;
                }
                for(i = 0; i < model->L;i++){
                        ft->emission[i][j] /= sum;
                        //fprintf(stdout,"%d %d: %f\n", i,j, ft->emission[i][j]);
                }

        }

        /* kind of important... */
        ft->last_state = last_state;

        //ft->root->print_tree(ft->root,stdout);


        MFREE(tmp_prob);
        return OK;
ERROR:
        return FAIL;
}


struct fhmm* build_finite_hmm_from_infinite_hmm(struct ihmm_model* model)
{
        struct fhmm* fhmm = NULL;
        struct fast_hmm_param* ft = NULL;

        double** s1_e = NULL;
        double** s2_e = NULL;

        double** s1_t = NULL;
        double** s2_t = NULL;

        int initial_states = 10;
        double sum;
        int iterations = 1000;
        int iter,c,i,j;

        int* used = NULL;
        int local_num_states = 0;

        //model->alloc_num_states;


        ASSERT(model != NULL, "No model");

        MMALLOC(used, sizeof(int)* model->num_states);

        for(i = 0; i < model->num_states;i++){
                used[i] = -1;
        }
        local_num_states = 0;
        used[0] = local_num_states;
        local_num_states++;
        used[1] = local_num_states;
        local_num_states++;

        for(j = 2; j < model->num_states-1;j++){
                sum = 0.0;
                for(i = 0; i < model->L;i++){
                        sum += model->emission_counts[i][j];
                }
                //fprintf(stdout,"%f ",sum);
                if(sum){
                        used[j] = local_num_states;
                        local_num_states++;
                }
        }
        //fprintf(stdout,"\n");

        /* for(i = 0; i < model->num_states;i++){ */
        /*         fprintf(stdout,"%d ",used[i]); */
        /* } */
        /* fprintf(stdout,"\n"); */

        RUNP(fhmm = alloc_fhmm());

        fhmm->alloc_K = model->alloc_num_states;
        fhmm->K = local_num_states;//model->num_states;
        fhmm->L = model->L;

        RUNP(ft = alloc_fast_hmm_param(initial_states,model->L));

        //RUN(fill_background_emission_from_model(ft,model));

        /* copy background probabilitys into fhmm */

        //MMALLOC(fhmm->background, sizeof(double) * fhmm->L);
        RUN(galloc(&fhmm->background, fhmm->L));
        for (i = 0; i < fhmm->L; i++){
                fhmm->background[i] = (double) model->background[i];//  ft->background_emission[i];
        }

        /* Note: there is a possibility that un-visited states exist
         * in the ihmm. These will have to be removed from the
         * fhmm. */

        /* Set alphabet and number of states. */

        /*RUNP(s1_e = galloc(s1_e, local_num_states, model->L, 0.0));
        RUNP(s2_e = galloc(s2_e, local_num_states, model->L, 0.0));

        RUNP(s1_t = galloc(s1_t, local_num_states, local_num_states, 0.0));
        RUNP(s2_t = galloc(s2_t, local_num_states, local_num_states, 0.0));
        */
        //LOG_MSG("States alloc (max) : %d", model->alloc_num_states);
        RUN(galloc(&s1_e, fhmm->alloc_K, model->L));
        RUN(galloc(&s2_e, fhmm->alloc_K, model->L));
        RUN(galloc(&s1_t, fhmm->alloc_K, fhmm->alloc_K));
        RUN(galloc(&s2_t, fhmm->alloc_K, fhmm->alloc_K));
        for(i = 0; i < fhmm->alloc_K;i++){
                for(j = 0;j < model->L;j++){
                        s1_e[i][j] = 0.0;
                        s2_e[i][j] = 0.0;
                }
                for(j = 0;j < fhmm->alloc_K;j++){
                        s1_t[i][j] = 0.0;
                        s2_t[i][j] = 0.0;
                }
        }


        for(iter=  0;iter < iterations;iter++){
                RUN(fill_fast_transitions_only_matrices(model,ft));
                for(i = 0;i < model->num_states;i++){
                        if(used[i] != -1){
                                for(c = 0; c < model->L;c++){
                                        s1_e[used[i]][c] +=  ft->emission[c][i];
                                        s2_e[used[i]][c] += (ft->emission[c][i] * ft->emission[c][i]);
                                }
                        }
                }
                for(i = 0;i < model->num_states;i++){
                        if(used[i] != -1){
                                for(c = 0;c < model->num_states;c++){
                                        if(used[c] != -1){
                                                s1_t[used[i]][used[c]] += ft->transition[i][c];
                                                s2_t[used[i]][used[c]] += (ft->transition[i][c] * ft->transition[i][c]);

                                        }
                                }
                        }
                }
        }

        for(i = 0; i < local_num_states;i++){
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
                for(j = 0;j < local_num_states;j++){
                        //fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_t[i][j], s2_t[i][j]);
                        s2_t[i][j] = sqrt(  ((double) iterations * s2_t[i][j] - s1_t[i][j] * s1_t[i][j])/ ((double) iterations * ((double) iterations -1.0)));
                        s1_t[i][j] = s1_t[i][j] / (double) iterations;
                        sum+= s1_t[i][j];

                }
                if(sum){
                        //fprintf(stdout,"transition:%d\n",i);
                        for(j = 0;j < local_num_states;j++){
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

        MFREE(used);
        gfree(s2_e);
        gfree(s2_t);
        free_fast_hmm_param(ft);
        return fhmm;
ERROR:
        gfree(s1_e);
        gfree(s1_t);

        gfree(s2_e);
        gfree(s2_t);

        free_fast_hmm_param(ft);

        return NULL;
}

int convert_ihmm_to_fhmm_models(struct model_bag* model_bag)
{
        int i;

        int miter;
        ASSERT(model_bag != NULL," No models");

        if(model_bag->finite_models){
                //WARNING_MSG("Finite hmm models exist - will over-write.");
                for(i = 0; i < model_bag->num_models;i++){
                        free_fhmm(model_bag->finite_models[i]);
                        model_bag->finite_models[i] = NULL;
                }
                MFREE(model_bag->finite_models);
                model_bag->finite_models = NULL;
        }
        //ASSERT(model_bag->finite_models == NULL, "Warning fhmm models already exist?");

        MMALLOC(model_bag->finite_models, sizeof(struct fhmm*)* model_bag->num_models);

#ifdef HAVE_OPENMP

#pragma omp parallel shared(model_bag) private(miter)
        {
#pragma omp for schedule(dynamic) nowait
#endif
                for(miter = 0; miter < model_bag->num_models;miter++){
                        //LOG_MSG("Looking at model: %d ",miter);
                        model_bag->finite_models[miter] = NULL;

                        model_bag->finite_models[miter] = build_finite_hmm_from_infinite_hmm(model_bag->models[miter]);
                }
#ifdef HAVE_OPENMP
        }
#endif


        return OK;
ERROR:
        return FAIL;
}

int run_build_fhmm_file(char* h5file, int allow_zero_counts)
{
        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;

        int initial_states = 10;
        int iter;
        int i,j,c;
        double** s1_e = NULL;
        double** s2_e = NULL;

        double** s1_t = NULL;
        double** s2_t = NULL;
        double sum;
        int iterations = 1000;

        ASSERT(h5file != NULL, "No parameters found.");

        //RUNP(model = read_model_hdf5(h5file));

        RUNP(ft = alloc_fast_hmm_param(initial_states,model->L));

        /* first index is state * letter ; second is sample (max = 100) */
        /*RUNP(s1_e = galloc(s1_e, model->num_states, model->L, 0.0));
        RUNP(s2_e = galloc(s2_e, model->num_states, model->L, 0.0));

        RUNP(s1_t = galloc(s1_t, model->num_states, model->num_states, 0.0));
        RUNP(s2_t = galloc(s2_t, model->num_states, model->num_states, 0.0));
        */
        RUN(galloc(&s1_e, model->num_states, model->L));
        RUN(galloc(&s2_e, model->num_states, model->L));
        RUN(galloc(&s1_t, model->num_states, model->num_states));
        RUN(galloc(&s2_t, model->num_states, model->num_states));
        for(i = 0; i < model->num_states;i++){
                for(j = 0;j < model->L;j++){
                        s1_e[i][j] = 0.0;
                        s2_e[i][j] = 0.0;
                }
                for(j = 0;j < model->num_states;j++){
                        s1_t[i][j] = 0.0;
                        s2_t[i][j] = 0.0;
                }
        }


        if(allow_zero_counts){
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
        }else{
                for(iter = 0;iter < iterations;iter++){
                        RUN(fill_fast_transitions_only_matrices(model,ft));
                        for(i = 0;i < model->num_states;i++){
                                for(c = 0; c < model->L;c++){
                                        s1_e[i][c] += ft->emission[c][i];
                                        s2_e[i][c] += (ft->emission[c][i] * ft->emission[c][i]);
                                }
                        }
                        for(i = 0;i < model->num_states;i++){
                                for(c = 0;c < model->num_states;c++){
                                        if(model->transition_counts[i][c]){
                                                s1_t[i][c] += ft->transition[i][c];
                                                s2_t[i][c] += (ft->transition[i][c] * ft->transition[i][c]);
                                        }
                                }
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
        LOG_MSG("Adding fhmm to %s",h5file);

        //RUN(add_fhmm(h5file,s1_t,s1_e, model->num_states, model->L  ));

        gfree(s1_e);
        gfree(s1_t);

        gfree(s2_e);
        gfree(s2_t);

        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        return OK;
ERROR:
        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        return FAIL;
}

/*
int fill_background_emission_from_model(struct fast_hmm_param*ft, struct ihmm_model* model)
{
        int i,j;
        double sum = 0.0;
        ASSERT(ft != NULL, "No parameters");
        ASSERT(model != NULL, "No model ");
        for(i = 0; i < ft->L;i++){
                ft->background_emission[i] = 0.0;
        }

        for(i = 0; i < model->L;i++){
                for(j = 0 ; j < model->num_states;j++){
                        //fprintf(stdout,"%f ", model->emission_counts[i][j]);
                        ft->background_emission[i] += model->emission_counts[i][j];
                        sum +=  model->emission_counts[i][j];
                }
                //fprintf(stdout,"\n");
        }

        ASSERT(sum != 0.0,"No sequence counts found:%f ", sum);
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
        double sum = 0.0;

        ASSERT(ft != NULL, "No parameters");
        ASSERT(sb != NULL, "No sequences");

        for(i = 0; i < ft->L;i++){
                ft->background_emission[i] = 0.0;
        }

        for(i = 0; i < sb->num_seq;i++){
                for(j = 0;j < sb->sequences[i]->seq_len;j++){
                        ft->background_emission[sb->sequences[i]->seq[j]]++;
                        sum++;
                }
        }
        ASSERT(sum != 0.0,"No sequence counts found");
        for(i = 0; i < ft->L;i++){
                ft->background_emission[i] /= sum;
                //fprintf(stdout,"BACK: %f\n", ft->background_emission[i]);
        }
        return OK;
ERROR:
        return FAIL;
}*/
