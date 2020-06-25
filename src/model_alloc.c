
#include "tldevel.h"
#include "randomkit.h"

#include "global.h"
#include "model_struct.h"
#include "null_model_emission.h"

#include "finite_hmm_alloc.h"

#define MODEL_ALLOC_IMPORT
#include "model_alloc.h"




struct model_bag* alloc_model_bag(int* num_state_array, int L, int num_models, int max_states, int seed)
{
        struct model_bag* b = NULL;

        int i;
        //unsigned long seed;
        unsigned long local_seed;
        ASSERT(num_models > 0,"need to allocate at least one model");


        //seed = rk_ulong(rndstate);
        MMALLOC(b, sizeof(struct model_bag));

        b->models = NULL;
        b->finite_models = NULL;
        b->min_u = NULL;
        b->num_models = num_models;
        b->max_num_states = max_states;
        b->best_model = -1;
        b->seed = seed;
        /* Set seed in model */
        if(seed){
                rk_seed(seed, &b->rndstate);
        }else{
                rk_randomseed(&b->rndstate);
        }
        MMALLOC(b->models, sizeof(struct ihmm_model*)* b->num_models);
        MMALLOC(b->min_u , sizeof(double) * b->num_models);
        for(i = 0; i < b->num_models;i++){
                /* set seed in each model based on RNG in main model bag */
                local_seed = rk_ulong(&b->rndstate);
                b->models[i] = NULL;
                b->min_u[i] = 0.0;
                RUNP(b->models[i] = alloc_ihmm_model(num_state_array[i],b->max_num_states, L,local_seed));
        }
        return b;
ERROR:
        free_model_bag(b);
        return NULL;
}

void free_model_bag(struct model_bag* b)
{
        int i;
        if(b){
                if(b->num_models){
                        for(i = 0; i < b->num_models;i++){
                                if(b->models[i]){
                                        free_ihmm_model(b->models[i]);
                                }
                                if(b->finite_models[i]){
                                        free_fhmm(b->finite_models[i]);
                                }
                        }
                        MFREE(b->finite_models);
                        MFREE(b->models);
                }

                if(b->min_u){
                        MFREE(b->min_u);
                }
                MFREE(b);
        }
}


struct ihmm_model* alloc_ihmm_model(int K, int maxK, int L, unsigned int seed)
{
        struct ihmm_model* model = NULL;
        int i,j;
        ASSERT(K>3, "No states requested");
        ASSERT(L>1, "No letters");

        ASSERT(K <= maxK,"Too many states");

        MMALLOC(model, sizeof(struct ihmm_model));


        model->transition_counts = NULL;
        model->emission_counts = NULL;
        model->beta = NULL;
        model->background = NULL;
        model->num_states = 0;
        model->alloc_num_states = maxK;
        model->training_iterations = 0;
        model->L = L;
        model->alpha = IHMM_PARAM_PLACEHOLDER;
        model->alpha_a = IHMM_PARAM_PLACEHOLDER;
        model->alpha_b = IHMM_PARAM_PLACEHOLDER;
        model->gamma = IHMM_PARAM_PLACEHOLDER;
        model->gamma_a = IHMM_PARAM_PLACEHOLDER;
        model->gamma_b = IHMM_PARAM_PLACEHOLDER;
        model->alpha_limit = DBL_MAX;
        model->gamma_limit = DBL_MAX;

        model->seed = seed;
        if(seed){
                rk_seed(seed, &model->rndstate);
        }else{
                rk_randomseed(&model->rndstate);
        }

        //while(K > model->alloc_num_states){
        //model->alloc_num_states = model->alloc_num_states << 1;
        //}
        model->num_states = K;

        /* RUNP(model->transition_counts = galloc(model->transition_counts, model->alloc_num_states, model->alloc_num_states, 0.0)); */
        /* RUNP(model->emission_counts = galloc(model->emission_counts , model->L, model->alloc_num_states, 0.0)); */
        model->background = NULL;

        RUN(get_null_model_emissions(&model->background, model->L));


        RUN(galloc(&model->transition_counts, model->alloc_num_states, model->alloc_num_states));
        for(i = 0; i < model->alloc_num_states;i++){
                for(j = 0; j < model->alloc_num_states;j++){
                        model->transition_counts[i][j] = 0.0;
                }
        }

        RUN(galloc(&model->emission_counts , model->L, model->alloc_num_states));
        for(i = 0; i < model->L;i++){
                for(j = 0; j < model->alloc_num_states;j++){
                        model->emission_counts[i][j] = 0.0;
                }

        }

        //MMALLOC(model->beta,sizeof(double) * model->alloc_num_states);
        RUN(galloc(&model->beta, model->alloc_num_states));
        for(i = 0; i < model->alloc_num_states;i++){
                model->beta[i] = 0.0;
        }

        return model;
ERROR:
        free_ihmm_model(model);
        return NULL;
}

int resize_ihmm_model(struct ihmm_model* ihmm, int K)
{
        int old_size;
        int i,j,c;
        ASSERT(ihmm != NULL, "No model");
        ASSERT(K > 3,"No states requested");

        old_size = ihmm->alloc_num_states;
        if(K >ihmm->alloc_num_states ){
                while(K > ihmm->alloc_num_states){
                        ihmm->alloc_num_states = ihmm->alloc_num_states << 1;
                }
                //LOG_MSG("Resizing model to %d states",ihmm->alloc_num_states);
                //RUNP(ihmm->transition_counts = galloc(ihmm->transition_counts, ihmm->alloc_num_states, ihmm->alloc_num_states, 0.0));
                //RUNP(ihmm->emission_counts = galloc(ihmm->emission_counts , ihmm->L, ihmm->alloc_num_states, 0.0));


                RUN(galloc(&ihmm->transition_counts, ihmm->alloc_num_states, ihmm->alloc_num_states));
                RUN(galloc(&ihmm->emission_counts , ihmm->L, ihmm->alloc_num_states));
                for(j = 0;j < ihmm->alloc_num_states;j++){
                        for(c = 0;c < ihmm->alloc_num_states;c++){
                                ihmm->transition_counts[j][c] = 0.0;
                        }
                }
                for(j = 0;j < ihmm->L;j++){
                        for(c = 0;c < ihmm->alloc_num_states;c++){
                                ihmm->emission_counts[j][c] = 0.0;
                        }
                }

                /* BUG - do I need to reset counts here? */

                RUN(galloc(&ihmm->beta, ihmm->alloc_num_states));
                //MREALLOC(ihmm->beta,sizeof(double) * ihmm->alloc_num_states);
                for(i = old_size;i < ihmm->alloc_num_states;i++){
                        ihmm->beta[i] = -1.0;
                }
        }
        return OK;
ERROR:
        free_ihmm_model(ihmm);
        return FAIL;
}

void free_ihmm_model(struct ihmm_model* ihmm)
{
        if(ihmm){
                gfree(ihmm->background);
                gfree(ihmm->transition_counts);
                gfree(ihmm->emission_counts);
                gfree(ihmm->beta);
                MFREE(ihmm);
        }
}
