#include "tldevel.h"

#include "randomkit_tl_add.h"
#include "model_struct.h"

#define MODEL_HELP_IMPORT
#include "model_help.h"

static int compare_model(struct ihmm_model* a, struct ihmm_model* b);

int print_counts(struct ihmm_model* ihmm)
{
        int i,j;
        int sum;
        fprintf(stdout,"States:%d\n",ihmm->num_states);

        fprintf(stdout,"Transition counts\n");
        for(i = 0; i < ihmm->num_states;i++){
                fprintf(stdout,"s%3d", i );
                for(j = 0; j < ihmm->num_states;j++){
                         fprintf(stdout," %0.2f", ihmm->transition_counts[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        fprintf(stdout,"Emission counts\n");
        sum = 0;
        for(i = 0; i < ihmm->L;i++){
                fprintf(stdout,"s%3d", i );
                for(j = 0; j < ihmm->num_states;j++){
                         fprintf(stdout," %0.5f", ihmm->emission_counts[i][j]);
                         sum+= ihmm->emission_counts[i][j];
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        fprintf(stdout,"SUM: %d\n", sum);
        return OK;
}

int print_model_parameters(struct ihmm_model* ihmm)
{
        int i;
        double sum = 0.0;
        fprintf(stdout,"NUMBER of STATES: %d\n", ihmm->num_states);
        fprintf(stdout,"%3.3f alpha\t%3.3f beta",ihmm->alpha , ihmm->gamma);
        for(i = 0; i < ihmm->num_states;i++){
                fprintf(stdout," %3.3f",ihmm->beta[i]);
                sum += ihmm->beta[i];
        }

        fprintf(stdout,"\tsum: %3.3f\n",sum);
        return OK;
}


int compare_model_bag(struct model_bag* a, struct model_bag* b)
{
        int i;
        ASSERT(a != NULL, "No model in a");
        ASSERT(b != NULL, "No model in b");
        if(a->best_model != b->best_model){
                WARNING_MSG("Models differ: best_model\t%d\t%d", a->best_model,b->best_model);
        }

        if(a->max_num_states != b->max_num_states){
                WARNING_MSG("Models differ: max_num_states\t%d\t%d", a->max_num_states,b->max_num_states);
        }

        if(a->num_models != b->num_models){
                WARNING_MSG("Models differ: num_models\t%d\t%d", a->num_models,b->num_models);
        }

        if(a->seed != b->seed){
                  WARNING_MSG("Models differ: seed\t%lud\t%lud", a->seed,b->seed);

        }
        for(i = 0; i < a->num_models;i++){
                if(a->min_u[i] != b->min_u[i]){
                        /* This should be ok. */
                        //WARNING_MSG("Models differ: min_u %d: \t%f\t%f",i, a->min_u[i],b->min_u[i]);
                }
        }
        for(i = 0; i < a->num_models;i++){
                RUN(compare_model(a->models[i],b->models[i]));
        }
        for( i = 0; i < a->num_models;i++){
                compare_rk_state(&a->models[i]->rndstate,&b->models[i]->rndstate);
        }


        return OK;
ERROR:
        return FAIL;
}

int compare_model(struct ihmm_model* a, struct ihmm_model* b)
{
        int i,j;
        ASSERT(a != NULL, "No model in a");
        ASSERT(b != NULL, "No model in b");

        if(a->alpha != b->alpha){
                WARNING_MSG("Models differ: alpha\t%f\t%f", a->alpha,b->alpha);
        }
        if(a->alpha_a != b->alpha_a){
                WARNING_MSG("Models differ: alpha_a\t%f\t%f", a->alpha_a,b->alpha_a);
        }
        if(a->alpha_b != b->alpha_b){
                WARNING_MSG("Models differ: alpha_b\t%f\t%f", a->alpha_b,b->alpha_b);
        }
        if(a->alpha_limit != b->alpha_limit){
                WARNING_MSG("Models differ: alpha_limit\t%f\t%f", a->alpha_limit,b->alpha_limit);
        }

        if(a->gamma != b->gamma){
                WARNING_MSG("Models differ: gamma\t%f\t%f", a->gamma,b->gamma);
        }
        if(a->gamma_a != b->gamma_a){
                WARNING_MSG("Models differ: gamma_a\t%f\t%f", a->gamma_a,b->gamma_a);
        }
        if(a->gamma_b != b->gamma_b){
                WARNING_MSG("Models differ: gamma_b\t%f\t%f", a->gamma_b,b->gamma_b);
        }
        if(a->gamma_limit != b->gamma_limit){
                WARNING_MSG("Models differ: gamma_limit\t%f\t%f", a->gamma_limit,b->gamma_limit);
        }

        if(a->num_states != b->num_states){
                WARNING_MSG("Models differ: num_states\t%f\t%f", a->num_states,b->num_states);
        }else{
                for(i = 0; i < a->num_states;i++){
                        if(a->beta[i] != b->beta[i]){
                                WARNING_MSG("Models differ: beta %d \t%f\t%f",i, a->beta[i],b->beta[i]);
                        }
                }
                for(i = 0; i < a->num_states;i++){
                        for(j = 0; j < a->num_states;j++){
                                if(a->transition_counts[i][j] != b->transition_counts[i][j]){
                                        WARNING_MSG("Models differ: trans counts \t%f\t%f",a->transition_counts[i][j] != b->transition_counts[i][j]);
                                }
                        }
                }

        }

        if(a->seed != b->seed){
                WARNING_MSG("Models differ: seed\t%f\t%f", a->seed,b->seed);
        }
        if(a->alloc_num_states != b->alloc_num_states){
                WARNING_MSG("Models differ: allocnum_states\t%f\t%f", a->alloc_num_states,b->alloc_num_states);
        }



        return OK;
ERROR:
        return FAIL;
}
