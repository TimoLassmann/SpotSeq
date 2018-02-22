
#include "beam_sample.h"

#include "fast_hmm_param_test_functions.h"


#ifdef ITESTBEAM
static int add_state_integration_test(void);
static int shrink_grow_integration_test(void);

int main(const int argc,const char * argv[])
{
        RUN(print_program_header((char * const*)argv,"Integration Test"));

        RUN(add_state_integration_test());

        RUN(shrink_grow_integration_test());
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


static int shrink_grow_integration_test(void)
{
       struct fast_hmm_param* ft = NULL;
       struct ihmm_model* model = NULL;
       struct ihmm_sequences* iseq = NULL;

       char *tmp_seq[4] = {
               "ACGT",
               "ACGT",
               "ACGT",
               "ACGT"};

       int initial_states = 4;
       int i;
       
       RUNP(ft = alloc_fast_hmm_param(initial_states,initial_states));

       RUNP(model = alloc_ihmm_model(initial_states, 4));
       RUNP(iseq = create_ihmm_sequences(tmp_seq ,4));
       RUN(random_label_ihmm_sequences(iseq, initial_states * 10));
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
       
       RUN(random_label_ihmm_sequences(iseq, 2));
       RUN(fill_counts(model,iseq));
       RUN(print_counts(model));
       
       RUN(iHmmHyperSample(model, 10));
       
       RUN(print_model_parameters(model));
       RUN(fill_fast_transitions(model,ft));
       RUN(print_fast_hmm_params(ft));

       for(i = 0;i < 10;i++){
               RUN(add_state_from_fast_hmm_param(model,ft));
       }
       RUN(print_model_parameters(model));
       RUN(print_fast_hmm_params(ft));

       RUN(random_label_ihmm_sequences(iseq, 6));
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
        struct ihmm_sequences* iseq = NULL;

        char *tmp_seq[4] = {
                "ACGT",
                "ACGT",
                "ACGT",
                "ACGT"};

        
        int initial_states = 4;
       
       

        /* Let's start allocating all structures  */
        RUNP(ft = alloc_fast_hmm_param(initial_states,initial_states));

        RUNP(model = alloc_ihmm_model(initial_states, 4));
        RUNP(iseq = create_ihmm_sequences(tmp_seq ,4));
        RUN(random_label_ihmm_sequences(iseq, initial_states));
      

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
        print_model_parameters(model); 
        RUN(fill_fast_transitions(model,ft));
        RUN(print_fast_hmm_params(ft));
        print_model_parameters(model); 
        RUN(add_state_from_fast_hmm_param(model,ft));
        RUN(add_state_from_fast_hmm_param(model,ft));
        RUN(print_fast_hmm_params(ft));
        print_model_parameters(model);       
        
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

int fill_fast_transitions(struct ihmm_model* model,struct fast_hmm_param* ft)
{
        struct fast_t_item** list = NULL;
        int i,j;
        int list_index;
        int last_index; 
        int last_state; 
        float sum;
        ASSERT(model != NULL, "No model");
        ASSERT(ft != NULL,"No fast_hmm_param structure");

        last_state = model->num_states -1;
        fprintf(stdout,"%d last state",last_state);
        /* check if there is enough space to hold new transitions... */
        /* This is slightly to generous as I am allocating memory for the
         * infinity state as well */
        RUN(expand_fast_hmm_param_if_necessary(ft, model->num_states *model->num_states  ));
        /* Empty previous transitions by setting the index to zero  */
        ft->num_items = 0;
        list_index = ft->num_items;

        
        list = ft->list;
        list_index = ft->num_items;
        

        /* Let's start with the transitions from the start state!  */
        last_index = list_index;
        sum = 0.0;
        
        /* Disallow Start to start transitions */
        list[list_index]->from = IHMM_START_STATE;
        list[list_index]->to   = IHMM_START_STATE;
        list[list_index]->t    = 0.0f;
        list_index++;
        
        /* Disallow Start to end transitions i.e. zero length sequences are not allowed*/
        list[list_index]->from = IHMM_START_STATE;
        list[list_index]->to   = IHMM_END_STATE; 
        list[list_index]->t    = 0.0f;
        list_index++;

        /* Now to the remaining existing transitions... */
        for(i = 2; i < last_state;i++){
                list[list_index]->from = IHMM_START_STATE;
                list[list_index]->to   = i;
                list[list_index]->t    = rk_gamma(&model->rndstate, model->transition_counts[IHMM_START_STATE][i] + model->beta[i] * model->alpha,1.0);
                sum += list[list_index]->t;
                list_index++;
        }
        /* the last to infinity transition (this is just used in stick breaking
         * when adding states ). here there should be no counts as this
         * possibility was not observed in the last transition. */
        list[list_index]->from = IHMM_START_STATE;
        list[list_index]->to   = last_state;
        list[list_index]->t    = rk_gamma(&model->rndstate, model->beta[last_state] * model->alpha,1.0);
        sum += list[list_index]->t;
        list_index++;

        /* Normalize!  */
        for(i = last_index; i < list_index;i++){
                list[i]->t = list[i]->t / sum; 
        }

        /* And now the stop state...  */
        /* There is no possibility escape the end state - all transitions from
         * end are zero. I am not sure if this matters in my dyn prig. code but
         * why not! */
        for(i = 0; i < last_state;i++){
                list[list_index]->from = IHMM_END_STATE;
                list[list_index]->to   = i;
                list[list_index]->t    = 0.0f;
                sum += list[list_index]->t;
                list_index++;
        }
        
        for(i = 2; i < last_state;i++){
                /* Remeber where I started filling...  */
                last_index = list_index; 
                sum = 0.0;
                
                for(j = 1; j < last_state;j++){
                        list[list_index]->from = i;
                        list[list_index]->to   = j;
                        list[list_index]->t    = rk_gamma(&model->rndstate, model->transition_counts[i][j] + model->beta[j] * model->alpha,1.0);
                        sum += list[list_index]->t;
                        list_index++;
                }
                list[list_index]->from = i;
                list[list_index]->to   = last_state;
                list[list_index]->t    = rk_gamma(&model->rndstate, model->beta[last_state] * model->alpha,1.0);
                sum += list[list_index]->t;
                list_index++;
                /* Normalize  */
                for(j = last_index; j < list_index;j++){
                        list[j]->t = list[j]->t / sum; 
                }
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
                        ft->emission[i][j] = rk_gamma(&model->rndstate, model->emission_counts[i][j]+ EMISSION_H, 1.0);
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
        ft->num_items = list_index;
        ft->last_state = last_state;

        
        
        return OK;
ERROR:
        return FAIL;
}

/* This function assumes (oh no!) that beta has space for an additional
g * element */
int add_state_from_fast_hmm_param(struct ihmm_model* ihmm,struct fast_hmm_param* ft)
{
        struct fast_t_item** list = NULL;
        float* beta;
        float alpha;
        float gamma;
        rk_state rndstate;
        
        float sum,be,bg,pe,pg, a,b;
        int i,new_k,list_index;
        int l,r;
        
        //int pg_hack;            /* I don't want add states that are not reachable. */
        //float* tmp_pg = NULL;


        ASSERT(ihmm != NULL, "No model");
        ASSERT(ft != NULL, "No ft.");

        
        rndstate = ihmm->rndstate;
        
        
        list_index = ft->num_items;
        
        /* First add empty space to host the newstate -> old state transitions. */
        if(list_index + ft->last_state + ft->last_state + 1 >= ft->alloc_num_states){
                RUN(expand_fast_hmm_param_if_necessary(ft, list_index + ft->last_state + ft->last_state + 1));
        }
        /* Check if model needs to be extended (mainly beta of course) */

        RUN(resize_ihmm_model(ihmm, ihmm->num_states + 1));

        ihmm->num_states = ihmm->num_states + 1;
        
        beta = ihmm->beta;
        alpha = ihmm->alpha;
        gamma = ihmm->gamma;

        
        new_k = ft->last_state;
        //fprintf(stdout,"LAST: %d\n",new_k);
        
        list = ft->list;
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
        }
        for(i = ft->num_items;i < list_index;i++){
                list[i]->t /= sum;
        }
        ft->num_items = list_index;
    
        
        
        //first get beta for new column
        be = beta[ft->last_state];
        bg = rk_beta(&rndstate, 1.0,gamma );
	
        beta[ft->last_state] = bg*be;
        beta[ft->last_state+1] = (1.0 - bg) *be;
        
        ihmm->beta = beta;
        //now split prob in last columns...
        a = alpha * beta[ft->last_state];
        b = 0.0;
        for(i = 0; i <= ft->last_state;i++){
                
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
        
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_asc);

        
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
  
        }


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
        ft->num_items = list_index;
        ft->last_state = new_k+1;
        ihmm->rndstate = rndstate;
        return OK;
ERROR:
        //if(tmp_pg){
        //        MFREE(tmp_pg);
        // }
        return FAIL;
}

