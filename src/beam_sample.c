
#include "beam_sample.h"

#include "fast_hmm_param_test_functions.h"



//static int set_u(struct seq_buffer* sb, struct ihmm_model* model, float* min_u);
static int set_u(struct seq_buffer* sb, struct ihmm_model* model, struct fast_hmm_param* ft, float* min_u);
static int get_max_to_last_state_transition(struct fast_hmm_param*ft,float* max);
static int check_if_ft_is_indexable(struct fast_hmm_param* ft, int num_states);


/* These are test funtions. */
        
static int add_state_integration_test(void);
static int shrink_grow_integration_test(void);
static int full_run_test(void);




int main(const int argc,const char * argv[])
{
        RUN(print_program_header((char * const*)argv,"Integration Test"));

        //UN(add_state_integration_test());

        //RUN(shrink_grow_integration_test());

        RUN(full_run_test());
        
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
                "AGTGGATATCACAGGCTAAAGGAGGGG"};
        int i;
        int numseq = 8;
        int initial_states = 8;
        
           
        /* First initialize beam param struct.  */

        RUNP(sb = create_ihmm_sequences_mem(tmp_seq ,numseq));

        RUNP(model = alloc_ihmm_model(initial_states, 4));
        /* Initial guess... */
        model->alpha0_a = 4.0f;
        model->alpha0_b = 2.0f;
        model->gamma_a = 3.0f;
        model->gamma_b = 6.0f;
        model->alpha = IHMM_PARAM_PLACEHOLDER;
        model->gamma = IHMM_PARAM_PLACEHOLDER;

        
        RUNP(ft = alloc_fast_hmm_param(initial_states,initial_states));

        RUN(run_beam_sampling( model, sb, ft, 100, 10));

        
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

//#endif

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
        //fprintf(stdout,"%d last state\n",last_state);
        /* check if there is enough space to hold new transitions... */
        /* This is slightly to generous as I am allocating memory for the
         * infinity state as well */

        RUN(expand_emission_if_necessary(ft, model->num_states));
        //RUN(expand_fast_hmm_param_if_necessary(ft, model->num_states *model->num_states  ));
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
        if(list_index == ft->alloc_items){
                RUN(expand_transition_if_necessary(ft));
                list = ft->list;
        }
        
        /* Disallow Start to end transitions i.e. zero length sequences are not allowed*/
        list[list_index]->from = IHMM_START_STATE;
        list[list_index]->to   = IHMM_END_STATE; 
        list[list_index]->t    = 0.0f;
        list_index++;
        if(list_index == ft->alloc_items){
                RUN(expand_transition_if_necessary(ft));
                list = ft->list;
        }

        /* Now to the remaining existing transitions... */
        for(i = 2; i < last_state;i++){
                list[list_index]->from = IHMM_START_STATE;
                list[list_index]->to   = i;
                list[list_index]->t    = rk_gamma(&model->rndstate, model->transition_counts[IHMM_START_STATE][i] + model->beta[i] * model->alpha,1.0);
                sum += list[list_index]->t;
                list_index++;
                if(list_index == ft->alloc_items){
                        RUN(expand_transition_if_necessary(ft));
                        list = ft->list;
                }
        }
        /* the last to infinity transition (this is just used in stick breaking
         * when adding states ). here there should be no counts as this
         * possibility was not observed in the last transition. */
        list[list_index]->from = IHMM_START_STATE;
        list[list_index]->to   = last_state;
        list[list_index]->t    = rk_gamma(&model->rndstate, model->beta[last_state] * model->alpha,1.0);
        sum += list[list_index]->t;
        list_index++;
        if(list_index == ft->alloc_items){
                RUN(expand_transition_if_necessary(ft));
        }

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
                if(list_index == ft->alloc_items){
                        RUN(expand_transition_if_necessary(ft));
                        list = ft->list;
                }
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

                        if(list_index == ft->alloc_items){
                                RUN(expand_transition_if_necessary(ft));
                                list = ft->list;
                        }
                }
                list[list_index]->from = i;
                list[list_index]->to   = last_state;
                list[list_index]->t    = rk_gamma(&model->rndstate, model->beta[last_state] * model->alpha,1.0);
                sum += list[list_index]->t;
                list_index++;
                if(list_index == ft->alloc_items){
                        RUN(expand_transition_if_necessary(ft));
                        list = ft->list;
                }
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
        sum = 0.0f;
        for(j = 0; j < list_index;j++){
                if(ft->list[j]->t == 0.0f){
                        sum += 1.0;
                }
        }
        //LOG_MSG("TOTAL: %d trans %f elements are blank!", list_index,  sum);
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
        /* Sorting is only strictly necessary if this is called after another function re-sorted it */
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_from_asc);
         
        rndstate = ihmm->rndstate;
        
        
        list_index = ft->num_items;
        
        /* First add empty space to host the newstate -> old state transitions. */
        //if(list_index + ft->last_state + ft->last_state + 1 >= ft->alloc_num_states){
        //        LOG_MSpG("requesting more memory in add state...");
                //RUN(expand_fast_hmm_param_if_necessary(ft, list_index + ft->last_state + ft->last_state + 1));
        //}
        /* Check if model needs to be extended (mainly beta of course) */

        RUN(resize_ihmm_model(ihmm, ihmm->num_states + 1));

        ihmm->num_states = ihmm->num_states + 1;
        RUN(expand_emission_if_necessary(ft, ihmm->num_states));
        
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
                if(list_index == ft->alloc_items){
                        RUN(expand_transition_if_necessary(ft));
                        list = ft->list;
                }
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


int run_beam_sampling(struct ihmm_model* model, struct seq_buffer* sb, struct fast_hmm_param* ft, int iterations, int num_threads)
{
        int i,j;
        int initial_number_of_states;
        int K;
        int iter;
        float min_u;
        float max;
        float** matrix = NULL;
        ASSERT(model != NULL, "nop model.");
        ASSERT(sb,"no sequence buffer");
        ASSERT(sb->num_seq > 0, "No sequences");
        ASSERT(ft != NULL, "No transition struct");
        ASSERT(iterations > 1, "No iterations");
        ASSERT(num_threads > 0, "No threads");

        /* First guess initial number of states. For now I simply pick K = average seq len  */
        
        K = 0;
        for(i = 0; i < sb->num_seq;i++){
                K += sb->sequences[i]->seq_len;
        }
        
        K = K / sb->num_seq;
        K = 100;
        
        LOG_MSG("Will start with %d states",K);
        
        RUN(random_label_ihmm_sequences(sb, K));
        RUN(fill_counts(model,sb));
        /* I am doing this as a pre-caution. I don't want the inital model
         * contain states that are not visited.. */
        RUN(remove_unused_states_labels(model, sb));
        RUN(fill_counts(model,sb));
        RUN(print_counts(model));
        
        for(i = 0;i < 10;i++){
                RUN(iHmmHyperSample(model, 20));
                
                RUN(print_model_parameters(model));
                for(j = 0; j < model->num_states;j++){
                        model->beta[j] = 1.0 / (float)(model->num_states);
                }
        }
        
        /* sample transitions / emission */
        RUN(fill_fast_transitions(model,ft));
        print_fast_hmm_params(ft);
        
        /* super important to make sure transitions are index-able!! */
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_from_asc);
        
        //just to make sure!
        RUN(check_if_ft_is_indexable(ft,ft->last_state));
        
        for(iter = 0;iter < 1000;iter++){//}iterations;iter++){
                /* Set U */
                //LOG_MSG("Iteration %d", iter);
                RUN(set_u(sb,model,ft, &min_u));
                //fprintf(stdout,"MIN_U:%f\n",min_u);
                //print_fast_hmm_params(ft);
                RUN(get_max_to_last_state_transition(ft, &max));
                //fprintf(stdout,"MAX:%f\n", max);
                while(max > min_u && model->num_states < sb->max_len){
                        //fprintf(stdout,"Add state! MAX:%f min_U:%f\n", max, min_u);
                        RUN(add_state_from_fast_hmm_param(model,ft));
                        RUN(get_max_to_last_state_transition(ft, &max));
                        //fprintf(stdout,"MAX:%f min_U:%f\n", max, min_u);
                        
                }
                //print_fast_hmm_params(ft);
                RUNP(matrix = malloc_2d_float(matrix,sb->max_len+1, ft->last_state, 0.0f));
                //dyn prog + labelling
                //LOG_MSG(" %d * %d = %d ",sb->max_len+1, ft->last_state, sizeof(float)*(sb->max_len+1) * ft->last_state );
                qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*), fast_hmm_param_cmp_by_t_desc);
                
               
                //print_labelled_ihmm_seq(sb->sequences[0]);
                for(i = 0; i < sb->num_seq;i++){
                        RUN(dynamic_programming(matrix,ft, sb->sequences[i]));
                }
                //print_labelled_ihmm_seq(sb->sequences[0]);
                //exit(0);
                
                //res = fast_hmm_param_binarySearch_t(ft, x);
                RUN(fill_counts(model,sb));
                //RUN(print_counts(model));
        
                /* I am doing this as a pre-caution. I don't want the inital model
                 * contain states that are not visited.. */
                RUN(remove_unused_states_labels(model, sb));
             
                //print_labelled_ihmm_seq(sb->sequences[0]);
                //remove unwantrd.
                RUN(fill_counts(model,sb));
                //hyper
                RUN(iHmmHyperSample(model, 20));
                // fill fast...
                RUN(fill_fast_transitions(model,ft));
                // print_fast_hmm_params(ft);
        }

        for(i = 0; i < sb->num_seq;i++){
                print_labelled_ihmm_seq(sb->sequences[i]);
        }
        free_2d((void**) matrix);
        
        return OK;
ERROR:
        if(matrix){
                free_2d((void**) matrix);
        }
        return FAIL;
}

int dynamic_programming(float** matrix,struct fast_hmm_param* ft, struct ihmm_sequence* ihmm_seq)
{
        int i,j,len,boundary;
        float* u = NULL;
        uint8_t* seq = NULL;
        int* label = NULL;
        int a,b;
        float sum;
        float* emission; 
        float* tmp_row;
        float r;
        int last_pick;
        struct fast_t_item** list = NULL;
        ASSERT(ft!= NULL, "no parameters");
        ASSERT(matrix != NULL,"No dyn matrix");

        u = ihmm_seq->u;
        len = ihmm_seq->seq_len;
        seq = ihmm_seq->seq;
        label = ihmm_seq->label;

        list = ft->list;
        emission = ft->emission;
        tmp_row = matrix[len];
        
        
        boundary = fast_hmm_param_binarySearch_t(ft, u[0]);
        for(i = 0; i < ft->last_state;i++){
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
        for(i = 0; i < ft->last_state;i++){
                
                matrix[0][i] *=  emission[i];

                sum += matrix[0][i];
        }
        for(i = 0; i < ft->last_state;i++){
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
                

                for(j = 0; j < ft->last_state;j++){
                        matrix[i][j] *=  emission[j];
                        sum += matrix[i][j];
                        
                }
                for(j = 0; j < ft->last_state;j++){
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
        last_pick = IHMM_END_STATE;
        /* First let's check if there is a path! i.e. end is reachable.  */

        sum = 0.0f;
        
        boundary = fast_hmm_param_binarySearch_t(ft, u[len]);
        for(j = 0; j < boundary;j++){
                a = list[j]->from;
                b = list[j]->to;
                if(list[j]->to == last_pick){
                        sum += matrix[len-1][a];
                }
        }
        
        if(sum != 0.0 && !isnan(sum)){
                last_pick = IHMM_END_STATE;
                for(i = len-1; i >= 0; i--){
                        //fprintf(stdout,"pick: %d %d\n",i,last_pick);
                        for(j = 0; j < ft->last_state;j++){
                                tmp_row[j] = -1.0;
                        }
                        sum = 0.0f;
                        boundary = fast_hmm_param_binarySearch_t(ft, u[i+1]);
                        for(j = 0; j < boundary;j++){
                                a = list[j]->from;
                                b = list[j]->to;
                                if(list[j]->to == last_pick){
                                        tmp_row[a] = matrix[i][a];
                                        sum += matrix[i][a];
                                }
                        }
                        r =  random_float_zero_to_x(sum);
                        for(j = 0; j < boundary;j++){
                                
                                a = list[j]->from;
                                b = list[j]->to;
                                if(list[j]->to == last_pick){
     
                                        r -= tmp_row[a];
                                        if(r <= 0.0f){
                                                last_pick = a;
                                        
                                                break;
                                        }
                                }
                        }
                        label[i] = last_pick;   
                }

        }else{
                LOG_MSG("No PATH!: %f",sum);
        }
        
        return OK;
ERROR:
        return FAIL;
}


int set_u(struct seq_buffer* sb, struct ihmm_model* model, struct fast_hmm_param* ft, float* min_u)
{
        int i,j,c,a,b;
        float* u = 0;
        int* label =0;
        int len;
        int last_state = 0;
        
        float local_min_u = 1.0;
        ASSERT(sb != NULL, "No sequences.");
        ASSERT(model != NULL, "No model.");
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_from_asc);
        last_state = ft->last_state;
        
        for(i = 0; i < sb->num_seq;i++){
                
                label = sb->sequences[i]->label;
                u = sb->sequences[i]->u;
                len = sb->sequences[i]->seq_len;
                c = IHMM_START_STATE * last_state + label[0];
                //c = a* (num_states-1) + b;
                u[0] =  rk_double(&model->rndstate) *(ft->list[c]->t);
                ASSERT(ft->list[c]->t != 0.0f,"BAD %d -> %d %f",ft->list[c]->from,ft->list[c]->to,ft->list[c]->t);
                
                local_min_u = MACRO_MIN(local_min_u, u[0]);
                for (j = 1; j < len;j++){
                        c = label[j-1] * last_state + label[j];
                        u[j] =  rk_double(&model->rndstate) *(ft->list[c]->t);//rk_double(&model->rndstate) *
                        //if(!i && j < 5){
                        //       fprintf(stdout,"%d->%d %f\n",label[j-1],label[j],ft->list[c]->t );
                        //}
                        //fprintf(stdout,"%d %d  ;; %d %d\n",label[j-1],label[j],ft->list[c]->from ,ft->list[c]->to);
                        local_min_u = MACRO_MIN(local_min_u, u[j]);
                        ASSERT(ft->list[c]->t != 0.0f,"BAD %d -> %d %f",ft->list[c]->from,ft->list[c]->to,ft->list[c]->t);

                        
                }

                c = label[len-1] * last_state + IHMM_END_STATE;
                u[len] =  rk_double(&model->rndstate) *(ft->list[c]->t);
                                 ASSERT(ft->list[c]->t != 0.0f,"BAD %d -> %d %f",ft->list[c]->from,ft->list[c]->to,ft->list[c]->t);
                                 //fprintf(stdout,"%d %d -> %d: %f  \n",label[len-1],ft->list[c]->from ,ft->list[c]->to, ft->list[c]->t );
                local_min_u = MACRO_MIN(local_min_u, u[len]);
                
        }
        
        *min_u = local_min_u;
        return OK;
ERROR:
        print_labelled_ihmm_seq(sb->sequences[i]);
        return FAIL;
}

int get_max_to_last_state_transition(struct fast_hmm_param*ft,float* max)
{

        int i,l,r;
        float local_max;
        struct fast_t_item** list  = NULL;
        ASSERT(ft != NULL, "No fast hmm parameters.");
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_asc);

        list = ft->list;
        
        l = fast_hmm_param_binarySearch_to_lower_bound(ft,ft->last_state);
        r = fast_hmm_param_binarySearch_to_upper_bound(ft,ft->last_state);
        local_max = -1.0f;
        for(i = l;i < r;i++){
                if(list[i]->t > local_max){
                        local_max = list[i]->t;
                }
        }
        *max = local_max;
        return OK;
ERROR:
        return FAIL;
}


/* Purpose is to check if can random access a particular transition.. */
int check_if_ft_is_indexable(struct fast_hmm_param* ft, int num_states)
{
       
        int i,a,b,c;
        ASSERT(ft != NULL, "No fast hmm parameters");       

        for(i = 0; i < 1000;i++){

                /* NOTE: below has to be -1 because a) the random function
                 * returns values including the last number. */


                /* NOTE: transitions from any state to the start state do not
                 * appear in FT; */
                /* NOTE: transitions from the end state do not appear in FT; */
                a = IHMM_END_STATE;
                b = IHMM_START_STATE;
                while(a == IHMM_END_STATE || b == IHMM_START_STATE){
                        a = random_int_zero_to_x(num_states-1);
                        b = random_int_zero_to_x(num_states-1);
                }
                
                c = a* (num_states) + b;
                //fprintf(stdout," %d -> %d ;; %d ->%d\n", a,b,ft->list[c]->from ,ft->list[c]->to);
                ASSERT(a == ft->list[c]->from,"ft is not index-able. Did you sort based on from then to?");
                ASSERT(b == ft->list[c]->to,"ft is not index-able. Did you sort based on from then to?");
        }
        return OK;
ERROR:
        return FAIL;
}
