
#include "beam_sample.h"

#include "fast_hmm_param_test_functions.h"

struct beam_thread_data{
        struct fast_hmm_param* ft;
        struct seq_buffer* sb;
        float** dyn;
        int thread_ID;
        int num_threads;
};


static struct beam_thread_data** create_beam_thread_data(int* num_threads, int max_len, int K);
static int resize_beam_thread_data(struct beam_thread_data** td,int* num_threads, int max_len, int K);
static void free_beam_thread_data(struct beam_thread_data** td, int num_threads);

void* do_dynamic_programming(void *threadarg);

//static int set_u(struct seq_buffer* sb, struct ihmm_model* model, float* min_u);
static int set_u(struct seq_buffer* sb, struct ihmm_model* model, struct fast_hmm_param* ft, float* min_u);
static int unset_u(struct seq_buffer* sb);


static int get_max_to_last_state_transition(struct fast_hmm_param*ft,float* max);
//static int check_if_ft_is_indexable(struct fast_hmm_param* ft, int num_states);


static int dynamic_programming(float** matrix,struct fast_hmm_param* ft, struct ihmm_sequence* ihmm_seq);

#ifdef ITESTBEAM
/* These are test funtions. */
        
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
        
           
        /* First initialize beam param struct.  */

        RUNP(sb = create_ihmm_sequences_mem(tmp_seq ,numseq));

        RUNP(model = alloc_ihmm_model(initial_states, sb->L));
        rk_randomseed(&model->rndstate);
        rk_randomseed(&model->rndstate);
        rk_randomseed(&model->rndstate);
        rk_randomseed(&model->rndstate);
        /* Initial guess... */
        model->alpha0_a = 6.0f;
        model->alpha0_b = 15.0f;
        model->gamma_a = 16.0f;
        model->gamma_b = 4.0f;
        model->alpha = IHMM_PARAM_PLACEHOLDER;
        model->gamma = IHMM_PARAM_PLACEHOLDER;

        
        RUNP(ft = alloc_fast_hmm_param(initial_states,sb->L));

        RUN(fill_background_emission(ft, sb));
        RUN(inititalize_model(model,sb,0));
        RUN(run_beam_sampling( model, sb, ft,NULL, 100, 10));

        
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
        
           
        /* First initialize beam param struct.  */

        RUNP(sb = create_ihmm_sequences_mem(tmp_seq ,numseq));

        RUNP(model = alloc_ihmm_model(initial_states, sb->L));
        /* Initial guess... */
        model->alpha0_a = 6.0f;
        model->alpha0_b = 15.0f;
        model->gamma_a = 16.0f;
        model->gamma_b = 4.0f;
        model->alpha = IHMM_PARAM_PLACEHOLDER;
        model->gamma = IHMM_PARAM_PLACEHOLDER;

        
        RUNP(ft = alloc_fast_hmm_param(initial_states,sb->L));

        RUN(fill_background_emission(ft, sb));
        
        RUN(inititalize_model(model,sb,0));
        RUN(run_beam_sampling( model, sb, ft,NULL, 100, 10));

        
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
       for(i = 0; i< 10;i++){
               fprintf(stdout,"%d->%d: %f\n", ft->list[i]->from,ft->list[i]->to,ft->list[i]->t);
       }

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
        //fprintf(stdout,"%d last state\n",last_state);
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
        //int l,r;
        
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
        
        for(i = 0;i < 10;i++){
                RUN(iHmmHyperSample(model, 20));
                //fprintf(stdout,"%f %f\n",model->alpha,model->gamma );
                //model->alpha = 1.0;
                //model->gamma = 1.0;
        }
        
        /* sample transitions / emission */
        RUN(fill_fast_transitions(model,ft));

        
        /* Threading setup...  */
        need_local_pool = 0;
        if(pool){
                local_pool = pool;
        }else{
                if((local_pool = thr_pool_create(num_threads,num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");
                need_local_pool =1;
        }
        
        RUNP(td = create_beam_thread_data(&num_threads,(sb->max_len+1)  , ft->last_state));
        LOG_MSG("Will use %d threads.", num_threads);

        
        
        for(iter = 0;iter < iterations;iter++){//}iterations;iter++){
                /* Set U */
                
                RUN(set_u(sb,model,ft, &min_u));
                        
                //fprintf(stdout,"MIN_U:%f\n",min_u);
                //print_fast_hmm_params(ft);
                RUN(get_max_to_last_state_transition(ft, &max));
                //fprintf(stdout,"MAX:%f\n", max);
                while(max > min_u && model->num_states < 300){//}sb->max_len){
                        //fprintf(stdout,"ITER: %d Add state! MAX:%f min_U:%f max_len: %d \n",iter , max, min_u,sb->max_len);
                        RUN(add_state_from_fast_hmm_param(model,ft));
                        RUN(get_max_to_last_state_transition(ft, &max));
                        //fprintf(stdout,"MAX:%f min_U:%f\n", max, min_u);
                        
                }

        
                RUN(make_flat_param_list(ft));
                LOG_MSG("Iteration %d (%d states)  alpha = %f, gamma = %f", iter, model->num_states, model->alpha ,model->gamma);
                //print_fast_hmm_params(ft);
                //RUNP(matrix = malloc_2d_float(matrix,sb->max_len+1, ft->last_state, 0.0f));
                RUN(resize_beam_thread_data(td, &num_threads,(sb->max_len+1)  , ft->last_state));
                
                
                //dyn prog + labelling
                //LOG_MSG(" %d * %d = %d ",sb->max_len+1, ft->last_state, sizeof(float)*(sb->max_len+1) * ft->last_state );
                //qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*), fast_hmm_param_cmp_by_t_desc);
                for(i = 0; i < num_threads;i++){
                        td[i]->ft = ft;


                        
                        td[i]->sb = sb;
                        if(thr_pool_queue(local_pool,do_dynamic_programming,td[i]) == -1){
                                fprintf(stderr,"Adding job to queue failed.");
                        }
                }
                
                //print_labelled_ihmm_seq(sb->sequences[0]);
                //for(i = 0; i < sb->num_seq;i++){
                //        RUN(dynamic_programming(matrix,ft, sb->sequences[i]));
                        
                //}
                thr_pool_wait(local_pool);

                no_path =0;
                for(i = 0; i < sb->num_seq;i++){
                        if(sb->sequences[i]->u[0] == -1){
                                no_path = 1;
                                LOG_MSG("weird split must have happened in seq %d",i);
                        }
                }
          
                //print_labelled_ihmm_seq(sb->sequences[0]);
                //exit(0);
                if(no_path){
                        LOG_MSG("weird split must have happened. %d",iter);
                        //print_fast_hmm_params(ft);
                        //RUN(print_counts(model));
                        
                        RUN(remove_unused_states_labels(model, sb));
                        RUN(fill_counts(model,sb));
                             
                        
                        //RUN(print_counts(model));
                        RUN(fill_fast_transitions(model,ft));
                        /*for(i = 0; i < sb->num_seq;i++){
                          print_labelled_ihmm_seq(sb->sequences[i]);
                          }
                          exit(0);*/
                        iterations++;
                }else{
                        //res = fast_hmm_param_binarySearch_t(ft, x);
                        //RUN(fill_counts(model,sb));
                        //RUN(print_counts(model));
        
                        /* I am doing this as a pre-caution. I don't want the inital model
                         * contain states that are not visited.. */
                        //RUN(remove_unused_states_labels(model, sb));
                        //LOG_MSG("%d states", model->num_states);
             
                        //print_labelled_ihmm_seq(sb->sequences[0]);
                        //remove unwantrd.
                        
                        RUN(remove_unused_states_labels(model, sb));
                        RUN(fill_counts(model,sb));
                       
                        //RUN(print_counts(model));
                        //hyper
                        //RUN(iHmmHyperSample(model, 20));
                        //fprintf(stdout,"%f %f\n",model->alpha,model->gamma );
                        //if(iter < 100){
                                /*     model->alpha = 1.0;
                                model->gamma = 1.0;
                                for(i = 0; i < model->num_states;i++){
                                        model->beta[i] = 2.0;//1.0 / (float)(model->num_states);
                                }

                                float sum = 0.0;
                                model->beta[0] = 0;
                                for(i = 1; i < model->num_states-1;i++){
                                        model->beta[i] = rk_gamma(&model->rndstate, 1.0, 1.0);
                                        sum += model->beta[i];
                                }
	
                                model->beta[model->num_states-1] =  rk_gamma(&model->rndstate, model->gamma, 1.0);
                                sum += model->beta[model->num_states-1] ;
                                for(i = 0; i < model->num_states;i++){
                                        model->beta[i] /= sum;
                                        }*/
                                
                                

                        //}else{
                                RUN(iHmmHyperSample(model, 1));
                                //}
                        // fill fast...
                        RUN(fill_fast_transitions(model,ft));
                
                        //print_fast_hmm_params(ft);
                }
                /* print out model  */
                if((iter+1) % 10 == 0){
                        LOG_MSG("print %d\n",iter);
                        char tmp_buffer[BUFFER_LEN];
                        snprintf(tmp_buffer,BUFFER_LEN,"model_at_%07d.h5",iter+1);
                        RUN(write_model_hdf5(model,tmp_buffer));
                        //RUN(add_annotation(param->output, "spotseq_model_cmd", param->cmd_line));
                        RUN(add_background_emission(tmp_buffer,ft->background_emission,ft->L));
                        RUN(run_build_fhmm(tmp_buffer));
                }

        }

        RUN(print_labelled_ihmm_buffer(sb));
        
        if(need_local_pool){
                 thr_pool_destroy(local_pool);
        }

        free_beam_thread_data(td, num_threads);
        //free_2d((void**) matrix);
        
        return OK;
ERROR:
        free_beam_thread_data(td, num_threads);
       
        //if(matrix){
        //       free_2d((void**) matrix);
        //}
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
                //  RUNP(matrix = malloc_2d_float(matrix,sb->max_len+1, ft->last_state, 0.0f));
                RUNP(td[i]->dyn = malloc_2d_float(td[i]->dyn, max_len, K, 0.0f));
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

        
        for(i = 0; i < local_num_treads;i++){
                RUNP(td[i]->dyn = malloc_2d_float(td[i]->dyn, max_len, K, 0.0f));
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
                        MFREE(td[i]);
                }
        
                MFREE(td);
        }
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
                }
        }
        return NULL;
ERROR:
        return NULL;
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
        emission = ft->emission[0];
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
                u[0] = -1.0f;
                //LOG_MSG("No PATH!: %f",sum);
        }
        
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
                //fprintf(stdout,"%d->%d %f\n", ft->infinity[j]->from, ft->infinity[j]->to, ft->infinity[j]->t);
        }
        *max = local_max;
        return OK;
ERROR:
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
 

int run_build_fhmm(char* h5file)
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
