
#include "model.h"
#include "ihmm_seq.h"


int inititalize_model(struct ihmm_model* model, struct seq_buffer* sb, int K)
{
        int i;
        if(K == 0){
                K = 0;
                for(i = 0; i < sb->num_seq;i++){
                        K += sb->sequences[i]->seq_len;
                }
        
                K = K / sb->num_seq;
        }
        LOG_MSG("Will start with %d states",K);
        //K = 10;
        RUN(random_label_ihmm_sequences(sb, K));
        RUN(fill_counts(model,sb));
        /* I am doing this as a pre-caution. I don't want the inital model
         * contain states that are not visited.. */
        RUN(remove_unused_states_labels(model, sb));
        RUN(fill_counts(model,sb));
        return OK;
ERROR:
        return FAIL;

}

int remove_unused_states_labels(struct ihmm_model* ihmm, struct seq_buffer* sb)
{
        int i,j;
        float sum;
        int len;
        int* relabel = NULL;
        int* used = NULL;
        int* lab = NULL;

        ASSERT(ihmm != NULL, "no model");
        ASSERT(sb != NULL, "no seq struct");

        MMALLOC(relabel, sizeof(int) * ihmm->num_states);
        MMALLOC(used, sizeof(int) * ihmm->num_states);
        for(i = 0; i < ihmm->num_states;i++){
                used[i] = 0;
                relabel[i] = -1;
        }
        used[IHMM_START_STATE] = 100;
        used[IHMM_END_STATE] = 100;
                
        for(i = 0; i < sb->num_seq;i++){
                lab = sb->sequences[i]->label;
                len = sb->sequences[i]->seq_len;
                for(j = 0; j < len;j++){
                        used[lab[j]]++;
                }
        }

        
        //fprintf(stdout,"ORG beta \n");
        //for(i = 0; i < ihmm->num_states;i++){
        //        fprintf(stdout,"%3.3f ", ihmm->beta[i]);
        //}
        //fprintf(stdout,"\n");
        
        j = 0;
        sum = 0.0;
        for(i = 0; i < ihmm->num_states;i++){
                if(used[i] != 0){
                        ihmm->beta[j] = ihmm->beta[i];
                        relabel[i] = j;
                        j++;
                }else{
                        relabel[i] = j;
                        sum += ihmm->beta[i];
                }
        }
        
        /*for(i = 0; i < ihmm->num_states;i++){
                fprintf(stdout,"%3d ",used[i]);
        }
        fprintf(stdout,"\n");
        for(i = 0; i < ihmm->num_states;i++){
                fprintf(stdout,"%3d ",relabel[i]);
        }
        fprintf(stdout,"\n");*/
        
        ihmm->beta[j] = sum;
        //for(i = j+1; i < ihmm->num_states;i++){
        //       ihmm->beta[i] = 0;
        //}
        ihmm->num_states = j+1; /* need to add one for the infinite stuff */
        RUN(resize_ihmm_model(ihmm, j+1));
        //fprintf(stdout,"CUR beta \n");
        sum = 0.0f;
        for(i = 0; i < ihmm->num_states;i++){
                //        fprintf(stdout,"%3.3f ", ihmm->beta[i]);
                sum+= ihmm->beta[i];
        }
        
        //fprintf(stdout,"\tsum: %f\n",sum);
        
        
        for(i = 0; i < sb->num_seq;i++){
                //       fprintf(stdout,"%3d",i);
                lab = sb->sequences[i]->label;
                len = sb->sequences[i]->seq_len;
                for(j= 0; j <  len;j++){
                        lab[j] = relabel[lab[j]];
                        //                fprintf(stdout," %d",lab[j]);
                }
                //       fprintf(stdout,"\n");
        }
        MFREE(used);
        MFREE(relabel);
        return OK;
ERROR:
        MFREE(used);
        MFREE(relabel);
        return FAIL;
}

int fill_counts(struct ihmm_model* ihmm, struct seq_buffer* sb)
{
        int i,j;
        int* label = NULL;
        uint8_t* seq;
        float** p = NULL;
        int max_state_ID;
        int len;
        ASSERT(ihmm != NULL,"No model.");
        ASSERT(sb != NULL,"No iseq struct");


        /* First I need to check what the largest state ID is and see if we have sufficient space allocated in the model.  */
        max_state_ID = -1;
        //fprintf(stdout,"%d numseq\n",sb->num_seq );
        for(i = 0; i < sb->num_seq;i++){
                label = sb->sequences[i]->label;
                len = sb->sequences[i]->seq_len;
                //fprintf(stdout,"%d len%d\n",i,len);
                for(j = 0; j < len;j++){
                        if(label[j] > max_state_ID){
                                max_state_ID = label[j];
                        }
                }
        }
        max_state_ID += 1; // for the infinity possibility; not observed in the current labeling
        max_state_ID += 1; // so I can use the < syntax rather than <= 
        ASSERT(max_state_ID > 2, "Not enough states found");

        RUN(resize_ihmm_model(ihmm, max_state_ID));
        ihmm->num_states = max_state_ID; 
        
        p = ihmm->transition_counts;
        /* clear transition counts */
        for(i = 0; i < ihmm->num_states ;i++){
                for(j = 0; j < ihmm->num_states ;j++){
                        p[i][j] = 0.0f;
                }
        }
        /* get transition counts */
        for(i = 0; i < sb->num_seq;i++){
                label = sb->sequences[i]->label;
                len = sb->sequences[i]->seq_len;
                
                 p[IHMM_START_STATE][label[0]] += 1;
                for(j = 1; j < len;j++){
                        p[label[j-1]][label[j]]++;
                }
		
                p[label[len-1]][IHMM_END_STATE] += 1;
        }

        
        /* clear emission counts */
        p = ihmm->emission_counts;
        for(i = 0; i < ihmm->L ;i++){
                for(j = 0; j < ihmm->num_states ;j++){
                        p[i][j] = 0.0;
                }
        }
        /* get emission counts */
        for(i = 0; i < sb->num_seq; i++){
                label = sb->sequences[i]->label;
                seq = sb->sequences[i]->seq;
                len = sb->sequences[i]->seq_len;
                for(j =0; j < len;j++){
                        p[(int)seq[j]][label[j]]++;
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int iHmmHyperSample(struct ihmm_model* model, int iterations)
{
        int i,j,c;
        int last_state;
        float** transition_counts = NULL;
        float** M = NULL;
        float** supp = NULL;
        float* sum_M = NULL;
        float* sum_N = NULL;
        float* w = NULL;
        float* p = NULL;
        float* s = NULL;
        float total_M = 0.0f;
        float alpha;
        float gamma;
        float sum, sum_s, sum_w, mu, pi_mu;

        ASSERT(model != NULL, "No model");
        ASSERT(iterations > 0, "No iterations");
        last_state = model->num_states-1;

        /* alloc auxillary data structures  */
        RUNP(M = malloc_2d_float(M, model->num_states,model->num_states, 0.0f));
        RUNP(supp = malloc_2d_float(supp, 5,model->num_states, 0.0f));
        
        /* am I in the first iteration? */
        if(model->alpha == IHMM_PARAM_PLACEHOLDER || model->gamma == IHMM_PARAM_PLACEHOLDER){
                /* What shall I do with an empty (just alloced) model with no
                 * _a/b variables ? Throw an error! */
                ASSERT(model->alpha0_a != IHMM_PARAM_PLACEHOLDER,"you need to set alpha0_a before calling iHmmHyperSample .");
                ASSERT(model->alpha0_b != IHMM_PARAM_PLACEHOLDER,"you need to set alpha0_b before calling iHmmHyperSample .");
                ASSERT(model->gamma_a != IHMM_PARAM_PLACEHOLDER,"you need to set gamma_a before calling iHmmHyperSample .");
                ASSERT(model->gamma_b != IHMM_PARAM_PLACEHOLDER,"you need to set gamma_b before calling iHmmHyperSample .");

                /* all good let's set the initial guess parameters...  */
                /* First of all initialize the random number generator */
                rk_randomseed(&model->rndstate);
                
                model->alpha = rk_gamma(&model->rndstate, model->alpha0_a,1.0 / model->alpha0_b);
                model->gamma = rk_gamma(&model->rndstate, model->gamma_a,1.0 / model->gamma_b);
                /* Note this also initializes the last (to infinity state) */
                for(i = 0; i < model->num_states;i++){
                        model->beta[i] = 1.0 / (float)(model->num_states);
                }
                /* Ok all set to re-estimate hyper parameters..  */
        }

        alpha = model->alpha;
        gamma = model->gamma;
        
        //WAS here - need to alloc auxillary arrays for calculations below!         
        transition_counts = model->transition_counts;
        sum_M = supp[0];
        sum_N = supp[1]; 
        
        total_M = 0.0f;
        for(i = 0; i < last_state;i++){
                for(j = 0; j < last_state;j++){
                        if(transition_counts[i][j] == 0){
                                M[i][j] = 0;
                        }else{
                                for(c = 1;c <= transition_counts[i][j];c++){
                                        if(rk_double( &model->rndstate) < (alpha *model->beta[j]) / (alpha * model->beta[j] + (float)c -1.0)){
                                                M[i][j] = M[i][j] + 1; // number of times state i generated color j...
                                                total_M = total_M + 1;
                                        }
                                        //M(j,k) = M(j,k) + (rand() < (ialpha0 * ibeta(k)) / (ialpha0 * ibeta(k) + l - 1));
                                }
                                
                        }
                        //fprintf(stdout," %d(%d)", M[i][j],transition_counts[i][j]);
                        sum_M[j] += M[i][j];
                        sum_N[i] += transition_counts[i][j];
                }
                //fprintf(stdout,"\n");
        }
        //fprintf(stdout,"\n");
        sum = 0.0;
        model->beta[0] = 0;
        for(i = 1; i < last_state;i++){
                model->beta[i] = rk_gamma(&model->rndstate, sum_M[i], 1.0);
                sum += model->beta[i];
        }
	
        model->beta[last_state] =  rk_gamma(&model->rndstate, gamma, 1.0);
        sum += model->beta[last_state] ;
        for(i = 0; i <= last_state;i++){
                model->beta[i] /= sum;
        }
        /* ok done with beta now gamma and alpha  */
        w = supp[2];
        p = supp[3];
        s = supp[4];
        
        for(j = 0;j < iterations;j++){
                sum_s = 0.0f;
                sum_w = 0.0f;
                for(i = 0; i < last_state;i++){
                        w[i] = rk_beta(&model->rndstate, alpha+1.0, sum_N[i]);
                        p[i] = sum_N[i] / alpha;
                        p[i] = p[i] / (p[i]+1.0);
                        s[i] = rk_binomial(&model->rndstate, 1, p[i]);
                        sum_s += s[i];
                        sum_w += log(w[i]);
                }
                alpha = rk_gamma(&model->rndstate, model->alpha0_a + total_M - sum_s, 1.0 / (model->alpha0_b -sum_w));
        }
        model->alpha = alpha;

        /* Let's do gamma now...    */
        mu = 0.0f;
        pi_mu = 0.0f;
	
        
        for(j = 0;j < iterations;j++){
                mu =  rk_beta(&model->rndstate, gamma+1.0, total_M);
                pi_mu = 1.0 / (1.0 + (total_M * ( model->gamma_b - log(mu) )) / (model->gamma_a + last_state -1.0));
                if(rk_double( &model->rndstate) < pi_mu){
                        gamma = rk_gamma(&model->rndstate,model->gamma_a + last_state, 1.0 / (model->gamma_b - log(mu)));
                }else{
                        gamma = rk_gamma(&model->rndstate, model->gamma_a + last_state-1.0, 1.0 / (model->gamma_b - log(mu)));
                }
        }
        model->gamma = gamma;
        
        free_2d((void**) M);
        free_2d((void**) supp);
	
        return OK;
ERROR:
        free_2d((void**) M);
        free_2d((void**) supp);
        return FAIL;
}


struct ihmm_model* alloc_ihmm_model(int K, int L)
{
        struct ihmm_model* ihmm = NULL;
        int i;
        ASSERT(K>3, "No states requested");
        ASSERT(L>1, "No letters");

        MMALLOC(ihmm, sizeof(struct ihmm_model));
        
        
        ihmm->transition_counts = NULL;
        ihmm->emission_counts = NULL;
        ihmm->beta = NULL;
        ihmm->num_states = 0;
        ihmm->alloc_num_states = 16;
        ihmm->L = L;
        ihmm->alpha = IHMM_PARAM_PLACEHOLDER;
        ihmm->alpha0_a = IHMM_PARAM_PLACEHOLDER;
        ihmm->alpha0_b = IHMM_PARAM_PLACEHOLDER;
        ihmm->gamma = IHMM_PARAM_PLACEHOLDER;
        ihmm->gamma_a = IHMM_PARAM_PLACEHOLDER;
        ihmm->gamma_b = IHMM_PARAM_PLACEHOLDER;

        while(K > ihmm->alloc_num_states){
                ihmm->alloc_num_states = ihmm->alloc_num_states << 1;
        }

        RUNP(ihmm->transition_counts = malloc_2d_float(ihmm->transition_counts, ihmm->alloc_num_states, ihmm->alloc_num_states, 0.0f));
        RUNP(ihmm->emission_counts = malloc_2d_float(ihmm->emission_counts , ihmm->L, ihmm->alloc_num_states, 0.0f));

        MMALLOC(ihmm->beta,sizeof(float) * ihmm->alloc_num_states);
        for(i = 0; i < ihmm->alloc_num_states;i++){
                ihmm->beta[i] = 0.0f;
        }
        
        return ihmm;
ERROR:
        free_ihmm_model(ihmm);
        return NULL;
}


int resize_ihmm_model(struct ihmm_model* ihmm, int K)
{
        int old_size;
        int i;
        ASSERT(ihmm != NULL, "No model");
        ASSERT(K > 3,"No states requested");

        old_size = ihmm->alloc_num_states;
        if(K >ihmm->alloc_num_states ){
                while(K > ihmm->alloc_num_states){
                        ihmm->alloc_num_states = ihmm->alloc_num_states << 1;
                }
                LOG_MSG("Resizing model to %d states",ihmm->alloc_num_states);
                RUNP(ihmm->transition_counts = malloc_2d_float(ihmm->transition_counts, ihmm->alloc_num_states, ihmm->alloc_num_states, 0.0f));
                RUNP(ihmm->emission_counts = malloc_2d_float(ihmm->emission_counts , ihmm->L, ihmm->alloc_num_states, 0.0f));
           
                MREALLOC(ihmm->beta,sizeof(float) * ihmm->alloc_num_states);
                for(i = old_size;i < ihmm->alloc_num_states;i++){
                        ihmm->beta[i] = -1;
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
                free_2d((void**) ihmm->transition_counts);
                free_2d((void**) ihmm->emission_counts);

                MFREE(ihmm->beta);
                MFREE(ihmm);
        }
}



int print_counts(struct ihmm_model* ihmm)
{
        int i,j;
        fprintf(stdout,"Transition counts\n"); 
        for(i = 0; i < ihmm->num_states;i++){
                fprintf(stdout,"s%3d", i ); 
                for(j = 0; j < ihmm->num_states;j++){
                         fprintf(stdout," %0.0f", ihmm->transition_counts[i][j]); 
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        fprintf(stdout,"Emission counts\n");

        for(i = 0; i < ihmm->L;i++){
                fprintf(stdout,"s%3d", i ); 
                for(j = 0; j < ihmm->num_states;j++){
                         fprintf(stdout," %0.0f", ihmm->emission_counts[i][j]); 
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        
        return OK;
}

int print_model_parameters(struct ihmm_model* ihmm)
{
        int i;
        float sum = 0.0;
        fprintf(stdout,"NUMBER of STATES: %d\n", ihmm->num_states);
        fprintf(stdout,"%3.3f alpha\t%3.3f gamma",ihmm->alpha , ihmm->gamma);
        for(i = 0; i < ihmm->num_states;i++){
                fprintf(stdout," %3.3f",ihmm->beta[i]);
                sum += ihmm->beta[i];
        }

        fprintf(stdout,"\tsum: %3.3f\n",sum);
        return OK;
}

#ifdef ITESTMODEL
int main(const int argc,const char * argv[])
{
        struct ihmm_model* ihmm = NULL;
        struct seq_buffer* iseq = NULL;
        int i;
        char *tmp_seq[4] = {
                "ACGT",
                "ACGT",
                "ACGT",
                "ACGT"};

        
        
        //119l
        RUN(print_program_header((char * const*)argv,"GAGA"));
        
        
        RUNP(iseq = create_ihmm_sequences_mem(tmp_seq ,4));
        RUN(random_label_ihmm_sequences(iseq, 10));

        
        RUNP(ihmm = alloc_ihmm_model(20, 4+3));
        /* Need to set alpha_[a/b] and gamma[a/b] manually before calling
         * hyper */
        RUN(resize_ihmm_model(ihmm, 16+3));
        RUN(resize_ihmm_model(ihmm, 2+3));

        /* At this stage there are no counts AND the model parameters are not
         * initialized. */  
        RUN(fill_counts(ihmm,iseq));
        RUN(print_counts(ihmm));

       
        
        /* I am doing this as a pre-caution. I don't want the inital model
         * contain states that are not visited.. */
        RUN(remove_unused_states_labels(ihmm, iseq));
        RUN(fill_counts(ihmm,iseq));
        RUN(print_counts(ihmm));
        /* Now there are counts but no model parameters. */
        ihmm->alpha0_a = 4.0f;
        ihmm->alpha0_b = 2.0f;
        ihmm->gamma_a = 3.0f;
        ihmm->gamma_b = 6.0f;
        ihmm->alpha = IHMM_PARAM_PLACEHOLDER;
        ihmm->gamma = IHMM_PARAM_PLACEHOLDER;
        RUN(iHmmHyperSample(ihmm, 10));
        /* Now I should have everything ready to go.  */
        RUN(print_model_parameters(ihmm));
        /* Just to verify if everything works..  */
        RUN(remove_unused_states_labels(ihmm, iseq));
        RUN(fill_counts(ihmm,iseq));
        RUN(print_counts(ihmm));
        
        /* Verify (by eye!) if the estimation of the hyperparameters works  */
        for(i = 0; i < 10;i++){
                ihmm->alpha0_a = 4.0f;
                ihmm->alpha0_b = 2.0f;
                ihmm->gamma_a = 3.0f;
                ihmm->gamma_b = 6.0f;
                ihmm->alpha = IHMM_PARAM_PLACEHOLDER;
                ihmm->gamma = IHMM_PARAM_PLACEHOLDER;
                RUN(iHmmHyperSample(ihmm, 10));
                RUN(print_model_parameters(ihmm));
        }
        RUN(print_counts(ihmm));
        RUN(write_model(ihmm, "test_model_file.txt"));
        free_ihmm_model(ihmm);

        RUNP(ihmm = read_model("test_model_file.txt"));
        LOG_MSG("After reading:");
        RUN(print_model_parameters(ihmm));
        RUN(print_counts(ihmm));
        
        free_ihmm_model(ihmm);
        free_ihmm_sequences(iseq);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
        
}
#endif
