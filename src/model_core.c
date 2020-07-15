

#include "global.h"
#include "tldevel.h"
#include "distributions.h"

#include "finite_hmm.h"
#include <math.h>
#include <stdint.h>

#include "sequence_struct.h"
#include "model_alloc.h"
//#include "ihmm_seq.h"
#include "finite_hmm.h"


#include "null_model_emission.h"


#define MODEL_CORE_IMPORT
#include "model_core.h"


static int fill_counts_i(struct ihmm_model* ihmm, struct ihmm_sequence* s, int model_index );

//static int label_seq_based_on_random_fhmm(struct seq_buffer* sb, int k, double alpha);






int esl_stats_LogGamma(double x, double *ret_answer)
{
        int i;
        double xx, tx;
        double tmp, value;
        static double cof[11] = {
                4.694580336184385e+04,
                -1.560605207784446e+05,
                2.065049568014106e+05,
                -1.388934775095388e+05,
                5.031796415085709e+04,
                -9.601592329182778e+03,
                8.785855930895250e+02,
                -3.155153906098611e+01,
                2.908143421162229e-01,
                -2.319827630494973e-04,
                1.251639670050933e-10
        };

        /* Protect against invalid x<=0  */
        //if (x <= 0.0)  ESL_EXCEPTION(eslERANGE, "invalid x <= 0 in esl_stats_LogGamma()");

        xx       = x - 1.0;
        tx = tmp = xx + 11.0;
        value    = 1.0;
        for (i = 10; i >= 0; i--)	/* sum least significant terms first */
        {
                value += cof[i] / tmp;
                tmp   -= 1.0;
        }
        value  = log(value);
        tx    += 0.5;
        value += 0.918938533 + (xx+0.5)*log(tx) - tx;
        *ret_answer = value;
        return OK;
}

int log_likelihood_model(struct ihmm_model* model, struct seq_buffer* sb)
{

        int i,j;
        double l;
        double* tmp = NULL;
        double sum = 0;
        double r = 0;
        ASSERT(model != NULL, "No model.");
        ASSERT(sb != NULL, "No sequences");


        for(i = 0; i < sb->num_seq;i++){
                RUN(clear_counts(model));
                RUN(fill_counts_i(model, sb->sequences[i],0));

        }
        RUN(fill_counts(model,sb,0));

        MMALLOC(tmp, sizeof(double) * model->num_states);
        l = 0;
        for(i = 0; i < model->num_states;i++){
                sum = 0.0;

                for(j = 0; j < model->num_states;j++){
                        tmp[j] = model->transition_counts[i][j] + model->beta[j] * model->alpha;
                        sum += tmp[j];



                }
                // + gammaln(alpha0) ...
                RUN(esl_stats_LogGamma(model->alpha,&r  ));
                l += r;

                // - gammaln(sum([N(k,:) 0]) + alpha0) ...

                RUN(esl_stats_LogGamma(sum+model->alpha,&r  ));

                //+ sum(gammaln( R(nzind)  )) ...
                l -= r;
                sum = 0.0;
                for(j = 0; j < model->num_states;j++){

                        if(tmp[j]){
                                RUN(esl_stats_LogGamma(tmp[j],&r  ));
                                sum += r;
                        }

                }
                l += r;

                //- sum(gammaln( ab(nzind) ));
                sum = 0.0;
                for(j = 0; j < model->num_states;j++){
                        if(model->beta[j] && model->alpha){
                                RUN(esl_stats_LogGamma(model->beta[j] * model->alpha,&r  ));
                                sum += r;
                        }
                }
                l -= r;

//- sum(gammaln( ab(nzind) ));

//this is simpkly looking at observation  (P(i ,j ,c)) based on counts + hyper d

                //model->alpha
                //model->alpha;





        }
        // for(i = 0; )
        // //  tmp_prob[j] = rk_gamma(&model->rndstate, model->transition_counts[i][j] + model->beta[j] * model->alpha,1.0);
        //

        return OK;
ERROR:
        return FAIL;
}

int remove_unused_states_labels(struct ihmm_model* ihmm, struct seq_buffer* sb, int model_index)
{
        int i,j;
        double sum;
        int len;
        int* relabel = NULL;
        int* used = NULL;
        uint16_t* lab = NULL;

        ASSERT(ihmm != NULL, "no model");
        ASSERT(sb != NULL, "no seq struct");
        MMALLOC(relabel, sizeof(int) * ihmm->num_states);
        MMALLOC(used, sizeof(int) * ihmm->num_states);
        for(i = 0; i < ihmm->num_states;i++){
                used[i] = 0;
                relabel[i] = -1;
        }
        used[START_STATE] = 100;
        used[END_STATE] = 100;
        int max = -1;
        for(i = 0; i < sb->num_seq;i++){
                lab = sb->sequences[i]->label_arr[model_index];
                len = sb->sequences[i]->seq_len;
                for(j = 0; j < len;j++){
                        //LOG_MSG("model:%d seq:%d pos:%d label:%d",model_index,i,j,lab[j]);
                        used[lab[j]]++;
                        if(lab[j] > max){
                                max= lab[j];
                        }
                }
        }



        /*fprintf(stdout,"ORG beta \n");
        for(i = 0; i < ihmm->num_states;i++){
                fprintf(stdout,"%3.3f ", ihmm->beta[i]);
        }
        fprintf(stdout,"\n");*/

        j = 0;
        sum = 0.0;
        for(i = 0; i < ihmm->num_states;i++){
                if(used[i] > sb->num_seq /2){  //} != 0){
                        ihmm->beta[j] = ihmm->beta[i];
                        relabel[i] = j;
                        j++;
                }else{
                        relabel[i] = j;
                        sum += ihmm->beta[i];
                }
        }
        /*fprintf(stdout,"USED::::\n");
        for(i = 0; i < ihmm->num_states;i++){
                fprintf(stdout,"%3d ",used[i]);
        }
        fprintf(stdout,"\n");
        for(i = 0; i < ihmm->num_states;i++){
                fprintf(stdout,"%3d ",relabel[i]);
        }
        fprintf(stdout,"\n");
        */
        ihmm->beta[j] = sum;
        //for(i = j+1; i < ihmm->num_states;i++){
        //       ihmm->beta[i] = 0;
        //}
        ihmm->num_states = j+1; /* need to add one for the infinite stuff */

        RUN(resize_ihmm_model(ihmm, j+1));
        /*fprintf(stdout,"CUR beta \n");
        sum = 0.0f;
        for(i = 0; i < ihmm->num_states;i++){
                        fprintf(stdout,"%3.3f ", ihmm->beta[i]);
                sum+= ihmm->beta[i];
        }

        fprintf(stdout,"\tsum: %f\n",sum);*/


        for(i = 0; i < sb->num_seq;i++){
                //fprintf(stdout,"%3d",i);
                lab = sb->sequences[i]->label_arr[model_index];
                len = sb->sequences[i]->seq_len;
                for(j= 0; j <  len;j++){
                        lab[j] = relabel[lab[j]];
                        //      fprintf(stdout," %d",lab[j]);
                }
                //fprintf(stdout,"\n");
        }
        MFREE(used);
        MFREE(relabel);
        return OK;
ERROR:
        MFREE(used);
        MFREE(relabel);
        return FAIL;
}

int fill_counts(struct ihmm_model* ihmm, struct seq_buffer* sb, int model_index)
{
        int i,j;
        uint16_t* label = NULL;
        int max_state_ID;
        int len;
        ASSERT(ihmm != NULL,"No model.");
        ASSERT(sb != NULL,"No iseq struct");

        /* First I need to check what the largest state ID is and see if we have sufficient space allocated in the model.  */
        max_state_ID = -1;
        //LOG_MSG("SEQ: %d",sb->)
        //fprintf(stdout,"%d numseq\n",sb->num_seq );
        for(i = 0; i < sb->num_seq;i++){
                label = sb->sequences[i]->label_arr[model_index];
                len = sb->sequences[i]->seq_len;
                for(j = 0; j < len;j++){
                        //fprintf(stdout,"%d ",label[j]);
                        if(label[j] > max_state_ID){
                                max_state_ID = label[j];
                        }
                }
                //fprintf(stdout,"\n");
        }
        max_state_ID += 1; // for the infinity possibility; not observed in the current labeling
        max_state_ID += 1; // so I can use the < syntax rather than <=
        ASSERT(max_state_ID > 2, "Not enough states found");

        //
        //LOG_MSG("MAX_STATE_ID:  %d", max_state_ID);
        //exit(0);
        //RUN(resize_ihmm_model(ihmm, max_state_ID));
        ihmm->num_states = max_state_ID;

        /* clear transition counts */
        /* clear emission counts */
        RUN(clear_counts(ihmm));
        //print_counts(ihmm);

        for(i = 0; i < sb->num_seq;i++){
                RUN(fill_counts_i(ihmm, sb->sequences[i],model_index));
        }
        return OK;
ERROR:
        return FAIL;
}



int add_pseudocounts_emission(struct ihmm_model* model, double alpha)
{
        double* background = NULL;
        int i,j;
        double sum;
        ASSERT(model != NULL, "No model.");


        background = model->background;

        for(i = 2; i < model->num_states;i++){
                sum = 0.0;
                for(j = 0; j < model->L;j++){
                        sum += model->emission_counts[j][i];
                }
                if(sum){
                        for(j = 0; j < model->L;j++){
                                model->emission_counts[j][i] += alpha * background[j];
                        }
                }

        }
        return OK;
ERROR:
        return FAIL;

}

int fill_counts_i(struct ihmm_model* ihmm, struct ihmm_sequence* s, int model_index)
{
        uint16_t* label = NULL;
        uint8_t* seq = NULL;
        double** e = NULL;
        double** m = NULL;
        double score = 0.0;
        //float* u = NULL;
        //float r;
        int len;
        int i;


        ASSERT(ihmm != NULL,"no model");

        label = s->label_arr[model_index];
        seq = s->seq;
        len = s->seq_len;
        score = s->score_arr[model_index];


        // 1.0;// - scaledprob2prob(s->score);
        //score =1.0 - scaledprob2prob ( s->score - logsum(s->score,  s->r_score));
        //u = sb->sequences[seq_ID]->u;
        //fprintf(stdout,"%f %f %f\n", s->score,s->r_score, scaledprob2prob ( s->score - logsum(s->score,  s->r_score)));

        e = ihmm->emission_counts;
        m = ihmm->transition_counts;

        m[START_STATE][label[0]]  += score;
        e[(int)seq[0]][label[0]]+= score;

        for(i = 1; i < len;i++){
                m[label[i-1]][label[i]] += score;
                e[(int)seq[i]][label[i]] += score;
                //e[(int)seq[i]][label[i]] += scaledprob2prob(u[i]);
                //r = rk_double(&ihmm->rndstate);
                //m[label[i-1]][label[i]] += (r - 0.5f);
                //r = rk_double(&ihmm->rndstate);
                //e[(int)seq[i]][label[i]]+=  (r - 0.5f);//r / 1.0;

        }
        m[label[len-1]][END_STATE] += score;
        // r = rk_double(&ihmm->rndstate);
        // m[label[len-1]][IHMM_END_STATE] += (r*2 - 1.0f);// r / 1.0;
        //LOG_MSG("SCORE: %f", score);
        return OK;
ERROR:
        return FAIL;
}

int iHmmHyperSample(struct ihmm_model* model, int iterations)
{
        int i,j,c;
        int last_state;
        double** transition_counts = NULL;
        double** M = NULL;
        double** supp = NULL;
        double* sum_M = NULL;
        double* sum_N = NULL;
        double* w = NULL;
        double* p = NULL;
        double* s = NULL;
        double total_M = 0.0;
        double alpha;
        double gamma;
        double sum, sum_s, sum_w, mu, pi_mu;

        ASSERT(model != NULL, "No model");
        ASSERT(iterations > 0, "No iterations");

        last_state = model->num_states-1;

        /* alloc auxillary data structures  */
        //RUNP(M = galloc(M, model->num_states,model->num_states, 0.0));
        RUN(galloc(&M, model->num_states,model->num_states));
        for(i = 0; i < model->num_states;i++){
                for(j = 0; j < model->num_states;j++){
                        M[i][j] = 0.0;
                }
        }
        //RUNP(supp = galloc(supp, 5,model->num_states, 0.0));
        RUN(galloc(&supp, 5,model->num_states));
        for(i = 0; i < 5;i++){
                for(j = 0; j < model->num_states;j++){
                        supp[i][j] = 0.0;
                }
        }
        //fprintf(stdout,"%f %f\n", model->alpha ,model->gamma );
        alpha = model->alpha;
        gamma = model->gamma;
        //gamma = 10;
        //WAS here - need to alloc auxillary arrays for calculations below!
        transition_counts = model->transition_counts;
        sum_M = supp[0];
        sum_N = supp[1];

        //LOG_MSG("%d %d  LAST:%d",  model->rndstate.pos, model->rndstate.key[model->rndstate.pos],last_state);

        total_M = 0.0f;
        for(i = 0; i < last_state;i++){
                for(j = 0; j < last_state;j++){
                        if(transition_counts[i][j] == 0){
                                M[i][j] = 0;
                        }else{
                                for(c = 1;c <= transition_counts[i][j];c++){
                                        if(rk_double( &model->rndstate) < (alpha * model->beta[j]) / (alpha * model->beta[j] + (double)c -1.0)){
                                                M[i][j] = M[i][j] + 1; // number of times state i generated color j...
                                                total_M = total_M + 1.0;
                                        }
                                }

                        }
                        //fprintf(stdout," %d(%d)", M[i][j],transition_counts[i][j]);
                        sum_M[j] += M[i][j];
                        sum_N[i] += transition_counts[i][j];
                }
                //fprintf(stdout,"\n");
        }
        //fprintf(stdout,"\n");
        //LOG_MSG("%d %d %f %d",  model->rndstate.pos, model->rndstate.key[model->rndstate.pos] , model->rndstate.gauss, model->rndstate.has_gauss);
        sum = 0.0;
        model->beta[0] = 0;
        for(i = 1; i < last_state;i++){
                //LOG_MSG("sumM%d %f ",i,sum_M[i]);
                model->beta[i] = rk_gamma(&model->rndstate, sum_M[i], 1.0);
                //LOG_MSG("%d %d %f %d",  model->rndstate.pos, model->rndstate.key[model->rndstate.pos] , model->rndstate.gauss, model->rndstate.has_gauss);
                sum += model->beta[i];
        }
        model->beta[last_state] =  rk_gamma(&model->rndstate, gamma, 1.0);
        //LOG_MSG("%d %d",  model->rndstate.pos, model->rndstate.key[model->rndstate.pos]);
        sum += model->beta[last_state] ;
        for(i = 0; i <= last_state;i++){
                model->beta[i] /= sum;
        }
        //LOG_MSG("%d %d",  model->rndstate.pos, model->rndstate.key[model->rndstate.pos]);
        //LOG_MSG("%f %f %f %f %f %f %f", total_M, gamma,model->gamma_a,model->gamma_b,alpha,model->alpha_a,model->alpha_b);
        /* Only re-estimate alpha and gamma if vague priors are set...  */
        if(model->alpha_a != IHMM_PARAM_PLACEHOLDER){
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
                        alpha = rk_gamma(&model->rndstate, model->alpha_a + total_M - sum_s, 1.0 / (model->alpha_b -sum_w));
                }
                if(alpha >= model->alpha_limit){
                        model->alpha = model->alpha_limit;
                }else{
                        model->alpha = alpha;
                }
        }
        if(model->gamma_a != IHMM_PARAM_PLACEHOLDER){
                /* Let's do gamma now...    */
                mu = 0.0f;
                pi_mu = 0.0f;
                //LOG_MSG("%d %f %f %f", total_M, gamma,model->gamma_a,model->gamma_b );
                for(j = 0;j < iterations;j++){
                        //LOG_MSG("%d %d",  model->rndstate.pos, model->rndstate.key[model->rndstate.pos]);
                        mu =  rk_beta(&model->rndstate, gamma+1.0, total_M);
                        //LOG_MSG("%f mu",mu);
                        pi_mu = 1.0 / (1.0 + (total_M * ( model->gamma_b - log(mu) )) / (model->gamma_a + (double) last_state -1.0));
                        if(rk_double( &model->rndstate) < pi_mu){
                                gamma = rk_gamma(&model->rndstate, model->gamma_a + (double) last_state, 1.0 / (model->gamma_b - log(mu)));
                        }else{
                                gamma = rk_gamma(&model->rndstate, model->gamma_a + (double) last_state-1.0, 1.0 / (model->gamma_b - log(mu)));
                        }
                }
                if(gamma >= model->gamma_limit){
                        model->gamma = model->gamma_limit;
                }else{
                        model->gamma = gamma;
                }

        }
        gfree(M);
        gfree(supp);

        return OK;
ERROR:
        gfree(M);
        gfree(supp);
        return FAIL;
}


int set_model_hyper_parameters(struct model_bag* b, double alpha, double gamma)
{
        struct ihmm_model* model = NULL;
        int i,j;
        ASSERT(b != NULL, "No model bag");

        for(i = 0; i < b->num_models;i++){
                model = b->models[i];
                //if(alpha == IHMM_PARAM_PLACEHOLDER){
                        model->alpha_a = 6.0;
                        model->alpha_b = 15.0;
                        model->alpha = rk_gamma(&model->rndstate, model->alpha_a,1.0 / model->alpha_b);
                        //}else{
                        //        model->alpha = alpha;
                        //}

                //if(gamma == IHMM_PARAM_PLACEHOLDER){
                        model->gamma_a = 16.0;
                        model->gamma_b = 4.0;
                        model->gamma = rk_gamma(&model->rndstate, model->gamma_a,1.0 / model->gamma_b);
                        //}else{
                        //        model->gamma = gamma;
                        //}
                        //model->gamma_limit = 50.0;
                        //model->alpha_limit = 2.0;
                if(alpha != IHMM_PARAM_PLACEHOLDER){
                        model->alpha_limit  = alpha;
                }

                if(gamma != IHMM_PARAM_PLACEHOLDER){
                        model->gamma_limit = gamma;
                }

                for(j = 0; j < model->num_states;j++){
                        model->beta[j] = (double)(model->num_states);
                }
        }
        return OK;
ERROR:
        return FAIL;
}





int clear_counts(struct ihmm_model* ihmm)
{
        int i,j;

        ASSERT(ihmm != NULL, "No model");
        for(i = 0;i < ihmm->num_states;i++){
                for(j = 0; j < ihmm->num_states;j++){
                        ihmm->transition_counts[i][j] = 0;
                }
        }
        for(i = 0;i < ihmm->L ;i++){
                for(j = 0; j < ihmm->num_states;j++){
                        ihmm->emission_counts[i][j] = 0;
                }
        }


        return OK;
ERROR:
        return FAIL;
}




