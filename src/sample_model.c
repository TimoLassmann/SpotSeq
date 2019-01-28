#include "model.h"



double** sample_transistion(double**t,struct ihmm_model* model)
{
        int last_state;

        double* tmp_prob = NULL;
        double sum;
        int i,j;
        ASSERT(model != NULL, "No model");


        MMALLOC(tmp_prob, sizeof(double) *(model->num_states));

        last_state = model->num_states -1;

        t = galloc(t,model->num_states,model->num_states,0.0);

        /* Let's start with the transitions from the start state!  */
        //last_index = list_index;
        sum = 0.0;

        /* Disallow Start to start transitions */
        // insert into transition matrix.
        t[IHMM_START_STATE][IHMM_START_STATE] = 0.0;
        /* Disallow Start to end transitions i.e. zero length sequences are not allowed*/
        t[IHMM_START_STATE][IHMM_END_STATE] = 0.0;
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
        for(i = 2; i <= last_state;i++){
                t[IHMM_START_STATE][i] = tmp_prob[i] / sum;
        }
        /* And now the stop state...  */
        /* There is no possibility escape the end state - all transitions from
         * end are zero. I am not sure if this matters in my dyn prig. code but
         * why not! */
        for(i = 0; i < last_state;i++){
                t[IHMM_END_STATE][i] = 0.0;
        }
        t[IHMM_END_STATE][last_state] = 0.0;


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

                for(j = 1; j <= last_state;j++){
                        t[i][j] = tmp_prob[j] / sum;
                }
        }
        MFREE(tmp_prob);
        return t;
ERROR:
        gfree(t);
        return NULL;
}

double** sample_emission(double**e,struct ihmm_model* model)
{
        int i,j;
        double sum;
        int last_state;
        ASSERT(model!= NULL, "No model");
        last_state = model->num_states -1;


        RUNP(e = galloc(e,  model->L, model->num_states, 0.0));

        for(i = 0; i < model->L;i++){
                for(j = 0; j < last_state;j++){
                        e[i][j] = rk_gamma(&model->rndstate, model->emission_counts[i][j] + EMISSION_H, 1.0);
                }
        }

        for(j = 0; j < last_state;j++){
                sum = 0.0;
                for(i = 0; i < model->L;i++){
                        sum += e[i][j];
                }
                if(sum == 0.0){
                        sum = 1.0;
                }
                for(i = 0; i < model->L;i++){
                        e[i][j] /= sum;
                }
        }
        return e;
ERROR:
        gfree(e);
        return NULL;
}
