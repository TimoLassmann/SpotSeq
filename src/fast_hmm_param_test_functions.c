#include "fast_hmm_param_test_functions.h"



int print_fast_hmm_params(struct fast_hmm_param* ft)
{
        int i,j;
        double** m = NULL;

        double sum = 0.0;



        ASSERT(ft != NULL, "No param");

        //RUNP(m = malloc_2d_float(m, ft->last_state+1,  ft->last_state+1, 0.0));
        fprintf(stdout,"States: %d\n",ft->last_state+1 );
        fprintf(stdout,"Transitions:\n");

        m = ft->transition;
        for(i = 0; i< ft->last_state+1;i++){
                fprintf(stdout,"S%d",i);
                sum = 0.0;
                for(j = 0; j< ft->last_state+1;j++){
                        fprintf(stdout," %0.3f",m[i][j]);
                        sum+= m[i][j];
                }
                fprintf(stdout,"\ts:%f\n",sum);
        }
        fprintf(stdout,"\n");

        fprintf(stdout,"Emission:\n");

        for(i = 0; i < ft->L;i++){

                sum = 0.0;
                for(j = 0; j< ft->last_state+1;j++){
                        fprintf(stdout," %0.3f",ft->emission[i][j]);
                        sum += ft->emission[i][j];
                }
                fprintf(stdout,"\ts:%f\n",sum);
        }
        LOG_MSG("Print infinity transitions.");
        for(j = 0; j< ft->last_state+1;j++){
                fprintf(stdout,"%d->%d %f\n", ft->infinity[j]->from, ft->infinity[j]->to, ft->infinity[j]->t);
        }
        LOG_MSG("Done.");
        return OK;
ERROR:
        return FAIL;
}

int fill_with_random_transitions(struct fast_hmm_param* ft, int k)
{
        //struct fast_t_item** list = NULL;

        struct fast_t_item* tmp = NULL;
        int i,j;
        //int num;
        float sum = 0;
        float* tmp_probs = NULL;
        ASSERT(ft != NULL, "No ft.");
        MMALLOC(tmp_probs, sizeof(float) * k);

        RUN(expand_ft_if_necessary(ft, k));

        rk_state rndstate;


        rk_randomseed(&rndstate);
        //num = ft->num_items;
        //list = ft->list;

        for(i = 0;i < k;i++){
                sum = 0.0;
                for(j = 0;j < k;j++){

                        tmp_probs[j] = (float) rk_double(&rndstate);// random_float_zero_to_x(1.0);


                        sum+=tmp_probs[j];

                }
                for(j = 0;j < k;j++){
                        tmp_probs[j] /= sum;
                        tmp = NULL;
                        MMALLOC(tmp, sizeof(struct fast_t_item));
                        tmp->from = i;
                        tmp->to = j;
                        tmp->t = tmp_probs[j];

                        ft->root->tree_insert(ft->root,tmp);
                        ft->transition[i][j] = tmp_probs[j];
                }
        }
        ft->last_state = k-1;
        MFREE(tmp_probs);

        return OK;
ERROR:
        return FAIL;
}
