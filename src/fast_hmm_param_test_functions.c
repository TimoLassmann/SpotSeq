#include "fast_hmm_param_test_functions.h"



int print_fast_hmm_params(struct fast_hmm_param* ft)
{
        int i,j;
        float** m = NULL;

        float sum = 0.0;
        RUNP(m = malloc_2d_float(m, ft->last_state+1,  ft->last_state+1, 0.0));
        fprintf(stdout,"Transitions:\n");
        for(i = 0; i < ft->num_items;i++){
                //fprintf(stdout, "%d) %d -> %d = %f\n", i,ft->list[i]->from,ft->list[i]->to,ft->list[i]->t);
                if(ft->list[i]->t != -1){
                        m[ft->list[i]->from][ft->list[i]->to] = ft->list[i]->t;
                }
        }

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

      
        for(j = 0; j< ft->last_state+1;j++){
                sum = 0.0;
                for(i = 0; i < 4;i++){
                        fprintf(stdout," %0.3f",ft->emission[i][j]);
                        sum += ft->emission[i][j];
                }
                fprintf(stdout,"\ts:%f\n",sum);
        }
    

        
        
        free_2d((void**)m);
        return OK;
ERROR:
        return FAIL;
}

int fill_with_random_transitions(struct fast_hmm_param* ft, int k)
{
        struct fast_t_item** list = NULL;
        int i,j;
        int num;
        float sum = 0;
        ASSERT(ft != NULL, "No ft.");

        RUN(expand_fast_hmm_param_if_necessary(ft, k));
        
        num = ft->num_items;
        list = ft->list;
 
        for(i = 0;i < k;i++){
                sum = 0.0;
                for(j = 0;j < k;j++){
                        list[num]->from = i;
                        list[num]->to = j;
                        list[num]->t = random_float_zero_to_x(1.0);
                        sum+=list[num]->t;
                        num++;
                }
                for(j = 0;j < k;j++){
                        list[num-k+j]->t /= sum;
                }
        }
        ft->num_items = num;
        ft->last_state = k-1;
        return OK;
ERROR:
        return FAIL;
}
