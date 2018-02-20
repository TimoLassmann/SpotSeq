
#include "fast_hmm_param.h"

/* The Idea here is to replace the two dimensional datastructure used for
 * transitions with a linear entity so that we can sort / add / remove states
 * and transitions easier. */

/* I wan function to do: */
/*         - stick breaking operation */
/*         - sort based on probability */
/*         - remove states  */



/* Housekeeping function  */
/* - alloc, relic , free*/
/* - sort */
/* - bin_search upper lower */





struct fast_hmm_param* alloc_fast_hmm_param(int k, int L)
{
        struct fast_hmm_param* ft = NULL;
        int i;
        
        ASSERT(L > 1, "Need more than one letter");
        MMALLOC(ft, sizeof(struct fast_hmm_param));
        ft->alloc_num_states = k;
        ft->alloc_items = ft->alloc_num_states* ft->alloc_num_states ;
        ft->last_state = 0;
        ft->active_states = NULL;
        ft->list = NULL;
        ft->num_items = 0;
        ft->emission = NULL;    /* This will be indexed by letter i.e. e['A']['numstate'] */
        ft->L = L;

        RUNP(ft->emission = malloc_2d_float(ft->emission, ft->L, ft->alloc_num_states, 0.0f));

        MMALLOC(ft->active_states, sizeof(int8_t)* ft->alloc_num_states);

        for(i = 0; i < ft->alloc_num_states;i++){
                ft->active_states[i] =0;
        }


        
        MMALLOC(ft->list, sizeof(struct fast_t_item*) * ft->alloc_items);

        for(i = 0; i < ft->alloc_items;i++){
                ft->list[i] = NULL;
                MMALLOC(ft->list[i], sizeof(struct fast_t_item));
                ft->list[i]->from = -1;
                ft->list[i]->to = -1;
                ft->list[i]->t = 0.0f;
        }        
        return ft;
ERROR:
        free_fast_hmm_param(ft);
        return NULL;
}

int expand_fast_hmm_param_if_necessary(struct fast_hmm_param* ft, int k)
{
        int i, cur_k;
        ASSERT(ft != NULL, "No ft struct!");
        ASSERT(k >2,"No states requested");
        cur_k = ft->alloc_num_states;
        
        if(k > ft->alloc_num_states){
                while(k > ft->alloc_num_states){
                        ft->alloc_num_states = ft->alloc_num_states + 64;
                }

                RUNP(ft->emission = malloc_2d_float(ft->emission, ft->L, ft->alloc_num_states, 0.0f));
                
                MREALLOC(ft->active_states, sizeof(int8_t)* ft->alloc_num_states);
                for(i = cur_k; i < ft->alloc_num_states;i++){
                        ft->active_states[i] =0;
                }
                
                ft->alloc_items = ft->alloc_num_states* ft->alloc_num_states ;
                
                MREALLOC(ft->list, sizeof(struct fast_t_item*) * ft->alloc_items);
                cur_k = cur_k * cur_k;
                for(i = cur_k; i < ft->alloc_items;i++){
                        ft->list[i] = NULL;
                        MMALLOC(ft->list[i], sizeof(struct fast_t_item));
                        ft->list[i]->from = -1;
                        ft->list[i]->to = -1;
                        ft->list[i]->t = 0.0f;
                }        
                
        }
        return OK;        
ERROR:
        free_fast_hmm_param(ft);
        return FAIL;
}

void free_fast_hmm_param(struct fast_hmm_param* ft)
{
        int i;
        if(ft){
                if(ft->list){
                        for(i = 0; i < ft->alloc_items;i++){
                                MFREE(ft->list[i]);
                        }
                        MFREE(ft->list);
                }
                if(ft->emission){
                        free_2d((void**) ft->emission);
                }
                if(ft->active_states){
                        MFREE(ft->active_states);
                }
                MFREE(ft);
        }
}






int fast_hmm_param_cmp_by_t_desc(const void *a, const void *b) 
{ 
    struct fast_t_item* const *one = a;
    struct fast_t_item* const *two = b;

    if((*one)->t > (*two)->t){
            return -1;            
    }else if((*one)->t == (*two)->t){
            return 0;
    }else{
            return 1;
    }
}


int fast_hmm_param_cmp_by_to_from_asc(const void *a, const void *b) 
{ 
    struct fast_t_item* const *one = a;
    struct fast_t_item* const *two = b;

    if((*one)->from > (*two)->from){
            return 1;            
    }else if((*one)->from == (*two)->from){
            if((*one)->to > (*two)->to){
                    return 1;            
            }else if((*one)->to == (*two)->to){
                    return 0;
            }else{
                    return -1;
            }
    }else{
            return -1;
    }
}


int fast_hmm_param_cmp_by_from_asc(const void *a, const void *b) 
{ 
    struct fast_t_item* const *one = a;
    struct fast_t_item* const *two = b;

    if((*one)->from > (*two)->from){
            return 1;            
    }else if((*one)->from == (*two)->from){
            
            return 0;
           
    }else{
            return -1;
    }
}




int fast_hmm_param_cmp_by_to_asc(const void *a, const void *b) 
{ 
    struct fast_t_item* const *one = a;
    struct fast_t_item* const *two = b;


    if((*one)->to > (*two)->to){
            return 1;            
    }else if((*one)->to == (*two)->to){
            return 0;
    }else{
            return -1;
    }

}

/* Selects item so that 0 .. return value is greater than x */
int fast_hmm_param_binarySearch_t(struct fast_hmm_param* ft, float x)
{
        struct fast_t_item** list = NULL;
        int l,r;

        l = 0;
        r = ft->num_items -1;
        list = ft->list;
        while (l <= r)
        {
                int m = l + (r-l)/2;
 
                // Check if x is present at mid
                if (list[m]->t == x){
                        return m;
                }
 
                // If x greater, ignore left half
                if (list[m]->t > x){
                        l = m + 1;
                }else{
                        r = m - 1;
                }
        }        
        return  MACRO_MAX(l,r);
}



int fast_hmm_param_binarySearch_to_lower_bound(struct fast_hmm_param* ft, int x)
{
        struct fast_t_item** list = NULL;
        int l,r;
        l = 0;
        r = ft->num_items -1;
        list = ft->list;
        while (l <= r)
        {
                int m = l + (r-l)/2;
                if (list[m]->to < x){
                        l = m + 1;
                }else{
                        r = m -1;
                }
        }        
        return  l;
}

int fast_hmm_param_binarySearch_to_upper_bound(struct fast_hmm_param* ft, int x)
{
        struct fast_t_item** list = NULL;
        int l,r;

        l = 0;
        r = ft->num_items -1;
        list = ft->list;
        while (l <= r)
        {
                int m = l + (r-l)/2;
                if (x < list[m]->to){
                        r = m -1;
                }else{
                        l = m + 1;
                }
        }        
        return  l;
}

int fast_hmm_param_binarySearch_from_lower_bound(struct fast_hmm_param* ft, int x)
{
        struct fast_t_item** list = NULL;
        int l,r;
        l = 0;
        r = ft->num_items -1;
        list = ft->list;
        while (l <= r)
        {
                int m = l + (r-l)/2;
                if (list[m]->from < x){
                        l = m + 1;
                }else{
                        r = m -1;
                }
        }        
        return  l;
}

int fast_hmm_param_binarySearch_from_upper_bound(struct fast_hmm_param* ft, int x)
{
        struct fast_t_item** list = NULL;
        int l,r;

        l = 0;
        r = ft->num_items -1;
        list = ft->list;
        while (l <= r)
        {
                int m = l + (r-l)/2;
                if (x < list[m]->from){
                        r = m -1;
                }else{
                        l = m + 1;
                }
        }        
        return  l;
}






#ifdef ITESTFASTPARAM


#include "fast_hmm_param_test_functions.h"



int main(const int argc,const char * argv[])
{
        fprintf(stdout,"Hello world\n");
        struct fast_hmm_param* ft = NULL;
        int i;
        int res = 0;
        float x;

        RUN(print_program_header((char * const*)argv,"Fast HMM param test"));

        RUNP(ft = alloc_fast_hmm_param(4,4));

        RUN(fill_with_random_transitions(ft, 8));
        RUN(print_fast_hmm_params(ft));
        fprintf(stdout,"%d items\n",ft->num_items);
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*), fast_hmm_param_cmp_by_t_desc);
        
        RUN(print_fast_hmm_params(ft));
        for(i =0; i < 10;i++){
                x = random_float_zero_to_x(1.0);
                
                res = fast_hmm_param_binarySearch_t(ft, x);
                fprintf(stdout,"search for %f: %d  \n",x, res);
        }
        x = ft->list[3]->t;

        res = fast_hmm_param_binarySearch_t(ft,x);
        fprintf(stdout,"search for %f: %d   \n",x, res);

        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_from_asc);
        RUN(print_fast_hmm_params(ft));


        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_asc);
        RUN(print_fast_hmm_params(ft));
        for(i =0; i < 4;i++){
                fprintf(stdout,"%d: %d -> %d\n",i, fast_hmm_param_binarySearch_to_lower_bound(ft,i), fast_hmm_param_binarySearch_to_upper_bound(ft,i));
        }
        
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_from_asc);
        RUN(print_fast_hmm_params(ft));
        for(i =0; i < 4;i++){
                fprintf(stdout,"%d: %d -> %d\n",i, fast_hmm_param_binarySearch_from_lower_bound(ft,i), fast_hmm_param_binarySearch_from_upper_bound(ft,i));
        }

        
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*), fast_hmm_param_cmp_by_to_from_asc);
        RUN(print_fast_hmm_params(ft));
        

         
        RUN(expand_fast_hmm_param_if_necessary(ft, 4));

        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*), fast_hmm_param_cmp_by_to_from_asc);
        
       

        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*), fast_hmm_param_cmp_by_to_from_asc);
        RUN(print_fast_hmm_params(ft));
        
        float* beta = NULL;
        float gamma = 6;
        //float alpha = 1000.2;
        float sum;
        rk_state rndstate;
                
        rk_randomseed(&rndstate);
        MMALLOC(beta, sizeof(float) * 64);
        
        sum = 0.0;
        beta[0] = 0;
        for(i = 1; i <  ft->last_state;i++){
                beta[i] = rk_gamma(&rndstate, (float)i*10, 1.0);
                sum += beta[i];
        }
	
        beta[ft->last_state] =  rk_gamma(&rndstate, gamma, 1.0);
        //fprintf(stdout,"BETA inf:%f\n",model->beta[model->infinityghost]  );
        sum += beta[ft->last_state] ;
        for(i = 0; i <= ft->last_state;i++){
		
                beta[i] /= sum;
                //fprintf(stdout,"BETA: %d %f\n",i,beta[i] );
        }
        for(i = 0;i < 4;i++){
                //RUN(add_state_from_fast_hmm_param(rndstate,ft,  beta, alpha, gamma));
        }
        for(i = 0; i <= ft->last_state;i++){
		
              
                fprintf(stdout,"BETA: %d %f  AFTER\n",i,beta[i] );
                }
        
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*), fast_hmm_param_cmp_by_to_from_asc);
        RUN(print_fast_hmm_params(ft));

        
        MFREE(beta);
               
        
        free_fast_hmm_param(ft);
        return EXIT_SUCCESS;
ERROR:
        free_fast_hmm_param(ft);
        return EXIT_FAILURE;
}


#endif 
