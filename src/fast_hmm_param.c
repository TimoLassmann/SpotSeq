
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
#define iHMM_START_STATE 0
static int fast_hmm_param_cmp_by_t_desc(const void *a, const void *b);
static int fast_hmm_param_cmp_by_to_from_asc(const void *a, const void *b);
static int fast_hmm_param_cmp_by_from_asc(const void *a, const void *b);
static int fast_hmm_param_cmp_by_to_asc(const void *a, const void *b);

/* return index of first element < x i.e. we can then do for(i =0; i < return;i++) */

static int binarySearch_t(struct fast_hmm_param* ft, float x);

/* These functions return the first and last+1 entry in list that has value of x */
static int binarySearch_to_lower_bound(struct fast_hmm_param* ft, int x);
static int binarySearch_to_upper_bound(struct fast_hmm_param* ft, int x);      
static int binarySearch_from_lower_bound(struct fast_hmm_param* ft, int x);
static int binarySearch_from_upper_bound(struct fast_hmm_param* ft, int x);



/* This function assumes (oh no!) that beta has space for an additional
 * element */
int add_state_from_fast_hmm_param(rk_state rndstate,struct fast_hmm_param* ft, float* beta, float alpha,float gamma)
{
        struct fast_t_item** list = NULL;
        float sum,be,bg,pe,pg, a,b;
        int i,new_k,list_index;
        int l,r;



        ASSERT(ft != NULL, "No ft.");
        list_index = ft->num_items;
        /* First add empty space to host the newstate -> old state transitions. */
        if(list_index + ft->last_state + ft->last_state + 1 >= ft->alloc_num_states){
                RUN(expand_fast_hmm_param_if_necessary(ft, list_index + ft->last_state + ft->last_state + 1));
        }
        
        new_k = ft->last_state +1;
       
        
        list = ft->list;
        list_index = ft->num_items;
        sum = 0.0;
        for(i = 0;i <= ft->last_state;i++){
                list[list_index]->from = new_k;
                list[list_index]->to = i;
                if(i!= iHMM_START_STATE){
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
        

        //now split prob in last columns...
	
        a = alpha * beta[ft->last_state];
        b = 0.0;
        for(i = 0; i <= ft->last_state;i++){
                b += beta[i];
        }
        
        b = alpha * (1.0 - b);

        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_asc);

        l = binarySearch_to_lower_bound(ft,ft->last_state);
        r = binarySearch_to_upper_bound(ft,ft->last_state);

        for(i = l;i < r;i++){
                if(a < 1e-2 || b < 1e-2){     // % This is an approximation when a or b are really small.
                        pg = rk_binomial(&rndstate, 1.0, a / (a+b));
                }else{
                        pg = rk_beta(&rndstate, a, b);
                }
                pe = list[i]->t;
                list[i]->t = pg * pe;

                list[list_index]->from = list[i]->from;
                list[list_index]->to = new_k;
                list[list_index]->t = (1.0-pg) * pe;
               
                list_index++;
  
        }
        ft->num_items = list_index;
        ft->last_state = new_k;

        return OK;
ERROR:
        return FAIL;
}




struct fast_hmm_param* alloc_fast_hmm_param(int k, int L)
{
        struct fast_hmm_param* ft = NULL;
        int i;
        
        ASSERT(L > 1, "Need more than one letter");
        MMALLOC(ft, sizeof(struct fast_hmm_param));
        ft->alloc_num_states = k+2;
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
static int binarySearch_t(struct fast_hmm_param* ft, float x)
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



static int binarySearch_to_lower_bound(struct fast_hmm_param* ft, int x)
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

static int binarySearch_to_upper_bound(struct fast_hmm_param* ft, int x)
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

static int binarySearch_from_lower_bound(struct fast_hmm_param* ft, int x)
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

static int binarySearch_from_upper_bound(struct fast_hmm_param* ft, int x)
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






#ifdef ITEST


/* for testing..  */
static int fill_with_random_transitions(struct fast_hmm_param* ft, int k);
static int print_fast_hmm_params(struct fast_hmm_param* ft);


int main(const int argc,const char * argv[])
{
        fprintf(stdout,"Hello world\n");
        struct fast_hmm_param* ft = NULL;
        int i;
        int res = 0;
        float x; 

        RUNP(ft = alloc_fast_hmm_param(4,4));

        RUN(fill_with_random_transitions(ft, 8));
        RUN(print_fast_hmm_params(ft));
        fprintf(stdout,"%d items\n",ft->num_items);
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*), fast_hmm_param_cmp_by_t_desc);
        
        RUN(print_fast_hmm_params(ft));
        for(i =0; i < 10;i++){
                x = random_float_zero_to_x(1.0);
                
                res = binarySearch_t(ft, x);
                fprintf(stdout,"search for %f: %d  \n",x, res);
        }
        x = ft->list[3]->t;

        res = binarySearch_t(ft,x);
        fprintf(stdout,"search for %f: %d   \n",x, res);

        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_from_asc);
        RUN(print_fast_hmm_params(ft));


        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_to_asc);
        RUN(print_fast_hmm_params(ft));
        for(i =0; i < 4;i++){
                fprintf(stdout,"%d: %d -> %d\n",i, binarySearch_to_lower_bound(ft,i), binarySearch_to_upper_bound(ft,i));
        }
        
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*),fast_hmm_param_cmp_by_from_asc);
        RUN(print_fast_hmm_params(ft));
        for(i =0; i < 4;i++){
                fprintf(stdout,"%d: %d -> %d\n",i, binarySearch_from_lower_bound(ft,i), binarySearch_from_upper_bound(ft,i));
        }

        
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*), fast_hmm_param_cmp_by_to_from_asc);
        RUN(print_fast_hmm_params(ft));
        

         
        RUN(expand_fast_hmm_param_if_necessary(ft, 4));

        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*), fast_hmm_param_cmp_by_to_from_asc);
        
       

        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*), fast_hmm_param_cmp_by_to_from_asc);
        RUN(print_fast_hmm_params(ft));
        
        float* beta = NULL;
        float gamma = 6;
        float alpha = 1000.2;
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
                RUN(add_state_from_fast_hmm_param(rndstate,ft,  beta, alpha, gamma));
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

static int print_fast_hmm_params(struct fast_hmm_param* ft)
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
                        fprintf(stdout," %0.2f",m[i][j]);
                        sum+= m[i][j];
                }
                fprintf(stdout,"\ts:%f\n",sum);
        }
        fprintf(stdout,"\n");
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
        
        srand48(time(0));
        
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

#endif 
