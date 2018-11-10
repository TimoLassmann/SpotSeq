
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

// auxiliary functions for RB tree...

static void* get_transition(void* ptr)
{
        struct fast_t_item* tmp = (struct fast_t_item*)  ptr;
        return &tmp->t;
}

static int resolve_default(void* ptr_a,void* ptr_b)
{
        return 1;
}

long int compare_transition(void* keyA, void* keyB)
{
        float* num1 = (float*)keyA;
        float* num2 = (float*)keyB;

        if(*num1 == *num2){
                return 0;
        }
        if(*num1 > *num2){
                return 1;
        }else{
                return -1;
        }
}


void print_fast_t_item_struct(void* ptr,FILE* out_ptr)
{
        struct fast_t_item* tmp = (struct fast_t_item*)  ptr;
        fprintf(out_ptr,"%d\t%d\t%f\n",tmp->from , tmp->to,tmp->t);
}

void free_fast_t_item_struct(void* ptr)
{
        //struct fast_t_item* tmp = (struct fast_t_item*)  ptr;
        if(ptr){
                MFREE(ptr);
        }
}



struct fast_param_bag* alloc_fast_param_bag(int num_models, int k, int L)
{
        struct fast_param_bag* b = NULL;
        int i;
        ASSERT(num_models > 0, "No models");

        MMALLOC(b, sizeof(struct fast_param_bag));

        b->fast_params = NULL;
        b->num_models = num_models;

        MMALLOC(b->fast_params,sizeof(struct fast_hmm_param*)* b->num_models);

        for(i = 0; i < b->num_models;i++){
                b->fast_params[i] = NULL;
                RUNP(b->fast_params[i] = alloc_fast_hmm_param(k,L));
        }

        return b;
ERROR:
        return NULL;
}

void free_fast_param_bag(struct fast_param_bag* b)
{
        int i;
        if(b){
                if(b->fast_params){
                        for(i = 0; i < b->num_models;i++){
                                if(b->fast_params[i]){
                                        free_fast_hmm_param(b->fast_params[i]);
                                }
                        }
                        MFREE(b->fast_params);
                }
                MFREE(b);
        }

}

/* Goal: efficient datastructure for double indexing transition: */
/* 1) by from and to */
/* 2) by accessing a list sorted by transition probability  */
/* Plan: have rb tree for sorted list (i.e. insert then 'flatten' tree)  */
/* AND a normal transition matrix. */
/* For efficiency, let's add a vector of transitin structs to avoid millions of malloc calls */
/* Actually let's not as nodes in RB-tree are allocated anyhow.. */

struct fast_hmm_param* alloc_fast_hmm_param(int k, int L)
{
        struct fast_hmm_param* ft = NULL;
        int i;
        /* function pointers for left leaning RB tree..  */
        void*  (*fp_get)(void* ptr) = NULL;
        long int (*fp_cmp)(void* keyA, void* keyB)= NULL;
        int (*fp_cmp_same)(void* ptr_a,void* ptr_b);
        void (*fp_print)(void* ptr,FILE* out_ptr) = NULL;
        void (*fp_free)(void* ptr) = NULL;

        ASSERT(L > 1, "Need more than one letter");
        MMALLOC(ft, sizeof(struct fast_hmm_param));
        ft->alloc_num_states = k;
        ft->alloc_items = 1024 ;
        ft->last_state = 0;
        ft->list = NULL;
        ft->infinity = NULL;
        ft->num_items = 0;
        ft->emission = NULL;    /* This will be indexed by letter i.e. e['A']['numstate'] */
        ft->transition = NULL;

        ft->L = L;
        ft->background_emission = NULL;

        ft->root = NULL;

        MMALLOC(ft->background_emission, sizeof(float) * L  );

        for(i = 0; i < L; i++){
                ft->background_emission[i] = 0.0f;
        }

        RUNP(ft->emission = malloc_2d_float(ft->emission, ft->L, ft->alloc_num_states, 0.0f));
        RUNP(ft->transition = malloc_2d_float(ft->transition,  ft->alloc_num_states,  ft->alloc_num_states, 0.0f));

        fp_get = &get_transition;
        fp_cmp = &compare_transition;
        fp_print = &print_fast_t_item_struct;
        fp_cmp_same = &resolve_default;
        fp_free = &free_fast_t_item_struct;

        ft->root = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free);

        MMALLOC(ft->infinity, sizeof(struct fast_t_item*) * ft->alloc_num_states);

        for(i = 0; i < ft->alloc_num_states;i++){
                ft->infinity[i] = NULL;
                MMALLOC(ft->infinity[i], sizeof(struct fast_t_item));
                ft->infinity[i]->from = -1;
                ft->infinity[i]->to = -1;
                ft->infinity[i]->t = 0.0f;
        }
        return ft;
ERROR:
        free_fast_hmm_param(ft);
        return NULL;
}

int expand_ft_if_necessary(struct fast_hmm_param* ft, int new_num_states)
{
        int i, num_old_item;
        ASSERT(ft != NULL, "No ft struct!");
        ASSERT(new_num_states >2,"No states requested");

        if(new_num_states > ft->alloc_num_states){
                num_old_item = ft->alloc_num_states;
                while(new_num_states > ft->alloc_num_states){
                        ft->alloc_num_states = ft->alloc_num_states + 64;
                }
                RUNP(ft->emission = malloc_2d_float(ft->emission, ft->L, ft->alloc_num_states, 0.0f));
                RUNP(ft->transition = malloc_2d_float(ft->transition,  ft->alloc_num_states,  ft->alloc_num_states, 0.0f));

                MREALLOC(ft->infinity, sizeof(struct fast_t_item*) * ft->alloc_num_states);
                for(i = num_old_item; i < ft->alloc_num_states;i++){
                        ft->infinity[i] = NULL;
                        MMALLOC(ft->infinity[i], sizeof(struct fast_t_item));
                        ft->infinity[i]->from = -1;
                        ft->infinity[i]->to = -1;
                        ft->infinity[i]->t = 0.0f;
                }
        }
        return OK;
ERROR:
        free_fast_hmm_param(ft);
        return FAIL;
}

/*int expand_transition_if_necessary(struct fast_hmm_param* ft)
  {
  int i, num_old_item;
  ASSERT(ft != NULL, "No ft struct!");
  num_old_item = ft->alloc_items;
  ft->alloc_items += 1024;
  LOG_MSG("expanding from: %d to %d",num_old_item, ft->alloc_items);
  MREALLOC(ft->infinity, sizeof(struct fast_t_item*) * ft->alloc_items);
  for(i = num_old_item; i < ft->alloc_items;i++){
  ft->infinity[i] = NULL;
  MMALLOC(ft->infinity[i], sizeof(struct fast_t_item));
  ft->infinity[i]->from = -1;
  ft->infinity[i]->to = -1;
  ft->infinity[i]->t = 0.0f;
  }
  return OK;
  ERROR:
  free_fast_hmm_param(ft);
  return FAIL;
  }*/

void free_fast_hmm_param(struct fast_hmm_param* ft)
{
        int i;
        if(ft){
                if(ft->root){
                        ft->root->free_tree(ft->root);
                }
                if(ft->infinity){
                        for(i = 0; i < ft->alloc_num_states;i++){
                                MFREE(ft->infinity[i]);
                        }
                        MFREE(ft->infinity);
                }
                if(ft->emission){
                        free_2d((void**) ft->emission);
                }

                if(ft->transition){
                        free_2d((void**) ft->transition);
                }
                if(ft->background_emission){
                        MFREE(ft->background_emission);
                }
                MFREE(ft);
        }
}


int make_flat_param_list(struct fast_hmm_param* ft)
{
        ASSERT(ft != NULL, "No parameters");
        if(ft->root->data_nodes){
                MFREE(ft->root->data_nodes);
                ft->root->data_nodes = NULL;
                ft->root->cur_data_nodes = 0;
        }
        RUN(ft->root->flatten_tree(ft->root));

        ft->num_items = ft->root->num_entries;
        ft->list = (struct fast_t_item**) ft->root->data_nodes;
        return OK;
ERROR:
        return FAIL;
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

        struct fast_hmm_param* ft = NULL;
        struct fast_t_item* tmp = NULL;
        int i,j;
        int res = 0;
        float x;

        RUN(print_program_header((char * const*)argv,"Fast HMM param test"));

        LOG_MSG("Allocating fast transition data structure.");
        RUNP(ft = alloc_fast_hmm_param(2,4));

        LOG_MSG("Fill with random transition (this will also expand the data structure if there isn't enough space..)");
        RUN(fill_with_random_transitions(ft, 4));


        RUN(print_fast_hmm_params(ft));

        RUN(make_flat_param_list(ft) );

        //ft->root->flatten_tree(ft->root);

        //ft->root->print_tree(ft->root,NULL);
        LOG_MSG("Sorted by RB tree.");
        fprintf(stdout,"%d items\n",ft->root->num_entries);
        for(i = 0; i < ft->root->num_entries;i++){
                tmp = (struct fast_t_item*) ft->root->data_nodes[i];
                fprintf(stdout,"%d %f\n",i , tmp->t);
        }
        LOG_MSG("Sorted by sort function.");
        ft->num_items = ft->root->num_entries;
        ft->list = (struct fast_t_item**) ft->root->data_nodes;
        qsort(ft->list, ft->num_items, sizeof(struct fast_t_item*), fast_hmm_param_cmp_by_t_desc);
        for(i = 0; i < ft->num_items;i++){
                fprintf(stdout,"%d %f\n",i , ft->list[i]->t);
        }

        RUN(print_fast_hmm_params(ft));

        for(i =0; i < 10;i++){
                x = random_float_zero_to_x(1.0);

                res = fast_hmm_param_binarySearch_t(ft, x);
                fprintf(stdout,"search for %f: %d  \n",x, res);
                for(j = 0; j < res;j++){
                        ASSERT(ft->list[j]->t >= x,"Warning - binary search seems to have failed");
                }
        }
        x = ft->list[3]->t;

        res = fast_hmm_param_binarySearch_t(ft,x);
        fprintf(stdout,"search for %f: %d   \n",x, res);
        for(j = 0; j < res;j++){
                ASSERT(ft->list[j]->t >= x,"Warning - binary search seems to have failed");
        }
        RUN(print_fast_hmm_params(ft));

        free_fast_hmm_param(ft);
        return EXIT_SUCCESS;
ERROR:

        free_fast_hmm_param(ft);
        return EXIT_FAILURE;
}


#endif
