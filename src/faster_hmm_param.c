#include "faster_hmm_param.h"

/* compare functions */
static int compare_transition (const void *pa, const void *pb, void *param);
static int compare_ab(const void *pa, const void *pb, void *param);
/* alloc functions... */
/* use my malloc wrapper in libavl ...  */
void* rb_malloc_tldevel (struct libavl_allocator *allocator, size_t size)
{
        void* p = NULL;
        ASSERT(allocator != NULL, "No allocator");

        MMALLOC(p,size);
        return p;
ERROR:
        return NULL;
}

void rb_free_tldevel (struct libavl_allocator *allocator, void *block)
{
        MFREE(block);
}

struct libavl_allocator rb_allocator_tldevel =
{
        rb_malloc_tldevel,
        rb_free_tldevel
};


static struct faster_hmm_param* alloc_faster_hmm_param(int k, int L);
static void free_faster_hmm_param(struct faster_hmm_param* ft);
static void free_faster_transition_node(void* item,void *param);
/* alloc functions...  done */


/* exploring new states  */
int get_max_to_infty_transition(struct faster_hmm_param*ft,double* max);
int add_state(struct faster_hmm_param* ft,int* num_states,double* beta, double alpha, double gamma, rk_state* rndstate);

/* setting u */

int sample_u_for_each_residue(struct faster_hmm_param*ft, int* label, int len)
{
        int i,j;
        double x;
        was here - need to get random number for each transition based on previous label, insert into u tree.
        x = ft->transition[IHMM_START_STATE][label[0]];
        return OK;
}


int purge_rb_trees(struct faster_hmm_param*ft)
{
        ASSERT(ft != NULL, "No fast parameters");
        rb_destroy(ft->root,  free_faster_transition_node);
        rb_destroy(ft->u,  free_faster_transition_node);
        rb_destroy(ft->boundary,  free_faster_transition_node);

        ft->root = NULL;
        ft->u = NULL;
        ft->boundary = NULL;

        RUNP(ft->root = rb_create(compare_transition, NULL,&rb_allocator_tldevel));
        RUNP(ft->u = rb_create(compare_transition, NULL,&rb_allocator_tldevel));
        RUNP(ft->boundary = rb_create(compare_ab, NULL,&rb_allocator_tldevel));
        return OK;
ERROR:
        return FAIL;
}


struct faster_param_bag* alloc_faster_param_bag(int num_models, int* K, int L)
{
        struct faster_param_bag* b = NULL;
        int i;
        ASSERT(num_models > 0, "No models");

        MMALLOC(b, sizeof(struct faster_param_bag));

        b->fast_params = NULL;
        b->num_models = num_models;
        b->max_last_state = -1;

        MMALLOC(b->fast_params,sizeof(struct faster_hmm_param*)* b->num_models);

        for(i = 0; i < b->num_models;i++){
                b->fast_params[i] = NULL;
                RUNP(b->fast_params[i] = alloc_faster_hmm_param(K[i],L));
        }

        return b;
ERROR:
        return NULL;
}

void free_faster_param_bag(struct faster_param_bag* b)
{
        int i;
        if(b){
                if(b->fast_params){
                        for(i = 0; i < b->num_models;i++){
                                if(b->fast_params[i]){
                                        free_faster_hmm_param(b->fast_params[i]);
                                }
                        }
                        MFREE(b->fast_params);
                }
                MFREE(b);
        }

}




struct faster_hmm_param* alloc_faster_hmm_param(int k, int L)
{
        struct faster_hmm_param* ft = NULL;
        int i;
        ASSERT(L > 1, "Need more than one letter");
        MMALLOC(ft, sizeof(struct faster_hmm_param));
        ft->alloc_num_states = k;
        ft->alloc_items = 1024 ;
        ft->last_state = 0;
        ft->infinity = NULL;
        ft->num_items = 0;
        ft->emission = NULL;    /* This will be indexed by letter i.e. e['A']['numstate'] */
        ft->transition = NULL;

        ft->L = L;
        ft->background_emission = NULL;

        ft->root = NULL;
        RUNP(ft->root = rb_create(compare_transition, NULL,&rb_allocator_tldevel));
         RUNP(ft->u = rb_create(compare_transition, NULL,&rb_allocator_tldevel));
        RUNP(ft->boundary = rb_create(compare_ab, NULL,&rb_allocator_tldevel));

        MMALLOC(ft->background_emission, sizeof(double) * L  );

        for(i = 0; i < L; i++){
                ft->background_emission[i] = 0.0;
        }

        //RUNP(ft->emission = galloc(ft->emission, ft->L, ft->alloc_num_states, 0.0));
        //RUNP(ft->transition = galloc(ft->transition,  ft->alloc_num_states,  ft->alloc_num_states, 0.0));

        MMALLOC(ft->infinity, sizeof(struct faster_t_item*) * ft->alloc_num_states);

        for(i = 0; i < ft->alloc_num_states;i++){
                ft->infinity[i] = NULL;
                MMALLOC(ft->infinity[i], sizeof(struct faster_t_item));
                ft->infinity[i]->a = -1;
                ft->infinity[i]->b = -1;
                ft->infinity[i]->t = 0.0;
        }
        return ft;
ERROR:
        free_faster_hmm_param(ft);
        return NULL;
}


/* Master function to explore new states  */

int explore_new_states(struct faster_hmm_param*ft, double min_u,int* num_states,double* beta, double alpha, double gamma, rk_state* rndstate)
{
        double max;

        RUN( get_max_to_infty_transition(ft, &max));
        while(max >= min_u && *num_states < MAX_NUM_STATES && max > 0.0 ){//}sb->max_len){
                RUN(add_state(ft, num_states,beta,alpha,gamma,rndstate));
                RUN( get_max_to_infty_transition(ft, &max));

        }
        return OK;
ERROR:
        return FAIL;
}

/* Function to add states  */

int get_max_to_infty_transition(struct faster_hmm_param*ft,double* max)
{

        int i;
        double local_max;

        ASSERT(ft != NULL, "No fast hmm parameters.");

        local_max = -1.0;


        for(i = 0; i< ft->last_state;i++){
                if(ft->infinity[i]->t > local_max){
                        local_max = ft->infinity[i]->t;
                }
                //fprintf(stdout,"%d->%d %f\n", ft->infinity[i]->from, ft->infinity[i]->to, ft->infinity[i]->t);
        }
        *max = local_max;
        return OK;
ERROR:
        return FAIL;
}

int add_state(struct faster_hmm_param* ft,int* num_states,double* beta, double alpha, double gamma, rk_state* rndstate)
{
        struct faster_t_item** infinity = NULL;
        struct faster_t_item* data = NULL;
        double* tmp_prob = NULL;
        //rk_state rndstate;

        double sum,be,bg,pe,pg, a,b;
        int i,new_k;


        ASSERT(ft != NULL, "No ft.");

        *num_states = *num_states +1;
        //RUN(expand_ft_if_necessary(ft, ihmm->num_states));
        /* +1 for infty + 1 for new state.. */
        MMALLOC(tmp_prob, sizeof(double) *(*num_states));


        new_k = ft->last_state;
        infinity = ft->infinity;
        //fprintf(stdout,"LAST: %d\n",new_k);
        /* fill out transition FROM new state  */
        sum = 0.0;
        for(i = 0;i <= new_k;i++){
                tmp_prob[i] =  rk_gamma(rndstate, beta[i] * alpha, 1.0);
                if(i == IHMM_START_STATE){
                        tmp_prob[i] = 0.0;
                }
                sum += tmp_prob[i];
        }
        for(i = 0;i < new_k;i++){
                data = NULL;
                MMALLOC(data, sizeof(struct faster_t_item));
                data->a = new_k;
                data->b = i;
                data->t = tmp_prob[i] / sum;
                data = rb_insert(ft->root, data);
                if(data){
                        ERROR_MSG("Insert has failed");
                }
        }
        infinity[new_k]->a = new_k;
        infinity[new_k]->b = new_k;
        infinity[new_k]->t = tmp_prob[new_k] / sum;

        //first get beta for new column
        be = beta[new_k];
        bg = rk_beta(rndstate, 1.0,gamma );

        beta[new_k] = bg*be;
        beta[new_k+1] = (1.0 - bg) *be;

        //ihmm->beta = beta;
        //now split prob in last columns...
        a = alpha * beta[new_k];
        b = 0.0;
        for(i = 0; i <= new_k;i++){
                b += beta[i];
        }
        b = alpha * (1.0 - b);

        // split last column - i.e. play with infinity.

        for(i = 0 ; i <= new_k;i++){
                if(a < 1e-2 || b < 1e-2){     // % This is an approximation when a or b are really small.
                        pg = rk_binomial(rndstate, 1.0, a / (a+b));
                }else{
                        pg = rk_beta(rndstate, a, b);
                }
                pe = infinity[i]->t;


                //transition to state just instantiated will go into the RB tree.
                data = NULL;
                MMALLOC(data, sizeof(struct faster_t_item));
                data->a = i;
                data->b =  new_k;
                data->t = pg * pe;
                data = rb_insert(ft->root, data);
                if(data){
                        ERROR_MSG("Insert has failed");
                }

                //transition into infinity will remain in the infinity array...
                infinity[i]->a = i;
                infinity[i]->b = new_k+1;
                infinity[i]->t = (1.0-pg) * pe;
        }


        /* add emission  */
        sum = 0.0;
        for(i = 0; i < ft->L;i++){
                ft->emission[i][new_k] = rk_gamma(rndstate, ft->background_emission[i], 1.0);
                sum += ft->emission[i][new_k];
        }
        for(i = 0; i < ft->L;i++){
                ft->emission[i][new_k] /= sum;
        }


        //MFREE(tmp_pg);
        //ft->num_items = list_index;
        ft->last_state = new_k+1;
        //ihmm->rndstate = rndstate;
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





int add_emission_to_frhmmp(struct faster_hmm_param* ft, double** e, int states)
{
        ASSERT(e != NULL, "No transitions");
        ASSERT(ft != NULL, "no faster parameters");
        ASSERT(states > 2, "too few states");
        ft->emission = e;

        ft->last_state = states - 1; /* is redundant because last_state is also set in add_transitions. */

        return OK;
ERROR:
        return FAIL;
}

int add_transitions_to_frhmmp(struct faster_hmm_param* ft, double** t, int states)
{
        struct faster_t_item** infinity = NULL;
        struct faster_t_item* data = NULL;
        int last_state;
        int i,j;


        ASSERT(t != NULL, "No transitions");
        ASSERT(ft != NULL, "no faster parameters");
        ASSERT(states > 2, "too few states");

        infinity = ft->infinity;
        last_state = states -1;
        data = NULL;
        MMALLOC(data, sizeof(struct faster_t_item));
        data->a = IHMM_START_STATE;
        data->b = IHMM_START_STATE;
        data->t = 0.0;
        data = rb_insert(ft->root , data);
        if(data){
                ERROR_MSG("Insert failed");
        }

        data = NULL;
        MMALLOC(data, sizeof(struct faster_t_item));
        data->a = IHMM_START_STATE;
        data->b = IHMM_END_STATE;
        data->t = 0.0;
        data = rb_insert(ft->root , data);
        if(data){
                ERROR_MSG("Insert failed");
        }

        /* transitions from start state... */
        for(i = 2; i < last_state;i++){
                data = NULL;
                MMALLOC(data, sizeof(struct faster_t_item));
                data->a = IHMM_START_STATE;
                data->b = i;
                data->t = t[IHMM_START_STATE][i];
                data = rb_insert(ft->root , data);
                if(data){
                        ERROR_MSG("Insert failed");
                }

        }
        infinity[IHMM_START_STATE]->a = IHMM_START_STATE;
        infinity[IHMM_START_STATE]->b = last_state;
        infinity[IHMM_START_STATE]->t = t[IHMM_START_STATE][last_state];

        /* transitions from stop state... */
        for(i = 0; i < last_state;i++){
                data = NULL;
                MMALLOC(data, sizeof(struct faster_t_item));
                data->a = IHMM_END_STATE;
                data->b = i;
                data->t = 0.0;

                data = rb_insert(ft->root , data);
                if(data){
                        ERROR_MSG("Insert failed");
                }
        }
        infinity[IHMM_END_STATE]->a = IHMM_END_STATE;
        infinity[IHMM_END_STATE]->b = last_state;
        infinity[IHMM_END_STATE]->t = 0.0;

        /* the rest... */

        for(i = 2; i < last_state;i++){

                for(j = 1; j < last_state;j++){
                        data = NULL;
                        MMALLOC(data, sizeof(struct faster_t_item));
                        data->a = i;
                        data->b = j;
                        data->t = t[i][j];
                        data = rb_insert(ft->root , data);
                        if(data){
                                ERROR_MSG("Insert failed");
                        }

                }
                infinity[i]->a = i;
                infinity[i]->b = last_state;
                infinity[i]->t =  t[i][last_state];
        }

        ft->last_state = last_state;
        return OK;
ERROR:
        return FAIL;
}

void free_faster_hmm_param(struct faster_hmm_param* ft)
{
        int i;
        if(ft){
                if(ft->root){
                        rb_destroy(ft->root,  free_faster_transition_node);
                }
                if(ft->infinity){
                        for(i = 0; i < ft->alloc_num_states;i++){
                                MFREE(ft->infinity[i]);
                        }
                        MFREE(ft->infinity);
                }
                if(ft->emission){
                        gfree(ft->emission);
                }

                if(ft->transition){
                        gfree(ft->transition);
                }
                if(ft->background_emission){
                        MFREE(ft->background_emission);
                }
                MFREE(ft);
        }
}

int compare_ab(const void *pa, const void *pb, void *param)
{
        const struct faster_t_item *a = pa;
        const struct faster_t_item *b = pb;
        int fa,fb;


        fa = a->a;      /* from a node name  */
        fb = b->a;       /* from b node name  */
        if(fa < fb){
                return -1;
        }else if(fa > fb){
                return 1;
        }else{ /* if transition probability and from index is equal   compare to index */
                int ta,tb;
                ta = a->b;
                tb = b->b;
                if(ta < tb){
                        return -1;
                }else if(ta > tb){
                        return 1;
                }else{
                        LOG_MSG("a == b");
                        return 0;
                }
        }
}

int compare_transition (const void *pa, const void *pb, void *param)
{
        const struct faster_t_item *a = pa;
        const struct faster_t_item *b = pb;

        double ta,tb;
        ta = a->t;
        tb = b->t;

        if(ta < tb){
                return -1;
        }else if(ta > tb){
                return 1;
        }else{                  /* if transition probability is equal compare from index */
                int fa,fb;


                fa = a->a;      /* from a node name  */
                fb = b->a;       /* from b node name  */
                if(fa < fb){
                        return -1;
                }else if(fa > fb){
                        return 1;
                }else{ /* if transition probability and from index is equal   compare to index */
                        int ta,tb;
                        ta = a->b;
                        tb = b->b;
                        if(ta < tb){
                                return -1;
                        }else if(ta > tb){
                                return 1;
                        }else{
                                LOG_MSG("a == b");
                                return 0;
                        }
                }
        }
}

void free_faster_transition_node(void* item,void *param)
{
        MFREE( item);
}



int main (int argc,char * argv[])
{
        fprintf(stdout,"Hello world\n");
        struct faster_param_bag*fb = NULL;
        int* start_states = NULL;
        int num_start_states = 10;
        int i;


        MMALLOC(start_states, sizeof(int) * num_start_states);

        for(i = 0; i < num_start_states;i++){
                start_states[i] = i + 10;

        }
        RUNP(fb = alloc_faster_param_bag(num_start_states,start_states   ,4));


        free_faster_param_bag(fb);

        MFREE(start_states);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;

}
