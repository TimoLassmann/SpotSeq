
#include "tldevel.h"
#include "tllogsum.h"
#include <string.h>

#include "null_model_emission.h"

#include "pst.h"


struct pst_node{
        struct pst_node* next[20];
        float probability[20];
        char* label;
        int label_len;

};

#include "pst_structs.h"

static int init_pst(struct pst** pst,float expected_error, int len, int L);


struct pst_node* build_pst(struct pst* pst,struct fpst*f,int curf, struct pst_node* n, struct count_hash* h);
static void print_pst(struct pst* pst,struct pst_node* n);
static void free_pst_node(struct pst_node* n);

static int prob2scaledprob_fpst(struct fpst* f);

static int alloc_node(struct pst_node** node,uint8_t* string,int len, int L);
static int alloc_fpst(struct fpst** fast_pst, int size, int L);
static int resize_fpst(struct fpst* f, int L);
static void free_fpst(struct fpst* f);

int score_pst(const struct pst* pst, const char* seq,const int len, float* P_M, float* P_R)
{
        int** s = pst->fpst_root->links;
        float** s_prob = pst->fpst_root->prob;
        const float* base_p = pst->fpst_root->prob[0];
        float a, b;
        register int i,c,l,pos,n;
        a = prob2scaledprob(1.0);
        b = a;
        for(i = 0; i < len; i++){
                l = seq[i];
                b = b + base_p[l];
                pos = i;
                n = 0;
                while(pos){
                        //LOG_MSG("Looking at pos: %d character: %d",pos, seq[pos]);
                        pos= pos-1;
                        c = seq[pos];
                        if(!s[n][c]){
                                break;
                        }
                        n = s[n][c];
                }
                a += s_prob[n][l];
        }
        *P_M = a;
        *P_R = b;
        return OK;
}


int run_build_pst(struct pst** pst, float expected_error, struct count_hash* h)//  struct tl_seq_buffer* sb)
{
        struct pst* p = NULL;
        struct pst_node* helper = NULL;
        float sum;
        int i;
        float c;
        uint64_t x;
        uint64_t key;
        khiter_t hi;

        init_logsum();

        RUN(init_pst(&p,expected_error,  h->len,h->L));

        p->p_min = MACRO_MAX(h->min_prop, p->p_min);

        RUN(alloc_node(&helper,NULL,0,p->L));

        sum = 0.0;
        for(i = 0;i < p->L;i++){
                key =  (1ULL<< 60ULL) | (uint64_t)i;
                hi= kh_get(exact, h->hash, key);
                c = 0;
                if(hi != kh_end(h->hash)){
                        //LOG_MSG("%d: %d", i,kh_value(h->hash, hi));
                        c = kh_value(h->hash, hi);
                }

                x = p->fpst_root->l;
                p->fpst_root->prob[x][i] = c;
                //p->fppt_root->prob[x][i] = c;
                sum += p->fpst_root->prob[x][i];
        }


        for(i = 0;i < p->L;i++){
                x = p->fpst_root->l;
                p->fpst_root->prob[x][i] /=sum;
                p->fpst_root->prob[x][i] = p->fpst_root->prob[x][i] * ( 1.0 -  p->gamma_min) + p->background[i]* p->gamma_min;
                helper->probability[i] = p->fpst_root->prob[x][i];
                p->fpst_root->links[x][i] = 0;
        }
        p->fpst_root->l = 0;

        helper = build_pst(p,p->fpst_root,0, helper,h);
        print_pst(p,helper);
        p->fpst_root->l++;
        free_pst_node(helper);
        helper = NULL;
        /*RUN(init_helper(&helper, k->counts[0], p->gamma_min));
        helper = build_ppt(p,p->fppt_root,0, helper,k);
        p->fppt_root->l++;
        free_pst_node(helper);
        helper = NULL;*/
        RUN(prob2scaledprob_fpst(p->fpst_root));
        //RUN(prob2scaledprob_fpst(p->fppt_root));
        //fprintf(stdout,"Size: %d\n", p->fppt_root->l + p->fpst_root->l * 8 *4);
//LOG_MSG("PPT:");
        //exit(0);
        *pst = p;
        return OK;
ERROR:
        return FAIL;
}


struct pst_node* build_pst(struct pst* pst,struct fpst*f,int curf, struct pst_node* n, struct count_hash* h)
{
        khiter_t hi;
        uint8_t ask[MAX_PST_LEN+20];

        int i;
        int j;
        int c;
        //int add;

        int len = n->label_len;
        float sum = 0.0f;
        float tmp_counts_s[20];
        float P_ask;
        float P_as;
        float Err;

        uint64_t ask_code;

        uint64_t key;

        int maxlen;

        maxlen = h->len;
        //step 2 test expansion
        //loop though letters at present node

        //LOG_MSG("%d < %d", len+1, maxlen);
        if(len + 1 < maxlen ){
                for(i = 0; i < pst->L;i++){

                        /* Let's calculate whether extending the suffix by letter A (indexed by i)
                           (i.e. adding it to the left of the existing label)
                           adds statistically significantly information to the PST.

p                           What I need is the probability of:
                           1) - P(ASK) probability of observing string ASK;
                           where 'A' is the letter by which I extend the suffix,
                           S is the existing suffix and K is any letter in the alphabet L
                           2) - P(SK) probability of observing string SK
                           3) - P(K | S) - probability of observing K after S

                           Stastistical error = (sum over all K) P(ASK) * log( P(ASK) / ( P(K|S)  * P (AS))

                        */
                        // here I start constructing the ASK by adding A and then copying the existing suffix S
                        ask_code = (uint64_t)i;
                        ask[0] = i;
                        for(j = 1; j < len+1;j++){
                                ask_code = (ask_code << 5ULL) | (uint64_t)n->label[j-1];
                                ask[j] = n->label[j-1];
                        }
                        c = 0;
                        key = ((uint64_t)(len+1) << 60ULL) | ask_code;
                        hi= kh_get(exact, h->hash, key);
                        if(hi != kh_end(h->hash)){
                                c  = kh_value(h->hash, hi);
                        }

                        P_as = (float) (c+1) / (float)(h->counts_l[len+1] + pst->L);

                        if(P_as > pst->p_min){
                                Err = 0.0f;
                                sum = 0.0f;
                                //to calculate the probabilities I add K
                                for(j = 0; j < pst->L;j++){
                                        ask_code = ask_code << 5ULL | (uint64_t)j;
                                        ask[len+1] = j;
                                        ask[len+2] = 0;
                                        /* now the ask code is complete - I can start looking up their counts to get to P(ask) and P(as) */
                                        c = 0;
                                        key = ((uint64_t)(len+2) << 60ULL) | ask_code;
                                        hi= kh_get(exact, h->hash, key);
                                        if(hi != kh_end(h->hash)){
                                                //LOG_MSG("%c: %f", "ACGT"[i],kh_value(h->hash, hi));
                                                c  = kh_value(h->hash, hi);
                                        }
                                        tmp_counts_s[j] = c;
                                        sum += c;

                                        P_ask =   (float) (c+1) / (float)(h->counts_l[len+2] + pst->L);

                                        ask_code = ask_code >> 5ULL;

                                        Err += P_ask * logf(P_ask /  ( n->probability[j] * P_as));

                                }
                                /* Shall we add this suffix ?  */
                                if (Err > pst->p_min){
                                        for(j = 0; j < pst->L;j++){
                                                tmp_counts_s[j] = tmp_counts_s[j] / sum;
                                                tmp_counts_s[j] = tmp_counts_s[j] * ( 1.0 -  pst->gamma_min) + pst->background[j]* pst->gamma_min;
                                        }
                                        f->links[curf][i] = f->l+1;

                                        f->l++;
                                        if(f->l== f->m){
                                                resize_fpst(f,pst->L);
                                        }
                                        RUN(alloc_node(&n->next[i] ,ask,len+1,pst->L));
                                        for(j = 0; j < pst->L;j++){
                                                f->prob[f->l][j] = tmp_counts_s[j];
                                                f->links[f->l][j] = 0;
                                                n->next[i]->probability[j] = tmp_counts_s[j];
                                        }
                                        n->next[i] = build_pst(pst,f, f->links[curf][i], n->next[i] ,h );
                                }
                        }
                }
        }
        return n;
ERROR:
        return NULL;
}

int init_pst(struct pst** pst,float expected_error, int len, int L)//, struct tl_seq_buffer* sb)
{
        struct pst* p = NULL;
        MMALLOC(p, sizeof(struct pst));
        p->len = len;
        p->L = L;
        p->gamma_min = 0.02f;
        if(expected_error){
                p->gamma_min = expected_error;
        }
        ASSERT(p->gamma_min * L <= 1.0,"expected error is too high.");
        p->p_min = 0.0001f;
        //p->r = 0.0001f;
        p->background = NULL;
        p->fpst_root = NULL;

        RUN(alloc_fpst(&p->fpst_root, 64,p->L));

        RUN(get_null_model_emissions(&p->background, p->L));
        *pst = p;
        return OK;
ERROR:
        free(p);
        return FAIL;
}

void free_pst(struct pst* p)
{
        if(p){
                free_fpst(p->fpst_root);
                gfree(p->background);
                MFREE(p);
        }
}

void free_pst_node(struct pst_node* n)
{
        int i;
        if(n){
                for(i = 0;i < 4;i++){
                        if(n->next[i]){
                                free_pst_node(n->next[i]);
                        }
                }
                MFREE(n->label);
                MFREE(n);
        }
}


int alloc_node(struct pst_node** node,uint8_t* string,int len, int L)
{
        struct pst_node* n =  NULL;
        int i;
        MMALLOC(n, sizeof(struct pst_node));
        n->label = NULL;
        MMALLOC(n->label, sizeof(char) * (len+1));
        //LOG_MSG("adding: ");
        for(i = 0; i < len;i++){
                n->label[i] = string[i];
                //fprintf(stdout,"%d ", n->label[i]);

        }
        //fprintf(stdout,"\n");
        n->label[len] = 0;
        n->label_len = len;
        for(i =0; i < L;i++){
                n->next[i] = NULL;
                n->probability[i] = 1.0f / (float) L;
        }
        *node = n;
        return OK;
ERROR:
        return FAIL;
}


void print_pst(struct pst* pst,struct pst_node* n)
{
        int i;
        int c = 0;

        for(i = 0;i < pst->L;i++){
                if(n->next[i]){
                        c++;
                }
        }
        if(!c){
                if(pst->L == 4){
                for(i = 0; i < n->label_len;i++){
                        fprintf(stdout,"%c","ACGT"[(int) n->label[i]]);
                }
                }else{
                        for(i = 0; i < n->label_len;i++){
                                fprintf(stdout,"%c","ACDEFGHIKLMNPQRSTVWY"[(int) n->label[i]]);
                        }


                }
                fprintf(stdout,"\t");
                for(i = 0;i < pst->L;i++){
                        fprintf(stdout,"%f ",n->probability[i]);
                }
                fprintf(stdout,"\n");
        }

        for(i = 0;i < pst->L;i++){
                if(n->next[i]){
                        //if(n->next[i]->in_T){
                        //fprintf(stderr,"Going:%d\n",i);
                        print_pst(pst,n->next[i]);
                        //}
                }
        }
}


int alloc_fpst(struct fpst** fast_pst, int size, int L)
{
        struct fpst*f = NULL;

        ASSERT(L > 3, "Strange alphabet...");
        ASSERT(size > 0, "no space!");
        MMALLOC(f, sizeof(struct fpst));
        f->prob = NULL;
        f->links = NULL;
        f->m = size;
        f->l= 0;

        RUN(galloc(&f->prob,f->m,L));
        RUN(galloc(&f->links,f->m,L));

        *fast_pst = f;
        return OK;
ERROR:
        return FAIL;
}

int resize_fpst(struct fpst* f, int L)
{
        int i,j, old;

        old = f->m;
        f->m = f->m + f->m /2;
        RUN(galloc(&f->prob,f->m,L));
        RUN(galloc(&f->links,f->m,L));

        for(i = old;i < f->m;i++){
                for(j = 0;j < L;j++){
                        f->prob[i][j] = 0.0f;
                        f->links[i][j] = 0;
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int prob2scaledprob_fpst(struct fpst* f)
{
        int i,j;

        ASSERT(f != NULL,"No F");

        for(i = 0; i < f->l;i++){
                for (j = 0; j < 4; j++) {
                        f->prob[i][j] = prob2scaledprob(f->prob[i][j]);

                }
        }
        /*fprintf(stdout,"\n");
        fprintf(stdout,"l: %d m: %d\n", f->l, f->m);
        i = f->l-1;
        for(j = 0; j < 4;j++){
                fprintf(stdout,"%f ",f->prob[i][j]);

        }
        fprintf(stdout,"\n");

        i = f->l;
        for(j = 0; j < 4;j++){
                fprintf(stdout,"%f ",f->prob[i][j]);

        }
        fprintf(stdout,"\n");
        */

        return OK;
ERROR:
        return FAIL;
}


void free_fpst(struct fpst* f)
{
        if(f){
                gfree(f->prob);
                gfree(f->links);
                MFREE(f);
        }
}


