
#include "tldevel.h"
#include "tllogsum.h"
#include <string.h>



#include "pst.h"


struct pst_node{
        struct pst_node* next[20];
        float probability[20];
        char* label;
        int label_len;

};

#include "pst_structs.h"
struct kmer_counts{
        uint32_t** counts;
        uint32_t* mask;
        int L;
};

//static int init_pst(struct pst** pst, int len);

static int init_pst(struct pst** pst,float expected_error, int len, int L);
static int init_helper(struct pst_node** helper, struct count_hash* ch,float gamma_min);
//static int init_helper(struct pst_node** helper, uint32_t* counts,float gamma_min);
static float get_pst_prob(struct pst_node* n, char* string,int target, int pos);

static float get_ppt_prob(struct pst_node* n, char* string,int target, int pos);

static float get_fpst_prob(struct fpst* f, char* string,int target, int pos);
static float get_fppt_prob(struct fpst* f, char* string,int target, int pos, int len);
//static struct pst_node* build_pst(struct pst* pst,struct pst_node* n );
//static struct pst_node* build_pst(struct pst* pst, struct fpst* f, struct pst_node* n,struct kmer_counts* k);
//struct pst_node* build_pst(struct pst* pst,struct fpst*f,int curf, struct pst_node* n,struct kmer_counts* k);

struct pst_node* build_pst(struct pst* pst,struct fpst*f,int curf, struct pst_node* n,struct kmer_counts* k, struct count_hash* h);
struct pst_node* build_ppt(struct pst* pst,struct fpst*f,int curf,struct pst_node* n ,struct kmer_counts* k);
//static struct pst_node* build_ppt(struct pst* pst,struct pst_node* n ,struct kmer_counts* k);
//static struct pst_node* build_ppt(struct pst* pst,struct pst_node* n );
static void print_pst(struct pst* pst,struct pst_node* n);

//static void print_pst(struct pst* pst,struct pst_node* n,int* num);
static void free_pst_node(struct pst_node* n);

//static struct pst_node* alloc_node(struct pst_node* n,char* string,int len);

static int alloc_node(struct pst_node** node,uint8_t* string,int len, int L);
static inline int nuc_to_internal(const char c);

static int alloc_fpst(struct fpst** fast_pst, int size, int L);
//static int alloc_fpst(struct fpst** fast_pst, int size);

static int resize_fpst(struct fpst* f, int L);
static int prob2scaledprob_fpst(struct fpst* f);
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
                l = nuc_to_internal(seq[i]);
                b = b + base_p[l];
                pos = i;
                n = 0;
                while(pos){
                        pos= pos-1;
                        c = nuc_to_internal(seq[pos]);
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

int score_pst_random(char* seq, int len,float* base_p, float*r)
{
        float P_R;
        register int i;
        P_R = prob2scaledprob(1.0);
        for(i = 0; i < len; i++ ){
                P_R = P_R + base_p[nuc_to_internal(seq[i])];
        }
        *r = P_R;
        return OK;
}

int scan_read_with_pst(struct pst* pst, char* seq, int len, float*r)
{
        int i;
        int  c;
        //float P_PT;
        float P_S;
        float P_P;
        float P_R;

        float P_F;
        //float* base_p = pst->pst_root->nuc_probability;

        float* base_p = pst->fpst_root->prob[0];

        float A;
        float B;

        P_S = prob2scaledprob(1.0);
        P_P = prob2scaledprob(1.0);
        P_R = prob2scaledprob(1.0);
        P_F = prob2scaledprob(1.0);
        //fprintf(stdout,"%s - query sequence\n", seq);
        for(i = 0; i < len; i++ ){
                c = nuc_to_internal(seq[i]);
                //fprintf(stdout,"%f %f %d\n", P,P_R,c);
                P_R = P_R + base_p[c];

                //A = get_pst_prob(pst->pst_root, seq,c, i);
                A = get_fpst_prob(pst->fpst_root, seq, c, i);
                P_S = P_S + A;

                //B = get_ppt_prob(pst->ppt_root, seq,c, i);
                //B = get_fppt_prob(pst->fppt_root, seq, c, i,len);
                //P_P = P_P + B;
                //LOG_MSG("%f %f  fowards (%d / %d)", A,B,i,len);


                //LOG_MSG("%f %f  fowards (%d / %d)", A,B,i,len);
                //P_F = P_F + prob2scaledprob(A*B);//MACRO_MAX(A,B));



//LOG_MSG("%c %f %f %f %f ",seq[i],A,B, P_S,P_P);
        }
        /*A = exp2f(P_S-P_R) / (1.0 + exp2f(P_S-P_R));
        fprintf(stdout,"%f\t%f,%s\n", P_S - P_R,A, seq);
        A = exp2f(P_P-P_R) / (1.0 + exp2f(P_P-P_R));
        fprintf(stdout,"%f\t%f,%s\n", P_P - P_R,A, seq);*/
        //LOG_MSG("%f %f %f %f", P_S,P_P,P_F,P_R);
        A = MACRO_MAX(P_S-P_R, P_P-P_R);
        A = exp2f(A) / (1.0 + exp2f(A));
        *r = A;
        *r = MACRO_MAX(P_S, P_P);
        //exit(0);
        return OK;
}

float get_pst_prob(struct pst_node* n, char* string,int target, int pos)
{
        if(pos == 0){
                return n->probability[target];
        }
        pos = pos -1;
        int c;
        c = nuc_to_internal(string[pos]);
        if(n->next[c]){
                return get_pst_prob(n->next[c], string, target,pos);
        }else{
                return n->probability[target];
        }
}

float get_fpst_prob(struct fpst* f, char* string,int target, int pos)
{
        int n = 0;
        int c;
        int**l = f->links;

        while(pos){
                pos= pos-1;
                c = nuc_to_internal(string[pos]);
                if(!l[n][c]){
                        break;
                }
                n = l[n][c];
        }
        return f->prob[n][target];
}


float get_ppt_prob(struct pst_node* n, char* string,int target, int pos)
{

        if(string[pos+1] == 0){
                return n->probability[target];
        }
        int c;
        pos = pos +1;
        c = nuc_to_internal(string[pos]);
        //c = nuc_code[(int)string[pos]];
        if(n->next[c]){
                return get_ppt_prob(n->next[c], string, target,pos);
        }else{

                return n->probability[target];
        }
}

float get_fppt_prob(struct fpst* f, char* string,int target, int pos, int len)
{
        int n = 0;
        int c;
        int**l = f->links;

        while(pos != len-1){
                pos= pos+1;
                c = nuc_to_internal(string[pos]);
                if(!l[n][c]){
                        break;
                }
                n = l[n][c];
        }
        //fprintf(stdout,"%d %d\n", n, target);
        return f->prob[n][target];
}

int run_build_pst(struct pst** pst, float expected_error, struct count_hash* h, struct kmer_counts* k)//  struct tl_seq_buffer* sb)
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

        //h->len = ;
        RUN(init_pst(&p,expected_error,  h->len,h->L));

        /*  */
        p->r = MACRO_MAX(h->min_prop, p->r);
        
        
        LOG_MSG("%d " , h->L);
        sum = 0.0;
        for(i = 0;i < p->L;i++){
                uint64_t x;
                key =  (1ULL<< 60ULL) | (uint64_t)i;


                hi= kh_get(exact, h->hash, key);
                if(hi != kh_end(h->hash)){
                        LOG_MSG("%d: %d", i,kh_value(h->hash, hi));
                }

                x = p->fpst_root->l;
                p->fpst_root->prob[x][i] = c+1;
                //p->fppt_root->prob[x][i] = c;
                sum += p->fpst_root->prob[x][i];
        }

        for(i = 0;i < p->L;i++){
                x = p->fpst_root->l;
                p->fpst_root->prob[x][i] /=sum;
                p->fpst_root->links[x][i] = 0;
        }
        //fprintf(stdout,"\n");
        //exit(0);
        p->fpst_root->l = 0;
        //p->fppt_root->l = 0;

        RUN(init_helper(&helper,h, p->gamma_min));
        helper = build_pst(p,p->fpst_root,0, helper,k,h);
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
        //print_pst(p,p->ppt_root);
        *pst = p;
        return OK;
ERROR:
        return FAIL;
}

int init_helper(struct pst_node** helper, struct count_hash* ch,float gamma_min)
{
        struct pst_node* h = NULL;
        int i;
        float c;
        float sum;
        khiter_t hi;
        uint64_t x;
        uint8_t* null = NULL;
        RUN(alloc_node(&h,null,0,ch->L));
        sum = 0.0;
        for(i = 0;i < ch->L;i++){
                x =  (1ULL<< 60ULL) | (uint64_t)i;

                c = 0;
                hi= kh_get(exact, ch->hash, x);
                if(hi != kh_end(ch->hash)){
                        //LOG_MSG("%c: %f", "ACGT"[i],kh_value(h->hash, hi));
                        c = kh_value(ch->hash, hi);
                }

                //c = (float) counts[i];
                h->probability[i] = c;
                sum += h->probability[i];
        }
        for(i = 0;i < ch->L;i++){
                h->probability[i] /= sum;
                //h->probability[i] = h->probability[i] *(1.0f  - (float) (ch->L) * gamma_min) + gamma_min;
                //fprintf(stdout,"%d %f\n", i, h->probability[i]);

        }

        *helper = h;
        return OK;
ERROR:
        return FAIL;

}

struct pst_node* build_pst(struct pst* pst,struct fpst*f,int curf, struct pst_node* n,struct kmer_counts* k, struct count_hash* h)
{
        khiter_t hi;
        //char alphabet[] = "ACGT";

        uint8_t tmp[MAX_PST_LEN+20];

        uint8_t ask[MAX_PST_LEN+20];


        int i;
        int j;
        int c;
        int add;

        int len = n->label_len;
        float sum = 0.0f;
        float tmp_counts_s[20];
        float P_ask;
        float P_as;
        float Err;
        uint64_t x;
        uint64_t ask_code;

        uint64_t key;
        int maxlen;

        maxlen = h->len;
        //step 2 test expansion
        //loop though letters at present node

        //LOG_MSG("%d < %d", len+1, maxlen);
        if(len + 1 < maxlen ){
                for(i = 0; i < pst->L;i++){

                        //LOG_MSG("%f , min: %f", n->probability[i], pst->p_min);
                        //if(n->probability[i] >= pst->p_min){ /// string (suffix  + letter is not rare compared to all other letters  ) // e.g. acgC = 0.2 and not acgC = 0.0000001. hence it makes sense to "extend the suffix - test [A,C,G,T] - acgC

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

                        //init longer suffix
                        // here I start constructing the ASK by adding A and then copying the existing suffix S
                        ask_code = (uint64_t)i;
                        ask[0] = i;//alphabet[i];
                        for(j = 1; j < len+1;j++){
                                ask_code = (ask_code << 5ULL) | (uint64_t)n->label[j-1];
                                ask[j] = n->label[j-1];

                        }

                        c = 0;
                        key = ((uint64_t)(len+1) << 60ULL) | ask_code;
                        hi= kh_get(exact, h->hash, key);
                        if(hi != kh_end(h->hash)){
                                //LOG_MSG("%c: %f", "ACGT"[i],kh_value(h->hash, hi));
                                c  = kh_value(h->hash, hi);
                        }
                        P_as = (float) (c+1) / (float)(h->counts_l[len+1] + pst->L);
                        //fprintf(stdout,"counts: %d %f   P_as: %f thres: %f\n",c, h->counts_l[len+1],   P_as,pst->r);

                        if(P_as > pst->r){
                                Err = 0.0f;
                                sum = 0.0f;
                                //to calculate the probabilities I add K
                                for(j = 0; j < pst->L;j++){
                                        ask_code = ask_code << 5ULL | (uint64_t)j;
                                        ask[len+1] = j;
                                        ask[len+2] = 0;

                                        /* now the ask code is complete - I can start looking up their counts to get to P(ask) and P(as) */

                                        //for(c = 0; c < len+2;c++){
                                        //      fprintf(stdout,"%d",ask[c]);
                                        //}

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
                                        //fprintf(stdout,"\t %d P_ask: %f P(K|S):%f \t",c,P_ask, n->probability[j]);
                                        ask_code = ask_code >> 5ULL;

                                        Err += P_ask * logf(P_ask /  ( n->probability[j] * P_as));
                                        //fprintf(stdout,"\t d: %f   error: %f\n",P_ask * logf(P_ask /  ( n->probability[j] * P_as)), Err);

                                }
                                //LOG_MSG("ERR: %f", Err);
                                /* Shall we add this suffix ?  */
                                if (Err > pst->r){
                                        for(j = 0; j < pst->L;j++){
                                                tmp_counts_s[j] = (tmp_counts_s[j] + 1.0f) / (sum + (float)pst->L);
                                        }
                                        f->links[curf][i] = f->l+1;

                                        f->l++;
                                        if(f->l== f->m){
                                                resize_fpst(f,pst->L);
                                        }
                                        RUN(alloc_node(&n->next[i] ,ask,len+1,pst->L));
                                        for(j = 0; j < pst->L;j++){
                                                f->prob[f->l][j] = tmp_counts_s[j];// *(1.0f  - (float)(pst->L) *  pst->gamma_min) + pst->gamma_min;
                                                f->links[f->l][j] = 0;
                                                n->next[i]->probability[j] = tmp_counts_s[j];// *(1.0f  - (float)(pst->L) *  pst->gamma_min) + pst->gamma_min;
                                        }
                                        //LOG_MSG("Going to %d, %d", i,f->links[curf][i]);
                                        n->next[i] = build_pst(pst,f, f->links[curf][i], n->next[i] ,k,h );


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
        p->p_min = 0.1f;
        p->r = 0.00000001f;

        p->fpst_root = NULL;
        p->mean = 0.0f;
        p->var = 0.0f;
        RUN(alloc_fpst(&p->fpst_root, 64,p->L));
        *pst = p;
        return OK;
ERROR:
        return FAIL;
}

void free_pst(struct pst* p)
{
        if(p){
                free_fpst(p->fpst_root);
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
                        fprintf(stdout,"%c","ACGT"[n->label[i]]);
                }
                }else{
                        for(i = 0; i < n->label_len;i++){
                                fprintf(stdout,"%c","ACDEFGHIKLMNPQRSTVWY"[n->label[i]]);
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

int alloc_kmer_counts(struct kmer_counts** k,  int len)
{
        struct kmer_counts* counts = NULL;
        int i,j,c;
        int total;

        MMALLOC(counts, sizeof(struct kmer_counts));

        counts->L = MACRO_MIN(MAX_PST_LEN, len);
        counts->counts = NULL;
        counts->mask = NULL;
        MMALLOC(counts->counts, sizeof(uint32_t*) * counts->L);
        MMALLOC(counts->mask, sizeof(uint32_t) * counts->L);
        c = 1;
        total = 1;
        counts->mask[0] = 0x3;

        for(i = 1;i < counts->L;i++){
                counts->mask[i] = (counts->mask[i-1] << 2) | 0x3;
        }

        for(i = 0;i < counts->L;i++){
                counts->counts[i] = NULL;
                c *= 4;
                MMALLOC(counts->counts[i], sizeof(uint32_t) * c);
                total += (int) (sizeof(uint32_t) * c);
                for(j = 0; j < c;j++){
                        counts->counts[i][j] = 0;
                }
        }
        //LOG_MSG("Total:%d", total);
        *k = counts;
        return OK;
ERROR:
        return FAIL;
}

int add_counts(struct kmer_counts* k , struct tl_seq_buffer* sb)
{
        uint32_t** counts = NULL;
        uint32_t* mask = NULL;
        char* seq;
        uint32_t x;
        int i,j,c,l,len;
        int L;
        counts = k->counts;
        mask = k->mask;
        L = k->L;
        //sb->num_seq = 1;
        for(i = 0; i < sb->num_seq;i++){
                x = 0;
                seq = sb->sequences[i]->seq;
                len = sb->sequences[i]->len;
                //fprintf(stdout,"%s\n", seq);
                l = 0;
                for(j = 0; j < len;j++){
                        x = (x << 2) | nuc_to_internal(seq[j]);
                        for(c = 0; c <= l;c++){
                                //LOG_MSG("Inserting len:%d seq: %d mask %x", c,x,x &mask[c]);
                                counts[c][x & mask[c]]++;
                        }
                        l++;
                        l = MACRO_MIN(l, L-1);
                }
        }
        return OK;
}

int rm_counts(struct kmer_counts* k , struct tl_seq_buffer* sb)
{
        uint32_t** counts = NULL;
        uint32_t* mask = NULL;
        char* seq;
        uint32_t x;
        int i,j,c,l,len;
        int L;
        counts = k->counts;
        mask = k->mask;
        L = k->L;
        //sb->num_seq = 1;

        for(i = 0; i < sb->num_seq;i++){
                x = 0;
                seq = sb->sequences[i]->seq;
                len = sb->sequences[i]->len;
                //fprintf(stdout,"%s\n", seq);
                l = 0;
                for(j = 0; j < len;j++){
                        x = (x << 2) | nuc_to_internal(seq[j]);
                        for(c = 0; c <= l;c++){
                                //LOG_MSG("Inserting len:%d seq: %d mask %x", c,x,x &mask[c]);
                                counts[c][x & mask[c]]--;
                        }
                        l++;
                        l = MACRO_MIN(l, L-1);
                }
        }
        return OK;
}


int test_kmer_counts(struct kmer_counts* k)
{
        int i,j,c;
        int len;

        len = k->L;
        c = 1;
        for(i = 0;i < len;i++){
                c *= 4;
                for(j = 0; j < c;j++){
                        ASSERT(k->counts[i][j] ==0,"%d %d is %d",i,j,k->counts[i][j]);
                }
        }
        return OK;
ERROR:
        return FAIL;
}

void free_kmer_counts(struct kmer_counts* k)
{
        if(k){
                int i;
                for(i = 0;i < k->L;i++){
                        MFREE(k->counts[i]);
                }
                MFREE(k->counts);
                MFREE(k->mask);
                MFREE(k);
        }
}


int alloc_fpst(struct fpst** fast_pst, int size, int L)
{
        struct fpst*f = NULL;
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
        //LOG_MSG("Resizing");
        f->m = f->m + f->m /2;
        RUN(galloc(&f->prob,f->m,L));
        RUN(galloc(&f->links,f->m,L));
        return OK;
ERROR:
        return FAIL;
}

int prob2scaledprob_fpst(struct fpst* f)
{
        int i,j;
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


static inline int nuc_to_internal(const char c)
{
        switch (c) {
        case 'A':
        case 'a':
                return 0;
                break;
        case 'C':
        case 'c':
                return 1;
                break;
        case 'G':
        case 'g':
                return 2;
                break;
        case 'T':
        case 't':
                return 3;
                break;
        case 'N':
        case 'n':
                return 0;
                break;
        default:
                return 0;
                break;
        }
        return -1;
}


/* error  correct stuff  */
/*
struct error_correct_seq{
        char** seq;
        int* error;
        float* P;
        int alloc_num_seq;
        int num_seq;
        int len;
        int pos;
};


static int init_error_correct_seq(struct error_correct_seq** err_seq,int n, int len);
static int resize_error_correct_seq(struct error_correct_seq* e);
static int free_error_correct_seq(struct error_correct_seq** err_seq);

int error_correct(struct pst* pst, struct error_correct_seq* e);


static struct pst_node* get_pst_node(struct pst_node* n,char* s,int pos);

struct pst_node* get_pst_node(struct pst_node* n,char* s,int pos)
{
        int c;
        if(pos == 0){
                return n;
        }
        pos = pos -1;
        c = nuc_to_internal(s[pos]);
        if(n->next[c]){
                return get_pst_node(n->next[c], s,pos);
        }else{
                return n;
        }
}


int init_error_correct_seq(struct error_correct_seq** err_seq,int n, int len)
{
        struct error_correct_seq* e = NULL;
        int i;
        MMALLOC(e,sizeof(struct error_correct_seq));
        e->seq = NULL;
        e->error = NULL;
        e->P = NULL;
        e->alloc_num_seq = n;
        e->num_seq = 0;
        e->len = len;
        e->pos = 0;

        RUN(galloc(&e->seq, e->alloc_num_seq ,e->len+1));
        RUN(galloc(&e->P,e->alloc_num_seq ));
        RUN(galloc(&e->error,e->alloc_num_seq));

        *err_seq = e;
        return OK;
ERROR:
        return FAIL;

}

int resize_error_correct_seq(struct error_correct_seq* e)
{
        e->alloc_num_seq = e->alloc_num_seq + e->alloc_num_seq /2;
        RUN(galloc(&e->seq, e->alloc_num_seq ,e->len+1));
        RUN(galloc(&e->P,e->alloc_num_seq ));
        RUN(galloc(&e->error,e->alloc_num_seq));

        return OK;
ERROR:
        return FAIL;
}

int free_error_correct_seq(struct error_correct_seq** err_seq)
{
        struct error_correct_seq* e = NULL;
        e = *err_seq;
        if(e){
                gfree(e->P);
                gfree(e->seq);
                gfree(e->error);
                MFREE(e);
                e = NULL;
        }
        *err_seq = e;
        return OK;
}

int error_correct(struct pst* pst, struct error_correct_seq* e)
{
        struct pst_node* n;
        float p;
        int i,j,c;

        int cur_numseq;

        cur_numseq = e->num_seq;

        for(i = 0; i < cur_numseq;i++){
                fprintf(stdout,"%d\t",e->pos);
                for(j = 0; j < e->len;j++){
                        fprintf(stdout,"%c", e->seq[i][j]);
                }
                fprintf(stdout,"\n");


                n = get_pst_node(pst->pst_root, e->seq[i], e->pos);

                c = nuc_to_internal(e->seq[i][e->pos]);
                p = n->nuc_probability[c];

                if(e->error[i] < 3){
                        for(j = 0; j < 5;j++){
                                if(j!= c){
                                if(n->nuc_probability[j] * 0.05f > p){
                                        //LOG_MSG("Adding pos%d because %f > %f (%d not %d),", e->pos, n->nuc_probability[j],p,j,c);
                                        strncpy(e->seq[e->num_seq], e->seq[i], e->len);
                                        e->seq[e->num_seq][e->len] = 0;
                                        e->seq[e->num_seq][e->pos] = "ACGTN"[j];
                                        e->error[e->num_seq] = e->error[i] +1;
                                        e->P[e->num_seq] = e->P[i] + prob2scaledprob(n->nuc_probability[j] * 0.05f);
                                        e->num_seq++;
                                        if(e->num_seq == e->alloc_num_seq){
                                                resize_error_correct_seq(e);
                                        }
                                }
                                }
                        }
                }
                p = prob2scaledprob(n->nuc_probability[c]);
                e->P[i] += p;

        }
        e->pos++;
        if(e->pos == e->len){
                return OK;

        }
        return error_correct(pst, e);
}
*/
