
#include "pst.h"
#include "ihmm_seq.h"


static int analyze_label_sequences_with_pst(struct seq_buffer* sb);

int qsort_lcs_cmp(const void *a, const void *b);

int count_int_string(int*p,int** suffix,int h,int len);
int binsearch_down(int*p,int** suffix,int h,int len);
int binsearch_up(int*p,int** suffix,int h,int len);

int cmp_int_array(int* a, int*b, int len);

int lcp_int_array(int* a, int*b);


static struct motif_struct* alloc_motif(int len);

static void* get_state_path(void* ptr);
static long int compare_state_path(void* keyA, void* keyB);
static int resolve_default(void* ptr_a,void* ptr_b);
static void print_motif_struct(void* ptr,FILE* out_ptr);
static void free_motif_struct(void* ptr);


static struct sa* build_sa(struct seq_buffer* sb);
/*
 fp_get = &get_transition;
        fp_cmp = &compare_transition;
        fp_print = &print_fast_t_item_struct;
        fp_cmp_same = &resolve_default;
        fp_free = &free_fast_t_item_struct;

        ft->root = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free);
*/

int analyze_label_sequences_with_pst(struct seq_buffer* sb)
{
        struct sa* sa = NULL;
        struct lcs_count* lcs_count = NULL;

        struct rbtree_root* root = NULL;
        struct motif_struct* motif = NULL;

        int i,c;
        int lower, upper;

        int min_lcp;
        int min_lcp_location;

        void*  (*fp_get)(void* ptr) = NULL;
        long int (*fp_cmp)(void* keyA, void* keyB)= NULL;
        int (*fp_cmp_same)(void* ptr_a,void* ptr_b);
        void (*fp_print)(void* ptr,FILE* out_ptr) = NULL;
        void (*fp_free)(void* ptr) = NULL;
        /* first challenge: build suffix array from labels  */

        ASSERT(sb != NULL, "No sequence Buffer");
        sa = build_sa(sb);


        fp_get = &get_state_path;
        fp_cmp = &compare_state_path;
        fp_print = &print_motif_struct;
        fp_cmp_same = &resolve_default;
        fp_free = &free_motif_struct;

        root = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free);

        MMALLOC(lcs_count, sizeof(struct lcs_count));
        lcs_count->counts = NULL;
        lcs_count->unique = 0;
        lcs_count->len = sb->num_seq;
        lcs_count->cur_longest = 0;

        MMALLOC(lcs_count->counts, sizeof(int) * lcs_count->len);
        for(i= 0 ; i < lcs_count->len;i++){
                lcs_count->counts[i] = 0;
        }


        for(i= 0 ; i < lcs_count->len;i++){
                lcs_count->counts[i] = 0;
        }
        lcs_count->unique = 0;
        lcs_count->cur_longest = 6;
        lcs_count->min_seq = 10; /* Look forstrings in at least 4 sequences. */


        /* sliding window to find lcs */

        lower = 0;
        upper = 0;


        for(i = lower; i <= upper;i++){
                c = sa->lcs[i]->seq_num; /* which sequence is in window */
                if(lcs_count->counts[c] == 0){
                        lcs_count->unique++;
                }
                lcs_count->counts[c]++;
        }


        while(lower <= sa->len){
                /* decide which window to move */
                /* we have enough sequences */
                if(lcs_count->unique >= lcs_count->min_seq){
                        c = sa->lcs[lower]->seq_num;
                        if(lcs_count->counts[c] == 1){
                                lcs_count->unique--;
                        }
                        lcs_count->counts[c]--;
                        lower++;
                }else{/* we don't have enough sequences */
                        upper++;
                        c = sa->lcs[upper]->seq_num;
                        if(lcs_count->counts[c] == 0){
                                lcs_count->unique++;
                        }
                        lcs_count->counts[c]++;
                        if(lcs_count->unique >= lcs_count->min_seq){
                                /* test for longer */
                                min_lcp = 100001;
                                min_lcp_location = -1;
                                for(i = lower+1;i <= upper;i++){
                                        if(min_lcp > sa->lcs[i]->lcp){
                                                min_lcp = sa->lcs[i]->lcp;
                                                min_lcp_location = i;
                                        }
                                }

                                /*if(min_lcp > lcs_count->cur_longest){

                                        //lcs_count->cur_longest = lcs[min_lcp_location]->lcp;
                                }else if(min_lcp == lcs_count->cur_longest){

                                }else{
                                        min_lcp_location = -1;
                                }*/

                                if(min_lcp >=  lcs_count->cur_longest){
                                        /*for(i = lower;i <= upper;i++){


                                                tmp= lcs[i]->str;

                                                fprintf(stdout,"i:%d %d s:%d p:%d\t",i ,lcs[i]->lcp,lcs[i]->seq_num, lcs[i]->pos );
                                                while(*tmp !=  -1){
                                                        fprintf(stdout," %2d", *tmp);
                                                        tmp++;
                                                }
                                                fprintf(stdout,"\n");
                                        }*/
                                        motif = alloc_motif(min_lcp+1);
                                        motif->start_in_sa = lower;
                                        motif->end_in_sa = upper;
                                        motif->len = min_lcp;
                                        for(i= 0; i < min_lcp;i++){
                                                motif->state_list[i] = sa->lcs[min_lcp_location]->str[i];
                                        }
                                        motif->state_list[min_lcp] = -1;

                                        root->tree_insert(root, motif);
                                        /*fprintf(stdout,"RESULT: %d looking between: %d %d (%d)\n",min_lcp_location,lower,upper,lcs_count->unique);
                                        for(i= 0; i < min_lcp;i++){
                                                fprintf(stdout," %3d", lcs[min_lcp_location]->str[i]);
                                        }
                                        fprintf(stdout,"\n");*/

                                }
                        }
                }
                if(upper == sa->len-1){
                        break;
                }

        }


        root->print_tree(root,stdout);
        root->free_tree(root);
        free_sa(sa);
        MFREE(lcs_count->counts);
        MFREE(lcs_count);

        return OK;
ERROR:
        return FAIL;
}

struct sa* build_sa(struct seq_buffer* sb)
{
        struct sa* sa = NULL;
        int total_sequence_len;
        int i,j;

        ASSERT(sb != NULL, "No sequence Buffer");

        MMALLOC(sa, sizeof(struct sa));

        sa->lcs = NULL;
        sa->len = 0;
        total_sequence_len = 0;
        for(i = 0 ;i < sb->num_seq;i++){
                total_sequence_len += sb->sequences[i]->seq_len;
                sb->sequences[i]->label[sb->sequences[i]->seq_len] = -1;

        }


        total_sequence_len += sb->num_seq;
        DPRINTF1("Total length is %d", total_sequence_len);

        sa->alloc_len = total_sequence_len;


        MMALLOC(sa->lcs, sizeof(struct lcs*) * sa->alloc_len);
        for(i = 0; i < sa->alloc_len;i++){
                sa->lcs[i] = NULL;
                MMALLOC(sa->lcs[i], sizeof(struct lcs));
                sa->lcs[i]->lcp = 0;
                sa->lcs[i]->seq_num = -1;
                sa->lcs[i]->str = NULL;
                sa->lcs[i]->pos = -1;
        }

        /* assign pointers  */
        sa->len = 0;
        for(i = 0 ;i < sb->num_seq;i++){
                for(j = 0; j < sb->sequences[i]->seq_len;j++){
                        sa->lcs[sa->len]->str =  &sb->sequences[i]->label[j];
                        sa->lcs[sa->len]->lcp = 0;
                        sa->lcs[sa->len]->seq_num = i;
                        sa->lcs[sa->len]->pos = j;
                        sa->len++;
                }
        }


        qsort( sa->lcs, sa->len, sizeof(struct lcs*), qsort_lcs_cmp);


        for(i = 1; i < sa->len;i++){
                struct lcs* g = NULL;

                g = sa->lcs[i];

                g->lcp = lcp_int_array(sa->lcs[i-1]->str,g->str);
                /*tmp= g->str;
                fprintf(stdout,"i:%d %d s:%d p:%d\t",i ,g->lcp,g->seq_num, g->pos);
                while(*tmp !=  -1){
                        fprintf(stdout," %2d", *tmp);
                        tmp++;
                }
                fprintf(stdout,"\n");*/

        }



        return sa;
ERROR:
        return NULL;

}


void free_sa(struct sa* sa)
{
        if(sa){
                int i;
                for(i = 0; i < sa->alloc_len;i++){
                        MFREE(sa->lcs[i]);
                }
                MFREE(sa->lcs);
                MFREE(sa);
        }
}
static void* get_state_path(void* ptr)
{
        struct  motif_struct* tmp = (struct motif_struct*)  ptr;
        // fprintf(stdout,"got: %p",)
        return tmp->state_list;
}




long int compare_state_path(void* keyA, void* keyB)
{
        int* num1 = (int*)keyA;
        int* num2 = (int*)keyB;
        int i;
        i = 0;
        while(1){
                //fprintf(stdout,"comp:%d %d  \n",num1[i] , num2[i]);
                if(num1[i] != num2[i]){
                        return num1[i] -num2[i];
                }
                if(num1[i] == -1){
                        return 0;
                }
                if(num2[i] == -1){
                        return 0;
                }
                i++;
        }
        return 0;
}

static int resolve_default(void* ptr_a,void* ptr_b)
{

        free_motif_struct(ptr_b);
        return 0;
}


void print_motif_struct(void* ptr,FILE* out_ptr)
{
        struct motif_struct* tmp = (struct motif_struct*)  ptr;
        int i;
        fprintf(out_ptr,"node %d %d %d \n", tmp->len, tmp->start_in_sa,tmp->end_in_sa);
        for(i = 0; i < tmp->len;i++){
                fprintf(out_ptr,"%3d", tmp->state_list[i]);
        }
        fprintf(out_ptr,"\n");
}

struct motif_struct* alloc_motif(int len)
{
        ASSERT(len != 0, "No len");
        struct motif_struct* motif = NULL;

        MMALLOC(motif, sizeof(struct motif_struct ));

        motif->start_in_sa = 0;
        motif->end_in_sa = 0;
        motif->len = len;
        motif->state_list = NULL;
        MMALLOC(motif->state_list, sizeof(int) * len);
        return motif;
ERROR:
        return NULL;
}

void free_motif_struct(void* ptr)
{
        struct motif_struct* tmp = (struct motif_struct*)  ptr;
        if(tmp){
                MFREE(tmp->state_list);
                MFREE(tmp);
        }
}

int count_int_string( int*p, int** suffix,int h,int len)
{
        int a,b;
        //for(i = 0; i < 1000000;i++){
        a = binsearch_down(p,suffix,h,len);
        b = binsearch_up(p,suffix,h,len);
        if(b == h && a != h){
                b++;
        }
        return b-a;
}

int binsearch_down( int*p, int** suffix,int h,int len)
{
        int m = 0;
        int l = 0;
        /*if (t_long_strncmp(p,text+suffix[l],len)<= 0){
          l = l;
          }else */
        if(cmp_int_array(p,suffix[h],len) >  0){
                return h;
        }else{
                while(h-l > 1){
                        //m = (l+h)/2;
                        m = (l + h) >> 1;
                        if(cmp_int_array(p,suffix[m],len) <= 0){
                                h = m;
                        }else{
                                l = m;
                        }
                }
        }
        return l+1;
}

int binsearch_up( int*p, int** suffix,int h,int len)
{
        int m = 0;
        int l = 0;
        /*if (t_long_strncmp(p,text+suffix[l],len)<= 0){
          l = l;
          }else*/
        if(cmp_int_array(p,suffix[h],len) >  0){
                return h;
        }else{
                while(h-l > 1){
                        //m = (l+h)/2;
                        m = (l + h) >> 1;
                        if(cmp_int_array(p,suffix[m],len) < 0){
                                h = m;
                        }else{
                                l = m;
                        }
                }
        }
        return l+1;
}

int lcp_int_array(int* a, int*b)
{
        int i = 0;

        while (a[i] == b[i]){
                if(a[i] == -1){
                        break;
                }
                if(b[i] == -1){
                        break;
                }
                i++;
        }
        return i;
}

int cmp_int_array(int* a, int*b, int len)
{
        int i;

        for(i =0; i < len;i++){
                if(a[i] != b[i]){
                        return a[i]-b[i];
                }

                if(a[i] == -1){
                        return 0;
                }
                if(b[i] == -1){
                        return 0;
                }
        }
        return 0;
}

int qsort_lcs_cmp(const void *a, const void *b)
{
        struct lcs* const *one = a;
        struct lcs* const *two = b;

        int i;
        for (i = 0; (*one)->str[i] == (*two)->str[i];i++){
                if((*one)->str[i] == -1){
                        return 0;
                }
                if((*two)->str[i] == -1){
                        return 0;
                }
        }
        return (*one)->str[i] - (*two)->str[i];
}


int main(const int argc,const char * argv[])
{
        struct seq_buffer* sb = NULL;



        char *tmp_seq[119] = {
                "ACAGGCTAAAGGAGGGGGCAGTCCCCA",
                "AGGCTAAAGGAGGGGGCAGTCCCCACC",
                "AGGCTAAAGGAGGGGGCAGTCCCCACC",
                "AGTCCCCACCATATTTGAGTCTTTCTC",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGT",
                "CTAAAGGAGGGGGCAGTCCCCACCATA",
                "GAGGCTAAAGGAGGGGGCAGTCCCCAT",
                "GAGGGGGCAGTCCCCACCATATTTGAG",
                "GAGGGGGCAGTCCCCACCATATTTGAG",
                "GAGGGGGCAGTCCCCACCATATTTGAG",
                "GAGGGGGCAGTCCCCACCATATTTGAT",
                "GAGGGGGCAGTCCCCACCATATTTGAT",
                "GAGGGGGCAGTCCCCACCATATTTGAT",
                "GAGTCTTTCTCCAAGTTGCGCCGGACA",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTGGATATCACAGGCTAAAGGAGT",
                "GCTAAAGGAGGGGGCAGTCCCCACCAT",
                "GGAGGGGGCAGTCCCCACCATATTTGA",
                "GGAGGGGGCAGTCCCCACCATATTTGA",
                "GGATATCACAGGCTAAAGGAGGGGGCA",
                "GGCAGTCCCCACCATATTTGAGTCTTC",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGGCAGTCACCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGGCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGGCAGTCCCCACCATATTTGAGTCT",
                "GGGGCAGTCCCCACCATATTTGAGTCT",
                "GGGGCAGTCCCCACCATATTTGAGTCT",
                "GGGGGCAGTCCCCACCATATTTGAGTC",
                "GGGGGCAGTCCCCACCATATTTGAGTC",
                "GTCCCCACCATATTTGAGTCTTTCTCT",
                "TCCCCACCATATTTGAGTCTTTCTCCA",
                "TCCCCACCATATTTGAGTCTTTCTCCA",
                "TGAGTCTTTCTCCAAGTTGCGCCGGAT",
                "TGGATATCACAGGCTAAAGGAGGGGGC"};

        //119l
        RUNP(sb = create_ihmm_sequences_mem(tmp_seq ,119));
        RUN(random_label_ihmm_sequences(sb,10,0.3f));
        RUN(analyze_label_sequences_with_pst(sb));
        free_ihmm_sequences(sb);



        return EXIT_SUCCESS;
ERROR:
        free_ihmm_sequences(sb);
        return EXIT_FAILURE;

}
