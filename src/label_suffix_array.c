#include "label_suffix_array.h"

static int qsort_lcs_cmp(const void *a, const void *b);
static int lcp_int_array(int* a, int*b);
static int cmp_int_array(int* a, int*b, int len);

static int binsearch_down(struct lcs** lcs,int h,int*p,int len);
static int binsearch_up(struct lcs** lcs,int h,int*p,int len);

int search_sa(struct sa* sa, int*p, int len, int* start, int*stop)
{
        *start = binsearch_down(sa->lcs, sa->len-1,p,len);
        *stop  = binsearch_up  (sa->lcs,sa->len-1,p,len);
        return OK;
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
                        //if(sb->sequences[i]->label[j] == -1 && i == 0){
                        //        fprintf(stdout,"Seq:%d pos: %d is -1\n", i,j);
                        //}
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
                //int* tmp = NULL;
                g = sa->lcs[i];

                g->lcp = lcp_int_array(sa->lcs[i-1]->str,g->str);
                /*tmp= g->str;
                fprintf(stdout,"i:%d %d s:%d p:%d\t",i ,g->lcp,g->seq_num, g->pos);
                for(j = 0; j < 20;j++){
                        fprintf(stdout," %2d", tmp[j]);
                        if(tmp[j] == -1){
                                break;
                        }
                }

                fprintf(stdout,"\n");*/

        }
        return sa;
ERROR:
        return NULL;

}

int binsearch_down(struct lcs** lcs,int h, int*p, int len)
{
        int m = 0;
        int l = 0;
        /*if (t_long_strncmp(p,text+suffix[l],len)<= 0){
          l = l;
          }else */
        if(cmp_int_array(p,lcs[h]->str  ,len) >  0){
                return h;
        }else{
                while(h-l > 1){
                        //m = (l+h)/2;
                        m = (l + h) >> 1;
                        if(cmp_int_array(p,lcs[m]->str,len) <= 0){
                                h = m;
                        }else{
                                l = m;
                        }
                }
        }
        return l+1;
}

int binsearch_up(struct lcs** lcs,int h,int*p,int len)
{
        int m = 0;
        int l = 0;
        /*if (t_long_strncmp(p,text+suffix[l],len)<= 0){
          l = l;
          }else*/
        if(cmp_int_array(p,lcs[h]->str,len) >  0){
                return h;
        }else{
                while(h-l > 1){
                        //m = (l+h)/2;
                        m = (l + h) >> 1;
                        if(cmp_int_array(p,lcs[m]->str,len) < 0){
                                h = m;
                        }else{
                                l = m;
                        }
                }
        }
        return l+1;
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
