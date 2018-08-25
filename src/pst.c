
#include "pst.h"
#include "ihmm_seq.h"
#include "finite_hmm.h"

//static int analyze_label_sequences_with_pst(struct seq_buffer* sb);
static struct rbtree_root* find_lcs(struct sa* sa, int total_sequences,int min_len, float min_seq);

static int calculate_relative_entropy(struct rbtree_root* root, struct fhmm* fhmm);
static int add_positional_distribution(struct rbtree_root* root, struct sa* sa,struct seq_buffer* sb, int pos_resolution);

static int write_lcs_motif_data_for_plotting(char* filename, struct sa* sa,struct seq_buffer*sb,struct rbtree_root* root, struct fhmm*  fhmm);

static int qsort_lcs_cmp(const void *a, const void *b);
static int qsort_motif_struct(const void *a, const void *b);
int count_int_string(int*p, struct lcs** lcs,int h,int len);
int binsearch_down(int*p,struct lcs** lcs,int h,int len);
int binsearch_up(int*p,struct lcs** lcs,int h,int len);

int cmp_int_array(int* a, int*b, int len);

int lcp_int_array(int* a, int*b);


static struct motif_struct* alloc_motif(int len);

static void* get_state_path(void* ptr);
static long int compare_state_path(void* keyA, void* keyB);
static int resolve_default(void* ptr_a,void* ptr_b);
static void print_motif_struct(void* ptr,FILE* out_ptr);
static void free_motif_struct(void* ptr);

static struct sa* build_sa(struct seq_buffer* sb);
static void free_sa(struct sa* sa);

int analyze_label_sequences_with_pst(char* filename, int min_pattern_len, float min_seq_occur, int pos_resolution)
{
        struct seq_buffer* sb = NULL;
        struct sa* sa = NULL;

        struct rbtree_root* root = NULL;
        struct fhmm* fhmm = NULL;


        /* first challenge: build suffix array from labels  */
        RUNP(sb = get_sequences_from_hdf5_model(filename));

        ASSERT(sb != NULL, "No sequence Buffer");
        sa = build_sa(sb);

        RUNP(root = find_lcs(sa, sb->num_seq,  min_pattern_len, min_seq_occur));

        /* check if we found motifs!!!! */

        if(root->num_entries == 0){
                WARNING_MSG("No motif found.");
                free_fhmm(fhmm);
                root->free_tree(root);
                free_sa(sa);
                free_ihmm_sequences(sb);
                return OK;
        }

        RUNP(fhmm = alloc_fhmm());
        /* get HMM parameters  */
        RUN(read_hmm_parameters(fhmm,filename));

        RUN(calculate_relative_entropy(root,fhmm));

        qsort(root->data_nodes, root->num_entries, sizeof(struct motif_struct*), qsort_motif_struct);

        RUN(add_positional_distribution( root, sa, sb, pos_resolution));

        RUN(write_lcs_motif_data_for_plotting("mysmalltest.h5",sa,sb, root, fhmm));
        /*int i;
         for(i = 0; i < root->num_entries;i++){
                print_motif_struct(root->data_nodes[i], stdout);
                }*/
        free_fhmm(fhmm);
//        root->print_tree(root,stdout);
        root->free_tree(root);
        free_sa(sa);
        free_ihmm_sequences(sb);

        return OK;
ERROR:
        return FAIL;

}

int write_lcs_motif_data_for_plotting(char* filename,struct sa* sa, struct seq_buffer*sb, struct rbtree_root* root, struct fhmm*  fhmm)
{
        char buffer[BUFFER_LEN];
        struct hdf5_data* hdf5_data = NULL;
        struct motif_struct* m = NULL;
        struct lcs* g = NULL;

        float** tmp = NULL;
        int i,j,c;
        uint8_t* seq;



        ASSERT(filename != NULL, "No file name");
        ASSERT(root != NULL, "No motifs");
        ASSERT(root->num_entries != 0, "No Motifs");


        /* determine max_len of motif - to alloc frequency matrix */

        RUNP(hdf5_data = hdf5_create());
        snprintf(buffer, BUFFER_LEN, "%s",PACKAGE_NAME );
        hdf5_add_attribute(hdf5_data, "Program", buffer, 0, 0.0f, HDF5GLUE_CHAR);
        snprintf(buffer, BUFFER_LEN, "%s", PACKAGE_VERSION);
        hdf5_add_attribute(hdf5_data, "Version", buffer, 0, 0.0f, HDF5GLUE_CHAR);


        RUN(hdf5_create_file(filename,hdf5_data));

        hdf5_write_attributes(hdf5_data, hdf5_data->file);


        RUN(hdf5_create_group("MotifData",hdf5_data));

        for(c = 0; c < root->num_entries;c++){
                m = root->data_nodes[c];


                RUNP(tmp = malloc_2d_float(tmp,m->len, fhmm->L, 0.0));


                for(i = m->start_in_sa; i < m->end_in_sa;i++){
                        g = sa->lcs[i];

                        seq = sb->sequences[g->seq_num]->seq + g->pos;
                        //g->possb->

                        for(j = 0; j < m->len;j++){
                                tmp[j][seq[j]]++;
                                //       c = sb->sequences[g->seq_num]->dd;

                        }
                        //l = sb->sequences[ g->seq_num]->seq_len;
                        //index =roundf((float)pos_resolution * ((float) g->pos / l));
                        //m->occ[index]++;
                        /*tmp= g->str;
                        fprintf(stdout,"i:%d %d s:%d p:%d\t",i ,g->lcp,g->seq_num, g->pos);
                        for(c = 0; c < m->len+2;c++){
                                fprintf(stdout," %2d", tmp[c]);
                                if(tmp[c] == -1){
                                        break;
                                }
                        }
                        fprintf(stdout,"\n");*/

                }

                        /*for(i = 0; i < fhmm->L;i++){
                        for(j = 0; j < m->len;j++){
                                tmp[j][i] = fhmm->e[m->state_list[j]][i];
                        }
                        }*/

                hdf5_data->rank = 2;
                hdf5_data->dim[0] = m->len;
                hdf5_data->dim[1] = fhmm->L;
                hdf5_data->chunk_dim[0] = m->len;
                hdf5_data->chunk_dim[1] = fhmm->L;
                hdf5_data->native_type = H5T_NATIVE_FLOAT;
                snprintf(buffer, BUFFER_LEN, "Motif%06d", c+1);
                hdf5_write(buffer,&tmp[0][0], hdf5_data);

                free_2d((void**) tmp);
                tmp = NULL;

        }

        RUN(hdf5_close_group(hdf5_data));
        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);
        return OK;
ERROR:
        if(hdf5_data){
                hdf5_close_file(hdf5_data);
                hdf5_free(hdf5_data);
        }
        return FAIL;
}


int add_positional_distribution(struct rbtree_root* root, struct sa* sa, struct seq_buffer*sb, int pos_resolution)
{
        struct motif_struct* m = NULL;
        struct lcs* g = NULL;
        //int* tmp = NULL;
        int i,j,c;
        int index;
        int l;
        for(i = 0; i < root->num_entries;i++){
                m = root->data_nodes[i];
                print_motif_struct(m, stdout);
                MMALLOC(m->occ, sizeof(int) * (pos_resolution+1));
                for(j = 0; j <  pos_resolution+1;j++){
                        m->occ[j] = 0;
                }
                for(j = m->start_in_sa; j < m->end_in_sa;j++){
                        g = sa->lcs[j];
                        l = sb->sequences[ g->seq_num]->seq_len;
                        index =roundf((float)pos_resolution * ((float) g->pos / l));
                        m->occ[index]++;
                        /*tmp= g->str;
                        fprintf(stdout,"i:%d %d s:%d p:%d\t",i ,g->lcp,g->seq_num, g->pos);
                        for(c = 0; c < m->len+2;c++){
                                fprintf(stdout," %2d", tmp[c]);
                                if(tmp[c] == -1){
                                        break;
                                }
                        }
                        fprintf(stdout,"\n");*/

                }
                c = 0;
                for(j = 0; j <  pos_resolution+1;j++){
                        c+= m->occ[j];
                        fprintf(stdout,"%4d",m->occ[j]);
                }
                fprintf(stdout," TOTAL:%d\n",c);

                /*c = (float) sb->sequences[i]->label[j];
                  index = roundf(1000.0f * ((float) j / l));
                  matrix[c][index] += 1.0f;
                  state_sums[c] += 1.0f;*/

        }
        return OK;
ERROR:
        return FAIL;
}

int calculate_relative_entropy(struct rbtree_root* root, struct fhmm* fhmm)
{
        int i,j,c,s;
        struct motif_struct* m = NULL;
        float* background = NULL;
        float x;

        ASSERT(root!= NULL, "No rbtree");
        ASSERT(fhmm != NULL, "No fhmm");

        background = fhmm->background;
        for(i = 0; i < root->num_entries;i++){
                m = root->data_nodes[i];
                x = 0.0f;
                for(j = 0; j < m->len;j++){
                        s = m->state_list[j];

                        for(c = 0; c < fhmm->L;c++){
                                if(fhmm->e[s][c] > 0.0f){
                                        x += fhmm->e[s][c] * log2f(fhmm->e[s][c] / background[c]);
                                }
                        }
                }
                x = x / (float) m->len;
                m->mean_rel_entropy = x;
        }
        return OK;
ERROR:
        return FAIL;
}


struct rbtree_root* find_lcs(struct sa* sa, int total_sequences,int min_len, float min_seq)
{
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

        ASSERT(sa != NULL, "No suffix array!");
        ASSERT(min_seq < total_sequences, "too many minseq");

        fp_get = &get_state_path;
        fp_cmp = &compare_state_path;
        fp_print = &print_motif_struct;
        fp_cmp_same = &resolve_default;
        fp_free = &free_motif_struct;

        root = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free);

        MMALLOC(lcs_count, sizeof(struct lcs_count));
        lcs_count->counts = NULL;
        lcs_count->unique = 0;
        lcs_count->len = total_sequences;
        lcs_count->cur_longest = 0;

        MMALLOC(lcs_count->counts, sizeof(int) * lcs_count->len);
        for(i= 0 ; i < lcs_count->len;i++){
                lcs_count->counts[i] = 0;
        }


        for(i= 0 ; i < lcs_count->len;i++){
                lcs_count->counts[i] = 0;
        }
        if(min_seq < 1.0){
                min_seq = (int) (min_seq * (float) total_sequences);
        }
        lcs_count->unique = 0;
        lcs_count->cur_longest = min_len;
        lcs_count->min_seq = min_seq; /* Look forstrings in at least 4 sequences. */


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
                                /*
                                if(min_lcp > lcs_count->cur_longest){

                                        lcs_count->cur_longest = sa->lcs[min_lcp_location]->lcp;
                                        //root->free_tree(root);

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

        MFREE(lcs_count->counts);
        MFREE(lcs_count);
        if(root->num_entries){

                /* flatten tree  */
                root->flatten_tree(root);

                /* count number of occurances...  */

                for(i = 0; i < root->num_entries;i++){
                        motif = root->data_nodes[i];

                        motif->start_in_sa = binsearch_down(motif->state_list, sa->lcs, sa->len-1, motif->len);
                        motif->end_in_sa =   binsearch_up  (motif->state_list,sa->lcs,sa->len-1, motif->len);

                        motif->n_occur =  motif->end_in_sa - motif->start_in_sa;
                        //print_motif_struct(root->data_nodes[i], stdout);
                }

        }
        return root;
ERROR:
        return NULL;

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
        fprintf(out_ptr,"node %d %d %d KL:%f (%d)\n", tmp->len, tmp->start_in_sa,tmp->end_in_sa, tmp->mean_rel_entropy,  tmp->n_occur);
        for(i = 0; i < tmp->len;i++){
                fprintf(out_ptr,"%4d", tmp->state_list[i]);
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
        motif->mean_rel_entropy = 0.0;
        motif->n_occur = 0;
        motif->occ = NULL;

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
                if(tmp->occ){
                        MFREE(tmp->occ);
                }
                MFREE(tmp->state_list);
                MFREE(tmp);
        }
}

int count_int_string( int*p, struct lcs** lcs,int h,int len)
{
        int a,b;

        //for(i = 0; i < 1000000;i++){
        a = binsearch_down(p,lcs,h,len);
        b = binsearch_up(p,lcs,h,len);
        if(b == h && a != h){
                b++;
        }
        return b-a;
}

int binsearch_down( int*p, struct lcs** lcs,int h,int len)
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

int binsearch_up( int*p, struct lcs** lcs,int h,int len)
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


int qsort_motif_struct(const void *a, const void *b)
{
        struct motif_struct* const *one = a;
        struct motif_struct* const *two = b;

        if((*one)->len == (*two)->len){
                if((*one)->mean_rel_entropy > (*two)->mean_rel_entropy){
                        return -1;
                }else{
                        return 1;
                }

        }
        return  (*two)->len - (*one)->len;
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
        RUN(analyze_label_sequences_with_pst("../upstream1000_selection.fa.h5", 12, 0.2, 100));
        free_ihmm_sequences(sb);



        return EXIT_SUCCESS;
ERROR:
        free_ihmm_sequences(sb);
        return EXIT_FAILURE;

}
