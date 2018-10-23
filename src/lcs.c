#include "lcs.h"
#include "ihmm_seq.h"
#include "finite_hmm.h"
#include <getopt.h>

struct parameters{
        char* input;
        char* out;
        float min_frac_occur;
        int min_pattern_len;
        int pos_resolution;
};

static int run_lcs(struct parameters* param);

//static int analyze_label_sequences_with_pst(struct seq_buffer* sb);
static struct rbtree_root* find_lcs(struct sa* sa, int total_sequences,int min_len, float min_seq, int all);

static int calculate_relative_entropy(struct motif_struct** list, int num_motif, struct fhmm* fhmm);
static int add_positional_distribution(struct motif_struct** list, int num_motif, struct sa* sa, struct seq_buffer*sb, int pos_resolution);

static int write_lcs_motif_data_for_plotting(char* filename,struct sa* sa, struct seq_buffer*sb,struct motif_struct** res_motif, int num_res_motif, struct fhmm*  fhmm);

static int find_overlapping_motifs(struct motif_struct** res_motif, int num_res_motif, struct fhmm*  fhmm);
static int compare_motif(struct motif_struct* a,struct motif_struct* b,struct fhmm*  fhmm, float* KLdivergence);

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


static int free_parameters(struct parameters* param);
static int print_help(char **argv);

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        tlog.echo_build_config();

        MMALLOC(param, sizeof(struct parameters));
        param->min_pattern_len = 10;
        param->pos_resolution = 100;
        param->out = NULL;
        param->input = NULL;
        param->min_frac_occur = 0.5;

        while (1){
                static struct option long_options[] ={
                        {"in",required_argument,0,'i'},
                        {"d",required_argument,0,'d'},
                        {"frac",required_argument,0,'f'},
                        {"len",required_argument,0,'l'},
                        {"out",required_argument,0,'o'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"i:o:l:d:hf:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 'i':
                        param->input = optarg;
                        break;
                case 'f':
                        param->min_frac_occur = atof(optarg);
                        break;

                case 'l':
                        param->min_pattern_len = atoi(optarg);
                        break;
                case 'd':
                        param->pos_resolution = atoi(optarg);
                        break;
                case 'o':
                        param->out = optarg;
                        break;

                case 'h':
                        RUN(print_help(argv));
                        MFREE(param);
                        exit(EXIT_SUCCESS);
                        break;
                default:
                        ERROR_MSG("not recognized");
                        break;
                }
        }

        LOG_MSG("Starting run");

        if(!param->pos_resolution){
                RUN(print_help(argv));
                ERROR_MSG("Numseq is 0! use --n 1 (or more!).");

        }
        if(!param->input){
                RUN(print_help(argv));
                ERROR_MSG("No input file! use -i <blah.h5>");

        }
        if(!param->out){
                RUN(print_help(argv));
                ERROR_MSG("No output file! use -o <blah.h5>");
        }

        ASSERT(param->min_pattern_len > 3, "Minimum pattern length has to be > 3");
        RUN(run_lcs(param));

        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}

int free_parameters(struct parameters* param)
{
        ASSERT(param != NULL, " No param found - free'd already???");

        MFREE(param);
        return OK;
ERROR:
        return FAIL;

}

int print_help(char **argv)
{
        const char usage[] = " -in <h5 model> -out <h5 out>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--len","Minimum pattern length." ,"[10]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--frac","Minimum fraction of seq containing pattern." ,"[0.5]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-d","Resolution." ,"[100]"  );
        return OK;
}

int run_lcs(struct parameters* param)
{
        struct seq_buffer* sb = NULL;
        struct sa* sa = NULL;

        struct rbtree_root* root = NULL;
        struct fhmm* fhmm = NULL;

        int i,j,c;
        int iter;
        int seq,pos;

        struct motif_struct** res_motif = NULL;
        struct motif_struct* motif_ptr = NULL;
        int num_alloc_res_motif = 1024;
        int num_res_motif = 0;

        MMALLOC(res_motif, sizeof(struct motif_struct*) * num_alloc_res_motif);
        for(i = 0; i < num_alloc_res_motif;i++){
                res_motif[i] = NULL;
        }


        RUNP(sb = get_sequences_from_hdf5_model(param->input));

        ASSERT(sb != NULL, "No sequence Buffer");

        /* try 512 times to get lcs -  */
        for(iter = 0;iter < 512;iter++){

                /* build / re-build suffix array */
                LOG_MSG("Building SA");
                sa = build_sa(sb);
                LOG_MSG("Done");

                LOG_MSG("LCA");
                RUNP(root = find_lcs(sa, sb->num_seq,  param->min_pattern_len, param->min_frac_occur,0));
                LOG_MSG("Done");

                if(root->num_entries){ /* copy motif instances  */
                        /* sort - note this will be done on length only. Rel entropy will be done later */
                        qsort(root->data_nodes, root->num_entries, sizeof(struct motif_struct*), qsort_motif_struct);
                        motif_ptr = root->data_nodes[0];
                        LOG_MSG("Found motif of length: %d",motif_ptr->len );
                        for(i = 0; i < root->num_entries;i++){
                                motif_ptr = root->data_nodes[i];
                                res_motif[num_res_motif] = alloc_motif(motif_ptr->len+1);
                                res_motif[num_res_motif]->len = motif_ptr->len;
                                res_motif[num_res_motif]->start_in_sa = motif_ptr->start_in_sa;
                                res_motif[num_res_motif]->end_in_sa = motif_ptr->end_in_sa;
                                res_motif[num_res_motif]->n_occur =  motif_ptr->n_occur;
                                for(j = 0; j < motif_ptr->len;j++){
                                        res_motif[num_res_motif]->state_list[j] = motif_ptr->state_list[j];
                                }
                                res_motif[num_res_motif]->state_list[motif_ptr->len] = -1;

                                /* Nuke motif in sb..  */
                                fprintf(stdout,"%d %d\n",motif_ptr->start_in_sa,motif_ptr->end_in_sa);
                                for(j = motif_ptr->start_in_sa; j <= motif_ptr->end_in_sa;j++){
                                        //fprintf(stdout,"Nuking: seq%d pos%d\n",sa->lcs[j]->seq_num,sa->lcs[j]->pos);
                                        seq = sa->lcs[j]->seq_num;
                                        pos = sa->lcs[j]->pos;
                                        for(c = 0; c < motif_ptr->len;c++){
                                                sb->sequences[seq]->label[pos+c] = -1;
                                        }
                                }
                                num_res_motif++;
                                if(num_res_motif == num_alloc_res_motif){
                                        j = num_alloc_res_motif;
                                        num_alloc_res_motif = num_alloc_res_motif << 1;
                                        MREALLOC(res_motif, sizeof(struct motif_struct*) * num_alloc_res_motif);
                                        for(c = j; c < num_alloc_res_motif;c++){
                                                res_motif[c] = NULL;
                                        }

                                }

                        }

                }else{          /* clean up and exit loop  */
                        LOG_MSG("Found nothing");
                        free_sa(sa);
                        sa = NULL;
                        root->free_tree(root);
                        root = NULL;
                        break;
                }
                root->free_tree(root);
                root = NULL;
                free_sa(sa);
                sa = NULL;

        }
        /* remove old sequence buffer and load a fresh one from file  */
        free_ihmm_sequences(sb);
        sb = NULL;
        RUNP(sb = get_sequences_from_hdf5_model(param->input));

        /* check if we found motifs!!!! */
        if(num_res_motif){
                sa = build_sa(sb);


                /* I can do this later...  */
                RUNP(fhmm = alloc_fhmm());
                /* get HMM parameters  */
                RUN(read_hmm_parameters(fhmm,param->input));

                RUN(calculate_relative_entropy(res_motif, num_res_motif,fhmm));

                qsort(res_motif, num_res_motif, sizeof(struct motif_struct*), qsort_motif_struct);

                RUN(find_overlapping_motifs(res_motif, num_res_motif, fhmm));

                RUN(add_positional_distribution(res_motif, num_res_motif , sa, sb, param->pos_resolution));

                RUN(write_lcs_motif_data_for_plotting(param->out,sa,sb,res_motif, num_res_motif, fhmm));

                for(i = 0; i < num_res_motif;i++){
                        free_motif_struct( (void*) res_motif[i]);
                }
                MFREE(res_motif);

                free_fhmm(fhmm);
                //root->free_tree(root);
                free_sa(sa);
        }
        free_ihmm_sequences(sb);

        return OK;
ERROR:
        return FAIL;

}


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

        RUNP(root = find_lcs(sa, sb->num_seq,  min_pattern_len, min_seq_occur,0));

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

        RUN(calculate_relative_entropy((struct motif_struct**)root->data_nodes, root->num_entries,fhmm));

        qsort(root->data_nodes, root->num_entries, sizeof(struct motif_struct*), qsort_motif_struct);

        /* look for overlapping motifs */
        RUN(add_positional_distribution((struct motif_struct**)root->data_nodes, root->num_entries, sa, sb, pos_resolution));

        RUN(write_lcs_motif_data_for_plotting("mysmalltest.h5",sa,sb,(struct motif_struct**)root->data_nodes, root->num_entries, fhmm));
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

int find_overlapping_motifs(struct motif_struct** res_motif, int num_res_motif, struct fhmm*  fhmm)
{
        int i,j;
        float x;
        ASSERT(res_motif != NULL, "No motifs");

        for(i = 0; i < num_res_motif;i++){
                for(j = i +1; j < num_res_motif;j++){
                        RUN(compare_motif(res_motif[i], res_motif[j],fhmm, &x));
                        if(x <= 0.5){
                                if(res_motif[i]->len > res_motif[j]->len){
                                        res_motif[j]->mean_rel_entropy = 0.0;
                                }else{
                                        res_motif[i]->mean_rel_entropy = 0.0;
                                }
                        }
                        //print_motif_struct(res_motif[i], stdout);
                        //print_motif_struct(res_motif[j], stdout);
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int compare_motif(struct motif_struct* a,struct motif_struct* b,struct fhmm*  fhmm, float* KLdivergence)
{
        struct motif_struct* tmp = NULL;
        int i,j,c;
        float x;
        float col;
        int sa;
        int sb;
        ASSERT(a != NULL, "No A motif");
        ASSERT(b != NULL, "No B motif");
        *KLdivergence = FLT_MAX;

        if(a->mean_rel_entropy <= 0.05 || b->mean_rel_entropy <= 0.05){
                return OK;
        }
        /* swap if b longer than a */
        if(a->len < b->len){
                tmp = a;
                a = b;
                b = tmp;
        }
        //print_motif_struct(a, stdout);
        //print_motif_struct(b, stdout);
        for(i = 0; i <= a->len - b->len;i++){

                x = 0.0f;
                for(j = 0; j < b->len;j++){
                        sb = b->state_list[j];
                        sa = a->state_list[i+j];

                        col = 0.0;
                        for(c = 0; c < fhmm->L;c++){
                                if(fhmm->e[sb][c] > 0.0f){
                                        col  += fhmm->e[sb][c] * log2f(fhmm->e[sb][c] / fhmm->e[sa][c]);
                                }
                        }
                        //LOG_MSG("%d vs %d : %f",  sb,sa,col);
                        x += col;
                }
                x = x / (float) b->len;
                //LOG_MSG("SCORE:%f",x);
                if(x < *KLdivergence){
                        *KLdivergence = x;
                }
        }
        //LOG_MSG("SCORE:%f",*KLdivergence);

        return OK;
ERROR:
        return FAIL;
}

int write_lcs_motif_data_for_plotting(char* filename,struct sa* sa, struct seq_buffer*sb,struct motif_struct** res_motif, int num_res_motif, struct fhmm*  fhmm)
{
        char buffer[BUFFER_LEN];
        struct hdf5_data* hdf5_data = NULL;
        struct motif_struct* m = NULL;
        struct lcs* g = NULL;

        float** tmp = NULL;
        int i,j,c;
        uint8_t* seq;



        ASSERT(filename != NULL, "No file name");
        ASSERT(res_motif != NULL, "No motifs");
        ASSERT(num_res_motif != 0, "No Motifs");


        /* determine max_len of motif - to alloc frequency matrix */

        RUNP(hdf5_data = hdf5_create());
        snprintf(buffer, BUFFER_LEN, "%s",PACKAGE_NAME );
        hdf5_add_attribute(hdf5_data, "Program", buffer, 0, 0.0f, HDF5GLUE_CHAR);
        snprintf(buffer, BUFFER_LEN, "%s", PACKAGE_VERSION);
        hdf5_add_attribute(hdf5_data, "Version", buffer, 0, 0.0f, HDF5GLUE_CHAR);


        RUN(hdf5_create_file(filename,hdf5_data));

        hdf5_write_attributes(hdf5_data, hdf5_data->file);


        RUN(hdf5_create_group("MotifData",hdf5_data));

        for(c = 0; c < num_res_motif;c++){
                m = res_motif[c];
                if(m->mean_rel_entropy > 0.5){


                        RUNP(tmp = malloc_2d_float(tmp,m->len, fhmm->L, 0.0));
                        m->start_in_sa = binsearch_down(m->state_list, sa->lcs, sa->len-1, m->len);
                        m->end_in_sa =   binsearch_up  (m->state_list,sa->lcs,sa->len-1, m->len);

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


int add_positional_distribution(struct motif_struct** list, int num_motif, struct sa* sa, struct seq_buffer*sb, int pos_resolution)
{
        struct motif_struct* m = NULL;
        struct lcs* g = NULL;
        //int* tmp = NULL;
        int i,j,c;
        int index;
        int l;
        for(i = 0; i < num_motif;i++){
                m = list[i];
                print_motif_struct(m, stdout);
                MMALLOC(m->occ, sizeof(int) * (pos_resolution+1));
                for(j = 0; j <  pos_resolution+1;j++){
                        m->occ[j] = 0;
                }
                for(j = m->start_in_sa; j < m->end_in_sa;j++){
                        g = sa->lcs[j];
                        l = sb->sequences[g->seq_num]->seq_len;
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

int calculate_relative_entropy(struct motif_struct** list, int num_motif, struct fhmm* fhmm)
{
        int i,j,c,s;
        struct motif_struct* m = NULL;
        float* background = NULL;
        float x;

        ASSERT(list!= NULL, "No rbtree");
        ASSERT(fhmm != NULL, "No fhmm");

        background = fhmm->background;
        for(i = 0; i < num_motif;i++){
                m = list[i];
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


struct rbtree_root* find_lcs(struct sa* sa, int total_sequences,int min_len, float min_seq, int all)
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
        if(min_seq <= 1.0){
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

                                if(min_lcp > lcs_count->cur_longest && !all){

                                        lcs_count->cur_longest = sa->lcs[min_lcp_location]->lcp;
                                        //root->free_tree(root);

                                }

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
                        if(sb->sequences[i]->label[j] == -1 && i == 0){
                                fprintf(stdout,"Seq:%d pos: %d is -1\n", i,j);
                        }
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
