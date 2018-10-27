#include "infoclust.h"
#include "ihmm_seq.h"
#include "finite_hmm.h"
#include "label_suffix_array.h"
#include "run_score.h"
#include <getopt.h>

struct parameters{
        char* input;
        char* out;
        float min_d_increase;
        int maxlen;
};

static int run_infoclust(struct parameters* param);

int calc_per_state_rel_entrophy(struct fhmm* fhmm, float* rel_entropy);

int calculate_relative_entropy(struct motif_list* m, struct fhmm* fhmm);

int calculate_log_likelihood(struct fhmm* fhmm,struct seq_buffer* sb, struct paraclu_cluster* motif,struct sa* sa);

int add_present_in_seq(struct paraclu_cluster* motif,struct sa* sa, int numseq);
int find_overlapping_hits(struct sa* sa, int start, int stop, int w);

int insert_motif(struct motif_list*m, struct paraclu_cluster* p,struct fhmm* fhmm);
int compare_motif(struct fhmm* fhmm, struct paraclu_cluster* a,struct paraclu_cluster* b, float* kl_div);

int merge_motif(struct paraclu_cluster* a,struct paraclu_cluster* b);
int write_para_motif_data_for_plotting(char* filename,struct motif_list* m, struct fhmm*  fhmm);


struct motif_list* init_motif_list(int initial_len);
int extend_motif_list(struct motif_list* m);
void free_motif_list(struct motif_list* m);

struct paraclu_cluster* init_paraclu_cluster(void);
void free_paraclu_cluster(struct paraclu_cluster* p);

int sort_paraclu_cluster_based_on_likelihood(const void *a, const void *b);
int max_score_segment(float* x , int start ,int end, int min_len,float min_density,struct motif_list* m);
float weakestPrefix(float* x , int start ,int end, int* min_prefix_pos, float* min_prefix);
void weakestSuffix(float* x , int start ,int end, int* min_suffix_pos, float* min_suffix);
int free_parameters(struct parameters* param);
int print_help(char **argv);

/* Here I plan to:
1) project information content of states onto labelled sequences
2) visualise
3) use a paraclu like algorithm to look for clusters of high information content
4) islands of high IC should somehow be compared across sequences to find shared patterns
 */

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        tlog.echo_build_config();

        MMALLOC(param, sizeof(struct parameters));
        param->out = NULL;
        param->input = NULL;
        param->maxlen = 20;
        param->min_d_increase = 1.00;


        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"out",required_argument,0,'o'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"m:o:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 'm':
                        param->input = optarg;
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

        if(!param->input){
                RUN(print_help(argv));
                ERROR_MSG("No input file! use -i <blah.h5>");

        }
        if(!param->out){
                RUN(print_help(argv));
                ERROR_MSG("No output file! use -o <blah.h5>");
        }

        RUN(run_infoclust(param));
        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}


int run_infoclust(struct parameters* param)
{
        struct seq_buffer* sb = NULL;
        struct ihmm_sequence* s = NULL;
        struct fhmm* fhmm = NULL;
        struct fhmm* fhmm_log = NULL;
        struct paraclu_cluster* p =  NULL;
        struct motif_list* m = NULL;
        struct motif_list* t = NULL;
        struct sa* sa = NULL;
        float* rel_entropy = NULL;
        float sum = 0;
        int i,j,c;
        int old_end;
        ASSERT(param != NULL, "No parameters.");

        /* read in sequences */
        RUNP(sb = get_sequences_from_hdf5_model(param->input));
        ASSERT(sb != NULL, "No sequence Buffer");


        RUNP(fhmm_log = alloc_fhmm());
        /* get HMM parameters  */
        RUN(read_hmm_parameters(fhmm_log,param->input));


        RUNP(fhmm = alloc_fhmm());
        /* get HMM parameters  */
        RUN(read_hmm_parameters(fhmm,param->input));

        RUN(setup_model(fhmm_log));

        RUN(run_label_sequences(fhmm_log,sb, 8));

        /* re-estimate emissions */


        for(i = 0; i < fhmm_log->K;i++){
                for(j = 0 ; j < fhmm_log->L;j++){
                        fhmm_log->e[i][j] = 0.0;
                        fhmm->e[i][j] = 0.0;
                }
        }
        for(i = 0;i < sb->num_seq;i++){
                s = sb->sequences[i];

                for(j = 0; j < s->seq_len ;j++){
                        fhmm->e[s->label[j]][s->seq[j]]++;
                        fhmm_log->e[s->label[j]][s->seq[j]]++;

                }
        }
        for(i = 0; i < fhmm_log->K;i++){
                sum = 0;
                for(j = 0 ; j < fhmm_log->L;j++){
                        sum += fhmm_log->e[i][j];

                }

                for(j = 0 ; j < fhmm_log->L;j++){
                        fhmm_log->e[i][j] =  prob2scaledprob(fhmm_log->e[i][j] / sum);
                        fhmm->e[i][j] =  fhmm->e[i][j] / sum;

                }

        }


        MMALLOC(rel_entropy, sizeof(float) * fhmm->K);

        /* Step one: calculate relative entropy for each state */
        RUN(calc_per_state_rel_entrophy(fhmm, rel_entropy));



        /* Step two: label sequences with rel entropy perstate */
        for(i = 0; i < sb->num_seq;i++){
                s = sb->sequences[i];
                for(j = 0; j < s->seq_len;j++){
                        s->u[j] = rel_entropy[s->label[j]];
                        //fprintf(stdout,"%d %d %f\n",j, s->seq[j],s->u[j]);
                }
                //fprintf(stdout,"\n");
        }

        RUNP(sa = build_sa(sb));
        /* Now I can set up */
        /* read in finite hmm */

        m = init_motif_list(sb->num_seq);


        for(i = 0; i < sb->num_seq;i++){
                s = sb->sequences[i];
                t = init_motif_list(sb->num_seq);

                p = init_paraclu_cluster();
                RUN(max_score_segment(sb->sequences[i]->u, 0, sb->sequences[i]->seq_len,6,0.1,t));

                /*for(j = 0; j < t->num_items;j++){
                        p = t->plist[j];
                        fprintf(stderr,"CLUSTER:%d:\t%d\t%d\t%d\t%f\t%f\t%f\n",j,p->start, p->stop,p->len,p->total,p->min_density,p->max_density);
                        }*/

                /* filter for max len */
                for(j = 0; j < t->num_items;j++){
                        p = t->plist[j];
                        if(p->total > 0){
                                if(p->len > param->maxlen){
                                        p->total = -1; /* marker to not include sequence */
                                }
                        }

                }
                /* filter for min increase in d */
                for(j = 0; j < t->num_items;j++){
                        p = t->plist[j];
                        if(p->total > 0){
                                if(p->min_density * param->min_d_increase > p->max_density){
                                        p->total = -2; /* marker to not include sequence */
                                }
                        }

                }

                /* filter clusters contained within clusters */
                old_end = -1;
                for(j = 0; j < t->num_items;j++){
                        p = t->plist[j];
                        if(p->total > 0){
                                if(p->stop <= old_end ){
                                        //   p->total = -3; /* marker to not include sequence */
                                }
                                old_end = p->stop;
                        }
                }
                /*fprintf(stdout,"\n");
                for(j = 0; j < t->num_items;j++){
                        p = t->plist[j];
                        fprintf(stderr,"CLUSTER:%d:\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n",j,p->start, p->stop,p->len,p->total,p->min_density,p->max_density,p->max_density /p->min_density  );
                }
                fprintf(stdout,"\n");*/
                //exit(0);

                /*if(p->start == -1){
                        free_paraclu_cluster(p);
                        break;
                }
                MMALLOC(p->state_sequence, sizeof(int) * p->stop-p->start);
                for(c = 0; c < p->stop-p->start;c++){
                        p->state_sequence[c] = s->label[p->start + c];
                }
                fprintf(stderr,"CLUSTER:%d:%d	%d	%d	%d	%f %f %f\n",i,j,p->start, p->stop,p->stop-p->start, p->kl_divergence,p->max_d ,p->kl_divergence / (float)(p->stop-p->start) );

                for(c = p->start;c < p->stop;c++){
                        sb->sequences[i]->u[c] = 0.0f;
                }*/
                /* print cluster  */
                for(j = 0; j < t->num_items;j++){
                        p = t->plist[j];
                        if(p->total > 0){
                                MMALLOC(p->state_sequence, sizeof(int) * p->stop-p->start);
                                for(c = 0; c < p->stop-p->start;c++){
                                        p->state_sequence[c] = s->label[p->start + c];
                                }


                                RUN(search_sa(sa,p->state_sequence, p->len,&p->start_in_sa,&p->end_in_sa));
                                RUN(add_present_in_seq(p,sa, sb->num_seq));
                                RUN(calculate_log_likelihood(fhmm_log,sb,p,sa));

                                //fprintf(stderr,"CLUSTER:%d	%d	%d	%d %f s:%d e:%d n:%d llr:%f in store:%d\n",i,p->start, p->stop,p->len,p->total, p->start_in_sa, p->end_in_sa,p->end_in_sa- p->start_in_sa,p->log_likelihood,m->num_items);
                                RUN(insert_motif(m,p,fhmm));
                                //fprintf(stderr,"%d\n", m->num_items);
                                t->plist[j] = NULL;

                        }
                }
                //exit(0);
                //free_paraclu_cluster(p);
                free_motif_list(t);

        }



        RUN(calculate_relative_entropy(m, fhmm));


        qsort(m->plist , m->num_items, sizeof(struct paraclu_cluster*), sort_paraclu_cluster_based_on_likelihood);

        for(i = 0 ; i < MACRO_MIN(1000000, m->num_items-1);i++){
                p = m->plist[i];

                if(p->num_present_in_seq < 50){
                        free_paraclu_cluster(p);
                        m->plist[i] = 0;
                }else{
                fprintf(stderr,"CLUSTER:%d	%d %f  llr:%f   %d in %d seq\n",i,p->len,p->total,p->log_likelihood,p->count, p->num_present_in_seq);
                }
        }
        RUN(write_para_motif_data_for_plotting(param->out, m, fhmm));
        free_fhmm(fhmm);
        free_fhmm(fhmm_log);
        free_motif_list(m);
        //run_paraclu_clustering(sb, "GAGA");
        /* clean up */
        free_ihmm_sequences(sb);
        free_sa(sa);
        MFREE(rel_entropy);
        return OK;
ERROR:
        if(fhmm){
                free_fhmm(fhmm);
        }
        if(sb){
                free_ihmm_sequences(sb);
        }
        MFREE(rel_entropy);
        return FAIL;
}

int add_present_in_seq(struct paraclu_cluster* motif,struct sa* sa, int numseq)
{
        int i,j;

        ASSERT(motif!= NULL, "No motif");
        ASSERT(sa != NULL, "No suffix array");
        motif->present_in_seq = NULL;
        motif->num_present_in_seq = 0;
        MMALLOC(motif->present_in_seq, sizeof(int) * numseq);
        motif->num_seq = numseq;

        for(i = 0; i < numseq;i++){
                motif->present_in_seq[i] = 0;
        }
        for(i = motif->start_in_sa; i < motif->end_in_sa;i++){
                motif->present_in_seq[ sa->lcs[i]->seq_num]++;
        }

        for(i = 0; i < numseq;i++){
                if(motif->present_in_seq[i]){
                        motif->num_present_in_seq++;
                }
        }


        return OK;
ERROR:
        return FAIL;

}
int calculate_log_likelihood(struct fhmm*  fhmm,struct seq_buffer* sb, struct paraclu_cluster* motif,struct sa* sa)
{
        int i,j,c;
        int w;
        int state;
        int occur;
        float a,b,cc;
        float sum;
        float lambda_0;         /* background prior */

        float lambda_1;         /* motif prior */

        float log_likelihood_0;

        float log_likelihood_1;

        struct ihmm_sequence* s = NULL;

        ASSERT(fhmm != NULL, "No hmm");
        ASSERT(sb != NULL, "No sequences");
        ASSERT(sb->num_seq != 0, " No sequences");
        ASSERT(motif != NULL, "No motif");

        w = motif->len;

        /* remove *some* overlapping motifs */

        RUN(find_overlapping_hits(sa, motif->start_in_sa,motif->end_in_sa,w));


        occur =  motif->end_in_sa - motif->start_in_sa;
        //fprintf(stdout,"%d occur \t",occur);
        occur = 0;
        for(i = motif->start_in_sa; i < motif->end_in_sa;i++){
                if(sa->lcs[i]->lcp){
                        occur++;
                }
        }
        //fprintf(stdout,"%d occur \n",occur);
        //w = 1;
        /* set z_ij */
        c= 0;
        for( i = 0; i < sb->num_seq;i++){
                c += sb->sequences[i]->seq_len - w;
        }
        //c /= 4;

        lambda_0 = (float)(c - occur) / (float)c;
        lambda_1 = (float)(occur) / (float)c;


        lambda_0 = prob2scaledprob(lambda_0);
        lambda_1 = prob2scaledprob(lambda_1);
         cc = prob2scaledprob(1.0f);
        for(i = motif->start_in_sa; i < motif->end_in_sa;i++){
                if(sa->lcs[i]->lcp){
                        s = sb->sequences[sa->lcs[i]->seq_num];
                        j = sa->lcs[i]->pos;


                        a = prob2scaledprob(1.0);
                        for(c = 0; c < w;c++){
                                state = motif->state_sequence[c];
                                a = a + fhmm->e[state][s->seq[j+c]];
                        }

                        b = prob2scaledprob(1.0);
                        for(c = 0; c < w;c++){
                                b = b + fhmm->background[s->seq[j+c]];
                        }
                        //fprintf(stdout,"%d %d %f  %f\n", i,j, (a - b) + (lambda_1 - lambda_0), exp((a - b) + (lambda_1 - lambda_0))/ (1.0f + exp((a - b) + (lambda_1 - lambda_0))));
                        cc += (a - b) + (lambda_1 - lambda_0);
                }

        }

        motif->log_likelihood  = cc;
        return OK;
ERROR:
        return FAIL;
}

/* Here I mark overlapping motifs using the lcp variable - this is a
 * hack there should be a 'use' variable in the lcs struct */

int find_overlapping_hits(struct sa* sa, int start, int stop, int w)
{
        int i,j;
        struct lcs* a;
        struct lcs* b;

        ASSERT(sa != NULL, "No suffix array");

        for(i = start; i < stop;i++){
                sa->lcs[i]->lcp = 1;
        }

        for(i = start; i < stop;i++){
                a = sa->lcs[i];
                for(j = i+1; j < stop;j++){
                        b = sa->lcs[j];
                        //if(a->lcp){
                                if(a->seq_num == b->seq_num){
                                        //b->lcp = 0;
                                        if(MACRO_MIN(a->pos, b->pos)+ w > MACRO_MAX(a->pos, b->pos)){
                                                        b->lcp = 0;
                                         }
                                }
                                //}

                }
        }
        return OK;
ERROR:
        return FAIL;

}


int calculate_relative_entropy(struct motif_list* m, struct fhmm* fhmm)
{
        int i,j,c,s;
        struct paraclu_cluster* p =  NULL;
        float* background = NULL;
        float x;

        ASSERT(m != NULL, "No rbtree");
        ASSERT(fhmm != NULL, "No fhmm");

        background = fhmm->background;
        for(i = 0; i < m->num_items;i++){
                p = m->plist[i];

                x = 0.0f;
                for(j = 0; j < p->len;j++){
                        s = p->state_sequence[j];

                        for(c = 0; c < fhmm->L;c++){
                                if(fhmm->e[s][c] > 0.0f){
                                        x += fhmm->e[s][c] * log2f(fhmm->e[s][c] / background[c]);
                                }
                        }
                }
                x = x / (float) p->len;
                p->total = x;

        }
        return OK;
ERROR:
        return FAIL;
}


int write_para_motif_data_for_plotting(char* filename,struct motif_list* m, struct fhmm*  fhmm)
{
        char buffer[BUFFER_LEN];
        struct hdf5_data* hdf5_data = NULL;
        struct paraclu_cluster* p =  NULL;
        float** tmp = NULL;
        int i,j,c;



        ASSERT(filename != NULL, "No file name");
        ASSERT(m != NULL, "No motifs");



        /* determine max_len of motif - to alloc frequency matrix */

        RUNP(hdf5_data = hdf5_create());
        snprintf(buffer, BUFFER_LEN, "%s",PACKAGE_NAME );
        hdf5_add_attribute(hdf5_data, "Program", buffer, 0, 0.0f, HDF5GLUE_CHAR);
        snprintf(buffer, BUFFER_LEN, "%s", PACKAGE_VERSION);
        hdf5_add_attribute(hdf5_data, "Version", buffer, 0, 0.0f, HDF5GLUE_CHAR);


        RUN(hdf5_create_file(filename,hdf5_data));

        hdf5_write_attributes(hdf5_data, hdf5_data->file);


        RUN(hdf5_create_group("MotifData",hdf5_data));


        for(i = 0 ; i < m->num_items;i++){
                p = m->plist[i];
                if(p){
                RUNP(tmp = malloc_2d_float(tmp, p->len, fhmm->L, 0.0));


                for(j = 0 ; j < p->len;j++){
                        for(c = 0; c < fhmm->L;c++){
                                tmp[j][c] = fhmm->e[p->state_sequence[j]][c];
                        }
                }


                hdf5_data->rank = 2;
                hdf5_data->dim[0] = p->len;
                hdf5_data->dim[1] = fhmm->L;
                hdf5_data->chunk_dim[0] = p->len;
                hdf5_data->chunk_dim[1] = fhmm->L;
                hdf5_data->native_type = H5T_NATIVE_FLOAT;
                snprintf(buffer, BUFFER_LEN, "M%03d;%4.2f", i+1,p->log_likelihood );
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



int calc_per_state_rel_entrophy(struct fhmm* fhmm, float* rel_entropy)
{
        int i,j;
        float* background = NULL;
        float x;
        ASSERT(rel_entropy != NULL, "No entropy array malloced");
        ASSERT(fhmm != NULL, "No fhmm");
        background = fhmm->background;
        for(i = 0; i < fhmm->K;i++){
                x = 0.0f;
                for(j = 0; j < fhmm->L;j++){
                        if(fhmm->e[i][j] > 0.0f){
                                x += fhmm->e[i][j] * log2f(fhmm->e[i][j] / background[j]);
                        }
                }
                rel_entropy[i] = x;
                /*fprintf(stdout,"State %d: %f\t",i, x);
                for(j = 0; j < fhmm->L;j++){
                        fprintf(stdout,"%f ", fhmm->e[i][j]);
                }
                fprintf(stdout,"\n");*/
        }
        return OK;
ERROR:
        return FAIL;
}

int insert_motif(struct motif_list* m, struct paraclu_cluster* p,struct fhmm* fhmm)
{
        int i,j;
        int new = 1;
        float kl;
        for(i = 0; i < m->num_items;i++){

                /* compare motifs - keep longer one... */
                RUN(compare_motif(fhmm, m->plist[i], p, &kl));
                //fprintf(stdout,"%d %f %f %f\n",i, kl,kl / MACRO_MIN(p->len, m->plist[i]->len),m->plist[i]->log_likelihood);
                if(kl / MACRO_MIN(p->len, m->plist[i]->len)  < 1.0f){
                        new = 0;

                        if(p->log_likelihood > m->plist[i]->log_likelihood){
                                //                fprintf(stdout,"replace existing cluster\n");
                                merge_motif(p,m->plist[i]);
                                free_paraclu_cluster(m->plist[i]);
                                m->plist[i] = p;
                        }else{
                                merge_motif(m->plist[i],p);
                                free_paraclu_cluster(p);

                        }
                        break;
                }

        }
        //fprintf(stdout,"Insert? %d\n",new);
        if(new){                /* still */
                m->plist[m->num_items] = p;
                m->num_items++;
                if(m->num_items == m->alloc_items){
                        RUN(extend_motif_list(m));
                }
        }
        qsort(m->plist , m->num_items, sizeof(struct paraclu_cluster*), sort_paraclu_cluster_based_on_likelihood);
        return OK;
ERROR:
        return FAIL;
}

int merge_motif(struct paraclu_cluster* a,struct paraclu_cluster* b)
{
        int i;

        ASSERT(a != NULL, "No a motif");
        ASSERT(b != NULL, "No b motif");
        a->count = a->count + b->count;
        for(i = 0; i < a->num_seq;i++){
                a->present_in_seq[i] +=  b->present_in_seq[i];
        }
        a->num_present_in_seq = 0;
        for(i = 0; i < a->num_seq;i++){
                if(a->present_in_seq[i]){
                        a->num_present_in_seq++;
                }
        }
        return OK;
ERROR:
        return FAIL;
}

/* just compare - decision what to do later */
int compare_motif(struct fhmm* fhmm, struct paraclu_cluster* a,struct paraclu_cluster* b, float* kl_div)
{
        float x;
        float y;
        int i,j,c;

        struct paraclu_cluster* temp = NULL;

        ASSERT(a != NULL, "No a");
        ASSERT(b != NULL, "No b");

        if(b->len > a->len){
                temp = a;
                a = b;
                b = temp;
        }

        *kl_div = 10000000;
        /* outer loop - slide shorter motif across longer one.. */

        for(i = 0; i <= a->len - b->len;i++){
                x = 0.0f;
                y = 0.0f;
                for(j = 0; j < b->len;j++){
                        for(c = 0; c < fhmm->L;c++){
                                if(fhmm->e[a->state_sequence[i+j]][c] > 0.0f && fhmm->e[b->state_sequence[j]][c] > 0.0f ){
                                        x += fhmm->e[a->state_sequence[i+j]][c] * log2f(fhmm->e[a->state_sequence[i+j]][c] / fhmm->e[b->state_sequence[j]][c]);
                                        y += fhmm->e[b->state_sequence[j]][c] * log2f(fhmm->e[b->state_sequence[j]][c] / fhmm->e[a->state_sequence[i+j]][c]);
                                }
                        }
                }
                if(x < *kl_div){
                        *kl_div = x;
                }
                if(y < *kl_div){
                        *kl_div = y;
                }
        }
        return OK;
ERROR:
        return FAIL;
}

struct motif_list* init_motif_list(int initial_len)
{
        struct motif_list* m = NULL;
        int i;
        ASSERT(initial_len > 0, "No length");

        MMALLOC(m, sizeof(struct motif_list));
        m->plist = NULL;
        m->num_items = 0;
        m->alloc_items = initial_len;
        MMALLOC(m->plist, sizeof(struct paraclu_cluster*)* m->alloc_items);
        for(i = 0; i < m->alloc_items;i++){
                m->plist[i] = NULL;
        }
        return m;
ERROR:
        return NULL;
}

int extend_motif_list(struct motif_list* m )
{
        int i,old_len;

        ASSERT(m != NULL, "No list");
        old_len = m->alloc_items;
        m->alloc_items = m->alloc_items << 1;
        MREALLOC(m->plist, sizeof(struct paraclu_cluster*)* m->alloc_items);
        for(i = old_len; i < m->alloc_items;i++){
                m->plist[i] = NULL;
        }
        return OK;
ERROR:

        return FAIL;
}

void free_motif_list(struct motif_list* m)
{
        int i;
        if(m){
                for(i = 0; i < m->num_items;i++){
                        if(m->plist[i]){
                                free_paraclu_cluster(m->plist[i]);
                        }
                }
                MFREE(m->plist);
                MFREE(m);
        }

}

struct paraclu_cluster* init_paraclu_cluster(void)
{
        struct paraclu_cluster* p = NULL;

        MMALLOC(p, sizeof(struct paraclu_cluster));
        p->state_sequence = NULL;
        p->seq_id = -1;
        p->start = -1;
        p->stop = -1;
        p->total = -1;
        p->min_density = 0.0f;
        p->max_density = 0.0f;
        p->log_likelihood = 0.0f;
        p->count = 1;
        p->start_in_sa = -1;
        p->end_in_sa = -1;
        p->present_in_seq = NULL;
        p->num_present_in_seq = 0;
        p->num_seq = 0;
        return p;
ERROR:
        return NULL;
}

void free_paraclu_cluster(struct paraclu_cluster* p)
{
        if(p){
                if(p->present_in_seq){
                        MFREE(p->present_in_seq);
                }
                if(p->state_sequence){
                        MFREE(p->state_sequence);
                }
                MFREE(p);
        }
}

int max_score_segment(float* x , int start ,int end, int min_len, float min_density,struct motif_list* m)
{
        struct paraclu_cluster* p;
        float max_density;
        float new_min_density;
        float min_prefix;
        float min_suffix;
        float total;
        int min_prefix_pos;
        int min_suffix_pos;
        int mid;
        //fprintf(stderr,"Looking at %d	%d\n",start,end);
        if (start == end){
                return OK;
        }
        total = weakestPrefix(x,start,end,&min_prefix_pos, &min_prefix);

        if (total < 0.0){
                return OK;
        }
        weakestSuffix(x,start,end,&min_suffix_pos,  &min_suffix);

        max_density =  (min_prefix < min_suffix) ? min_prefix : min_suffix;


        if (max_density > min_density &&  end-start >= min_len){
                //       if(max_density /  min_density >= p->max_d){
                p = init_paraclu_cluster();
                p->start = start;
                p->stop = end;
                p->len = end - start;
                p->total = total;
                p->min_density = min_density;
                p->max_density = max_density;
                m->plist[m->num_items] = p;
                m->num_items++;
                if(m->num_items == m->alloc_items){
                        RUN(extend_motif_list(m));
                }
        }

        if (max_density < 1e100) {
                mid = (min_prefix < min_suffix) ?  min_prefix_pos: min_suffix_pos;
                new_min_density =  (max_density > min_density) ? max_density : min_density;
                //Ci mid = (minPrefixDensity < minSuffixDensity) ? minPrefix : minSuffix;
                //fprintf(stderr,"test:%f	%f\n",new_min_density,max_density);
                RUN(max_score_segment(x ,start , mid , min_len,new_min_density, m));
                //fprintf(stderr,"test:%f	%f\n",new_min_density,max_density);
                //writeClusters(beg, mid, newMinDensity);
                RUN(max_score_segment(x ,mid , end , min_len,new_min_density, m));
                //writeClusters(mid, end, newMinDensity);
        }
        return OK;
ERROR:
        return FAIL;
}

float weakestPrefix(float* x , int start ,int end, int* min_prefix_pos, float* min_prefix)
{
        int origin = start;
        float density  = 0.0;
        //minPrefix = beg;
        *min_prefix = 1e100;
        float totalValue = x[start];
        density = totalValue / (float)(start  - origin);
        if (density < *min_prefix) {
                *min_prefix_pos = start;
                *min_prefix = density;
        }
        ++start;

        while (start < end) {
                //fprintf(stderr,"%d	%d	%f\n",origin, start,  totalValue / (float)(start  - origin));
                density = totalValue / (float)(start  - origin);
                if (density < *min_prefix) {
                        *min_prefix_pos = start;
                        *min_prefix = density;
                }
                totalValue += x[start];
                start++;
        }
        //fprintf(stderr,"%f	return\n",totalValue);
        return totalValue;
}

void weakestSuffix(float* x , int start ,int end, int* min_suffix_pos, float* min_suffix)
{

        --end;
        int origin = end;
        float density  = 0.0;
        //minSuffix = end + 1;
        *min_suffix = 1e100;
        float totalValue = x[end];

        while (end > start) {
                --end;
                //fprintf(stderr,"%d	%d	%f\n",origin, end,  totalValue / (float)( origin-end));
                density = totalValue / (origin - end);
                if (density < *min_suffix) {
                        *min_suffix_pos = end + 1;
                        *min_suffix = density;
                }
                totalValue += x[end];//end->value;
        }
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
        const char usage[] = " -m <h5 model> -out <h5 out>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");
        /*fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--len","Minimum pattern length." ,"[10]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--frac","Minimum fraction of seq containing pattern." ,"[0.5]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-d","Resolution." ,"[100]"  );*/
        return OK;
}


float* shuffle_arr_r(float* arr,int n,unsigned int* seed)
{
        int i,j;
        float tmp;
        for (i = 0; i < n - 1; i++) {
                j = i +  (int) (rand_r(seed) % (int) (n-i));
                tmp = arr[j];
                arr[j] = arr[i];
                arr[i] = tmp;
        }
        return arr;
}



int sort_paraclu_cluster_based_on_likelihood(const void *a, const void *b)
{
    struct paraclu_cluster* const *one = a;
    struct paraclu_cluster* const *two = b;

    if((*one)->log_likelihood <  (*two)->log_likelihood){
            return 1;
    }else{
            return -1;
    }
}
