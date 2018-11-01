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
        unsigned long seed;
        int maxlen;
};
/* global rnd state */
rk_state rndstate;

static int run_infoclust(struct parameters* param);


int clean_sequence_labelling(struct fhmm* fhmm_log,struct fhmm* fhmm, struct seq_buffer* sb);


int calc_per_state_rel_entrophy(struct fhmm* fhmm, float* rel_entropy);

int calculate_relative_entropy(struct motif_list* m, struct fhmm* fhmm);

int calculate_log_likelihood(struct fhmm* fhmm,struct seq_buffer* sb, struct paraclu_cluster* motif,struct sa* sa,struct seq_buffer* sb_temp);

int calculate_log_likelihood_based_on_hits(struct fhmm*  fhmm,struct seq_buffer* sb, struct paraclu_cluster* motif,struct sa* sa, struct seq_buffer* sb_temp);

int add_present_in_seq(struct paraclu_cluster* motif,struct sa* sa, int numseq);

int add_motif_count_matrix(struct seq_buffer* sb, struct paraclu_cluster* motif,struct sa* sa);
int find_overlapping_hits(struct sa* sa, int start, int stop, int w);

int insert_motif(struct motif_list*m, struct paraclu_cluster* p,struct fhmm* fhmm,struct seq_buffer *sb);

static int randomize_freq_matrix(float** m,int len, int L);


int compare_motif(struct fhmm* fhmm, struct paraclu_cluster* a,struct paraclu_cluster* b, float* kl_div);
int align_motif(struct fhmm* fhmm, struct paraclu_cluster* a,struct paraclu_cluster* b,float* score);
int traceback_and_merge_motifs(struct fhmm* fhmm,struct paraclu_cluster* a,struct paraclu_cluster* b);
static int motif_dyn_programming(struct fhmm* fhmm,float** e_a, float** e_b, int len_a, int len_b,int min_aln_len);
int pick_state(float* a, float* b, int len, int* pick);
int merge_motif(struct paraclu_cluster* a,struct paraclu_cluster* b);
int write_para_motif_data_for_plotting(char* filename,struct motif_list* m, struct fhmm*  fhmm);
int write_meme_output(char* filename,struct motif_list* m, struct fhmm*  fhmm);
int compare_hit_positions(struct hit** list_a, struct hit** list_b, int len_a,int len_b, int w,float* jac);
static int add_motif_count_matrix_based_on_hits(struct seq_buffer* sb, struct paraclu_cluster* motif);
int add_hit_locations(struct sa* sa, struct paraclu_cluster* p);
static int remove_overlapping_hits( struct paraclu_cluster* p);
struct motif_list* init_motif_list(int initial_len);
int extend_motif_list(struct motif_list* m);
void free_motif_list(struct motif_list* m);

struct paraclu_cluster* init_paraclu_cluster(void);
void free_paraclu_cluster(struct paraclu_cluster* p);

int sort_paraclu_cluster_based_on_likelihood(const void *a, const void *b);
int sort_paraclu_cluster_based_on_hits(const void *a, const void *b);
int float_cmp(const void *a, const void *b);
int max_score_segment(float* x , int start ,int end, int min_len,float min_density,struct motif_list* m);
float weakestPrefix(float* x , int start ,int end, int* min_prefix_pos, float* min_prefix);
void weakestSuffix(float* x , int start ,int end, int* min_suffix_pos, float* min_suffix);
int free_parameters(struct parameters* param);
int print_help(char **argv);

int shuffle_arr_r(int* arr,int n,unsigned int* seed);

static int sort_hit_positions(const void *a, const void *b);
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
        param->maxlen = 30;
        param->min_d_increase = 1.00;
        param->seed = 0;
        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"out",required_argument,0,'o'},
                        {"seed",required_argument,0,'s'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"m:o:s:",long_options, &option_index);

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
                case 's':
                        param->seed = atoi(optarg);
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

        if(param->seed){
                rk_seed(param->seed, &rndstate);
        }else{
                rk_randomseed(&rndstate);
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
        struct seq_buffer* sb_temp = NULL;
        struct ihmm_sequence* s = NULL;
        struct fhmm* fhmm = NULL;
        struct fhmm* fhmm_log = NULL;
        struct paraclu_cluster* p =  NULL;
        struct motif_list* m = NULL;
        struct motif_list* t = NULL;
        struct sa* sa = NULL;

        float* rel_entropy = NULL;

        int total_len = 0;
        int i,j,c;
        int old_end;
        ASSERT(param != NULL, "No parameters.");

        /* read in sequences */
        RUNP(sb = get_sequences_from_hdf5_model(param->input));
        ASSERT(sb != NULL, "No sequence Buffer");

        /* read sequences in again but this time only as temp storage... */
        RUNP(sb_temp = get_sequences_from_hdf5_model(param->input));
        ASSERT(sb_temp != NULL, "No sequence Buffer");

        RUNP(fhmm_log = alloc_fhmm());
        /* get HMM parameters  */
        RUN(read_hmm_parameters(fhmm_log,param->input));

        RUNP(fhmm = alloc_fhmm());

        /* get HMM parameters  */
        RUN(read_hmm_parameters(fhmm,param->input));

        RUN(alloc_dyn_matrices(fhmm));
        RUN(realloc_dyn_matrices(fhmm, param->maxlen));

        RUN(setup_model(fhmm_log));


        RUN(clean_sequence_labelling(fhmm_log,fhmm, sb));


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
                total_len += s->seq_len;
                //fprintf(stdout,"\n");
        }

        RUNP(sa = build_sa(sb));
        /* Now I can set up */

        int start, stop;
        int occ_cutoff = (int) sqrtf((float) sb->num_seq);
        int target_len = 12;

        m = init_motif_list(sb->num_seq);
        p = init_paraclu_cluster();
        MMALLOC(p->state_sequence, sizeof(int) * 50);
        for(target_len = 25; target_len >= 6;target_len--){
                LOG_MSG("Targeting %d.",target_len);
                start = 0;
                stop = -1;


                while(1){
                        int go = 1;
                        for(i = 0; i < target_len;i++){
                                if(sa->lcs[start]->str[i] == -1){
                                        start++;
                                        go = 0;
                                        break;
                                }
                        }

                        if(go){
                                RUN(search_sa(sa,sa->lcs[start]->str, target_len,&start,&stop));
                                if(stop- start >= occ_cutoff){
                                        p->start_in_sa = start;
                                        p->end_in_sa = stop;
                                        p->len = target_len;
                                        if(p->present_in_seq){
                                                MFREE(p->present_in_seq);
                                                p->present_in_seq = NULL;
                                                p->num_present_in_seq = 0;
                                        }
                                        RUN(add_present_in_seq(p,sa, sb->num_seq));

                                        for(i = 0; i < target_len;i++){
                                                p->state_sequence[i] =sa->lcs[start]->str[i];
                                        }
                                        RUN(calculate_log_likelihood(fhmm_log,sb,p,sa,sb_temp));
                                        float measure = (float) p->num_present_in_seq/ (float)(stop-start) * (float) p->num_present_in_seq;
                                        if(p->log_likelihood > 1.0 && measure >= occ_cutoff){

                                                RUN(add_hit_locations(sa,p));
                                                RUN(add_motif_count_matrix_based_on_hits(sb, p));

                                                fprintf(stdout,"occ:%d %f in %d seq\n", stop-start,p->log_likelihood,p->num_present_in_seq);




                                                RUN(insert_motif(m,p,fhmm,sb));
                                                p = NULL;
                                                p = init_paraclu_cluster();
                                                MMALLOC(p->state_sequence, sizeof(int) * 50);
                                        }
                                }
                                start = stop;
                                if(start == sa->len-1){
                                        break;
                                }
                        }
                }
        }


        RUN(write_para_motif_data_for_plotting(param->out, m, fhmm));

        RUN(write_meme_output("motif.meme",m,fhmm));

        free_fhmm(fhmm);
        free_fhmm(fhmm_log);
        free_motif_list(m);

        exit(0);
        /* read in finite hmm */



        for(i =0 ; i < sb->num_seq;i++){
                //LOG_MSG("Working on seq%d",i);
                s = sb->sequences[i];
                t = init_motif_list(sb->num_seq);

                p = init_paraclu_cluster();
                RUN(max_score_segment(sb->sequences[i]->u, 0, sb->sequences[i]->seq_len,8,0.01,t));

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
                                        //       p->total = -3; /* marker to not include sequence */
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
                /*for(j = 0; j < t->num_items;j++){
                        p = t->plist[j];
                        MMALLOC(p->state_sequence, sizeof(int) * p->stop-p->start);
                        for(c = 0; c < p->stop-p->start;c++){
                                p->state_sequence[c] = s->label[p->start + c];
                        }
                        RUN(search_sa(sa,p->state_sequence, p->len,&p->start_in_sa,&p->end_in_sa));
                        RUN(add_present_in_seq(p,sa, sb->num_seq));
                        RUN(calculate_log_likelihood(fhmm_log,sb,p,sa,sb_temp));
                        float t = sqrtf((float) sb->num_seq);
                        if(t <= p->end_in_sa- p->start_in_sa){
                                fprintf(stderr,"CLUSTER:%d\t%d\t%f\tinseq:%d s:%d e:%d n:%d llr:%f in store:%d %f\n",i,p->len,p->total, p->num_present_in_seq, p->start_in_sa, p->end_in_sa,p->end_in_sa- p->start_in_sa,p->log_likelihood,m->num_items,t);
                        }
                        }*/

                /* print cluster  */
                for(j = 0; j < t->num_items;j++){
                        p = t->plist[j];
                        if(p->total > 0 ){
                                MMALLOC(p->state_sequence, sizeof(int) * p->stop-p->start);
                                for(c = 0; c < p->stop-p->start;c++){
                                        p->state_sequence[c] = s->label[p->start + c];
                                }

                                RUN(search_sa(sa,p->state_sequence, p->len,&p->start_in_sa,&p->end_in_sa));
                                //RUN(calculate_log_likelihood_based_on_hits(fhmm_log,sb,p,sa,sb_temp));
                                //fprintf(stdout,"%f ", p->log_likelihood);
                                RUN(calculate_log_likelihood(fhmm_log,sb,p,sa,sb_temp));
                                
                                //fprintf(stdout,"%f len:%d\n", p->log_likelihood,p->len);
                                if(p->log_likelihood >= 2.0){
                                        //RUN(add_motif_count_matrix(sb, p,sa));

                                        RUN(add_present_in_seq(p,sa, sb->num_seq));
                                        RUN(add_hit_locations(sa,p));
                                        RUN(add_motif_count_matrix_based_on_hits(sb, p));
                                        //RUN(add_hit_locations(sa,p));
                                        //RUN(add_motif_count_matrix_based_on_hits(sb, p));

                                        fprintf(stderr,"CLUSTER:%d\t%d-%d\t%d\t%f\tinseq:%d\ts:%d\te:%d\tn:%d llr:%f in store:%d\n",i,p->start,p->stop, p->len,p->total, p->num_present_in_seq, p->start_in_sa, p->end_in_sa,p->end_in_sa- p->start_in_sa,p->log_likelihood,m->num_items);

                                        RUN(insert_motif(m,p,fhmm,sb));
                                        t->plist[j] = NULL;
                                        }else{
                                                free_paraclu_cluster(p);
                                                t->plist[j] = NULL;
                                        }
                        }
                }
                free_motif_list(t);
        }

        RUN(calculate_relative_entropy(m, fhmm));

        qsort(m->plist , m->num_items, sizeof(struct paraclu_cluster*), sort_paraclu_cluster_based_on_likelihood);
        c = 0;
        for(i = 0 ; i < MACRO_MIN(1000000, m->num_items);i++){
                p = m->plist[i];

                if(p->num_present_in_seq < (int)sqrtf((float) sb->num_seq)){
                        free_paraclu_cluster(p);
                        m->plist[i] = 0;
                }else{
                        fprintf(stderr,"CLUSTER:%d	%d %f  llr:%f   %d in %d seq\n",i,p->len,p->total,p->log_likelihood,p->count, p->num_present_in_seq);


                        RUN(add_motif_count_matrix_based_on_hits(sb,p));
                        //m->plist[c] = p;
                        //c++;
                        // m->plist[i] = 0;
                }
        }



        RUN(write_para_motif_data_for_plotting(param->out, m, fhmm));

        RUN(write_meme_output("motif.meme",m,fhmm));

        free_fhmm(fhmm);
        free_fhmm(fhmm_log);
        free_motif_list(m);

        //run_paraclu_clustering(sb, "GAGA");
        /* clean up */
        free_ihmm_sequences(sb_temp);
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
        if(sb_temp){
                free_ihmm_sequences(sb_temp);
        }

        MFREE(rel_entropy);
        return FAIL;
}



int clean_sequence_labelling(struct fhmm* fhmm_log,struct fhmm* fhmm, struct seq_buffer* sb)
{
        struct ihmm_sequence* s = NULL;
        int* array_rename = NULL;

        int i,j,c;
        float sum;

        ASSERT(fhmm_log != NULL, "No finite hmm model.");
        ASSERT(fhmm != NULL, "No finite hmm model.");
        /* Step 1: run full posterior decoding to label sequences */
        RUN(run_label_sequences(fhmm_log,sb, 8));
        /* re-estimate emissions */
        for(i = 0; i < fhmm_log->K;i++){
                for(j = 0 ; j < fhmm_log->L;j++){
                        fhmm_log->e[i][j] = 1.0 * fhmm->background[j]; /* pseudocount */
                        fhmm->e[i][j] = 1.0 * fhmm->background[j];
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


        /* Step 2: merge states with very similar emission probabilities */
        MMALLOC(array_rename, sizeof(int) *fhmm->K);
        for(i = 0; i < fhmm->K;i++){
                array_rename[i] = i;
        }
        for(i = 2; i < fhmm->K-1;i++){
                for(j = i+1; j < fhmm->K-1;j++){
                        float x,y;
                        x = 0.0;
                        y = 0.0;
                        for(c = 0; c < fhmm->L;c++){
                                if(fhmm->e[i][c] > 0.0f && fhmm->e[j][c] > 0.0f){
                                        x += fhmm->e[i][c] * log2f(fhmm->e[i][c] / fhmm->e[j][c]);
                                        y += fhmm->e[j][c] * log2f(fhmm->e[j][c] / fhmm->e[i][c]);
                                }
                        }


                        if((x+y) / 2.0f < 0.01){
                                //fprintf(stdout,"%d - %d : %f\n", i,j, (x+y) / 2.0f);
                                //for(c = 0; c < fhmm->L;c++){
                                //       fprintf(stdout,"%c %f %f\n","ACGT"[c],  fhmm->e[i][c],fhmm->e[j][c] );
                                        //}
                                array_rename[MACRO_MAX(i, j)] = MACRO_MIN(i,j);
                        }
                }
        }

        /* remove now redundant labels from sequences */
        for(i = 0; i < sb->num_seq;i++){
                for(j = 0; j < sb->sequences[i]->seq_len;j++){
                        sb->sequences[i]->label[j] = array_rename[sb->sequences[i]->label[j]];

                }
        }

        /* We don't want to lose emission counts - therefore we
         * re-count everything based on re-labelled sequences */
        /* re-estimate emissions */
        for(i = 0; i < fhmm_log->K;i++){
                for(j = 0 ; j < fhmm_log->L;j++){
                        fhmm_log->e[i][j] = 1.0 * fhmm->background[j]; /* pseudocount */
                        fhmm->e[i][j] = 1.0 * fhmm->background[j];
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
        MFREE(array_rename);
        return OK;
ERROR:
        return FAIL;
}
/*int align_output_motifs(struct motif_list* m, struct fhmm* fhmm)
{

        return OK;
ERROR:
        return FAIL;
        }*/

int add_hit_locations(struct sa* sa, struct paraclu_cluster* p)
{
        int i,j,c = 0;
        ASSERT(p != NULL, "No cluster");
        p->num_hits = p->end_in_sa - p->start_in_sa;
        MMALLOC(p->hits, sizeof(struct hit*) * p->num_hits);
        for(i = 0; i < p->num_hits;i++){
                p->hits[i] = NULL;
                MMALLOC(p->hits[i], sizeof(struct hit));
        }
        c = 0;
        for(i = p->start_in_sa; i < p->end_in_sa;i++){
                p->hits[c]->pos = sa->lcs[i]->pos;
                p->hits[c]->seq = sa->lcs[i]->seq_num;
                c++;
        }

        qsort(p->hits, p->num_hits, sizeof(struct hit*), sort_hit_positions);
        /*if(p->num_hits > 10){
                for(i = 0; i < p->num_hits;i++){
                        fprintf(stdout,"%d %d\n",p->hits[i]->seq ,p->hits[i]->pos);
                }
                }*/
        //RUN(remove_overlapping_hits(p));
        //qsort(p->hits, p->num_hits, sizeof(struct hit*), sort_hit_positions);

        c = p->num_hits;
        for(j = 0; j < p->num_hits;j++){
                if(p->hits[j]->seq == 1000000){
                        c = j;
                        break;
                }
                //fprintf(stdout,"%d %d\n", a->hits[j]->seq,a->hits[j]->pos);
        }
        //fprintf(stdout,"%d\n",c);

        for(j = c; j < p->num_hits;j++){
                MFREE(p->hits[j]);
                p->hits[j] = NULL;
        }
        p->num_hits = c;

        MREALLOC(p->hits, sizeof(struct hit*) *(p->num_hits));



        return OK;
ERROR:
        return FAIL;
}

int remove_overlapping_hits( struct paraclu_cluster* p)
{
        int i,j;
        struct hit* a;
        struct hit* b;
        if(p->num_hits > 1){
                for(i = 0;i < p->num_hits;i++){
                        a = p->hits[i];
                        for(j = i + 1; j < p->num_hits;j++){
                                b = p->hits[j];
                                if(a->seq == b->seq){
                                        if(a->pos + p->len >= b->pos){
                                                b->pos = 1000000;
                                                b->seq = 1000000;
                                                break;
                                        }
                                }
                        }

                }

        }

        return OK;
}

int add_present_in_seq(struct paraclu_cluster* motif,struct sa* sa, int numseq)
{
        int i;

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

int add_motif_count_matrix_based_on_hits(struct seq_buffer* sb, struct paraclu_cluster* motif)
{
        float beta = 0.1;
        float sum = 0.0;
        int i,j;
        int c;
        int seq;
        int pos;

        RUNP(motif->count_matrix = malloc_2d_float(motif->count_matrix, motif->len, sb->L, 0.0f));
        RUNP(motif->freq_matrix = malloc_2d_float(motif->freq_matrix, motif->len, sb->L, 0.0f));
        /* add dirichlet priors based on background frequencies */
        for(i =0; i < motif->len;i++){
                sum = 0.0;
                for(j = 0; j < sb->L;j++){
                        motif->count_matrix[i][j] =sb->background[j] * beta;// rk_gamma(&rndstate, sb->background[j], 1.0);
                        sum += motif->count_matrix[i][j];
                }
                for(j = 0; j < sb->L;j++){
                        motif->count_matrix[i][j] = beta * motif->count_matrix[i][j]/sum;
                        //fprintf(stdout,"%c:%0.3f ","ACGT"[j],prob2scaledprob(motif->matrix[i][j]));
                }
                //fprintf(stdout," SUM:%f\n",sum);
        }

        for(i = 0; i < motif->num_hits;i++){
                seq = motif->hits[i]->seq;
                pos = motif->hits[i]->pos;
                for(c = 0;c < motif->len;c++){
                        if(c + pos < 0){
                                //(stdout,"-");
                        }else if(c + pos >= sb->sequences[seq]->seq_len){
                                //fprintf(stdout,"-");
                        }else{
                                motif->count_matrix[c][sb->sequences[seq]->seq[c + pos]] += 1.0;
                                //fprintf(stdout,"%c","ACGT"[sb->sequences[seq]->seq[c + pos]]);
                        }
                }
                //fprintf(stdout, "\n");
        }
        //fprintf(stdout, "\n");
        for(i =0; i < motif->len;i++){
                sum = 0.0;
                for(j = 0; j < sb->L;j++){
                        sum += motif->count_matrix[i][j];
                }
                for(j = 0; j < sb->L;j++){
                        motif->freq_matrix[i][j] = motif->count_matrix[i][j] / sum;
                        //        fprintf(stdout,"%c:%0.3f ","ACGT"[j],motif->freq_matrix[i][j]);
                }
                //fprintf(stdout," SUM:%f\n",sum);
        }


        //ft->emission[i][j] = rk_gamma(&model->rndstate, model->emission_counts[i][j] + EMISSION_H, 1.0);

        return OK;

ERROR:
        return FAIL;
}

int add_motif_count_matrix(struct seq_buffer* sb, struct paraclu_cluster* motif,struct sa* sa)
{
        struct ihmm_sequence* s = NULL;
        int i,j;
        int offset;
        float sum = 0;

        ASSERT(sb != NULL, "No sequences");
        ASSERT(motif!= NULL, "No motif");
        ASSERT(sa != NULL, "No suffix array");

        RUNP(motif->count_matrix = malloc_2d_float(motif->count_matrix, motif->len, sb->L, 0.0f));

        for(i = motif->start_in_sa; i < motif->end_in_sa;i++){
                s = sb->sequences[sa->lcs[i]->seq_num];
                offset= sa->lcs[i]->pos;
                for(j = 0; j < motif->len;j++){
                        motif->count_matrix[j][s->seq[j+offset]] += 1.0;
                }
        }

        RUNP(motif->freq_matrix = malloc_2d_float(motif->freq_matrix, motif->len, sb->L, 0.0f));
        for(i = 0; i < motif->len;i++){
                sum = 0.0;
                for(j = 0; j < sb->L;j++){
                        sum += motif->count_matrix[i][j];
                }
                for(j = 0; j < sb->L;j++){
                        motif->freq_matrix[i][j] = prob2scaledprob(motif->count_matrix[i][j] / sum);
                }
        }

        return OK;
ERROR:
        return FAIL;
}

int calculate_log_likelihood_based_on_hits(struct fhmm*  fhmm,struct seq_buffer* sb, struct paraclu_cluster* motif,struct sa* sa, struct seq_buffer* sb_temp)
{
        int i,j,c;
        int offset;
        int w;
        float occur;
        float a,b,cc;
        float local_p;

        float lambda_0;         /* background prior */
        float lambda_1;         /* motif prior */

        struct ihmm_sequence* s = NULL;

        ASSERT(fhmm != NULL, "No hmm");
        ASSERT(sb != NULL, "No sequences");
        ASSERT(sb->num_seq != 0, " No sequences");
        ASSERT(motif != NULL, "No motif");

        w = motif->len;

        /* remove *some* overlapping motifs */
        for( i = 0; i < sb_temp->num_seq;i++){
                s = sb_temp->sequences[i];
                for(j = 0; j < s->seq_len;j++){
                        s->u[j] = 0.0;
                }
        }
        for(i = 0; i < motif->num_hits;i++){
                sb_temp->sequences[motif->hits[i]->seq]->u[motif->hits[i]->pos] = 1.0;
        }
        /* smooth function taken from: Unsupervised Learning of
         * Multiple Motifs in Biopolymers Using Expectation
         * Maximization. Elkan,Bailey 1993 */

        /* This helps in reducing the number of times a low complexity
         * motif occurs due to overlapping instances. */

        for(i = 0; i < sb_temp->num_seq;i++){
                s = sb_temp->sequences[i];
                for(offset = 0;offset < w;offset++){
                        for(j = offset; j < s->seq_len - (w*2);j += w){
                                local_p = 0.0;
                                for(c = 0; c < w; c++){
                                        local_p += s->u[j+c];
                                }

                                if(local_p > 1.0f){
                                        for(c = 0; c < w; c++){
                                                if(s->u[j+c]){
                                                        s->u[j+c] /= local_p;
                                                }
                                        }
                                }
                        }
                }
        }

        //

        occur = 0.0f;
        for( i = 0; i < sb_temp->num_seq;i++){
                s = sb_temp->sequences[i];
                for(j = 0; j < s->seq_len-w;j++){
                        occur += s->u[j];
                }
        }
        //occur = sb->num_seq;
        //occur =  motif->end_in_sa - motif->start_in_sa;
        //fprintf(stdout,"%f occur \n",occur);


        //exit(0);
        //w = 1;
        /* set z_ij */
        cc = 0;
        for( i = 0; i < sb->num_seq;i++){
                cc += sb->sequences[i]->seq_len - w;
        }
        //c /= 4;


        //occur = sqrtf(sb_temp->num_seq);

        lambda_1 = (occur) / cc;
        lambda_0 = 1.0 - lambda_1;

        lambda_0 = prob2scaledprob(lambda_0);
        lambda_1 = prob2scaledprob(lambda_1);
        cc = 0.0;


        for(i = 0; i < motif->num_hits;i++){
                s = sb_temp->sequences[motif->hits[i]->seq];
                j = motif->hits[i]->pos;
                a = prob2scaledprob(1.0);
                b = prob2scaledprob(1.0);

                for(c = 0; c < w;c++){
                        //state = motif->state_sequence[c];
                        a = a + prob2scaledprob( motif->freq_matrix[c][s->seq[j+c]]);
                        //fprintf(stdout,"%f %f %f %d \n", a,b,fhmm->e[state][s->seq[j+c]],s->seq[j+c]);
                }
                a = s->u[j] * (a + lambda_1);

                for(c = 0; c < w;c++){
                        b = b + fhmm->background[s->seq[j+c]];
                }
                a += (1.0 - s->u[j])* (b+ lambda_0);
                //fprintf(stdout,"%f %f \n", a,b);
                cc += a - b;
        }
        motif->log_likelihood  = cc;

        //fprintf(stdout,"occ:%d\tllr:%f\n",motif->end_in_sa - motif->start_in_sa,cc);
        return OK;
ERROR:
        return FAIL;


}
int calculate_log_likelihood(struct fhmm*  fhmm,struct seq_buffer* sb, struct paraclu_cluster* motif,struct sa* sa, struct seq_buffer* sb_temp)
{
        int i,j,c;
        int offset;
        int w;
        int state;
        float occur;
        float a,b,cc;
        float local_p;

        float lambda_0;         /* background prior */
        float lambda_1;         /* motif prior */

        struct ihmm_sequence* s = NULL;

        ASSERT(fhmm != NULL, "No hmm");
        ASSERT(sb != NULL, "No sequences");
        ASSERT(sb->num_seq != 0, " No sequences");
        ASSERT(motif != NULL, "No motif");

        w = motif->len;

        /* remove *some* overlapping motifs */
        for( i = 0; i < sb_temp->num_seq;i++){
                s = sb_temp->sequences[i];
                for(j = 0; j < s->seq_len;j++){
                        s->u[j] = 0.0;
                }
        }
        for(i = motif->start_in_sa; i < motif->end_in_sa;i++){
                sb_temp->sequences[sa->lcs[i]->seq_num]->u[sa->lcs[i]->pos] = 1.0;
        }
        /* smooth function taken from: Unsupervised Learning of
         * Multiple Motifs in Biopolymers Using Expectation
         * Maximization. Elkan,Bailey 1993 */

        /* This helps in reducing the number of times a low complexity
         * motif occurs due to overlapping instances. */

        for(i = 0; i < sb_temp->num_seq;i++){
                s = sb_temp->sequences[i];
                for(offset = 0;offset < w;offset++){
                        for(j = offset; j < s->seq_len - (w*2);j += w){
                                local_p = 0.0;
                                for(c = 0; c < w; c++){
                                        local_p += s->u[j+c];
                                }

                                if(local_p > 1.0f){
                                        for(c = 0; c < w; c++){
                                                if(s->u[j+c]){
                                                        s->u[j+c] /= local_p;
                                                }
                                        }
                                }
                        }
                }
        }

        //

        occur = 0.0f;
        for( i = 0; i < sb_temp->num_seq;i++){
                s = sb_temp->sequences[i];
                for(j = 0; j < s->seq_len-w;j++){
                        occur += s->u[j];
                }
        }
        //occur = sb->num_seq;
        //occur =  motif->end_in_sa - motif->start_in_sa;
        //fprintf(stdout,"%f occur \n",occur);


        //exit(0);
        //w = 1;
        /* set z_ij */
        cc = 0;
        for( i = 0; i < sb->num_seq;i++){
                cc += sb->sequences[i]->seq_len - w;
        }
        //c /= 4;


        occur = sqrtf(sb_temp->num_seq);

        lambda_1 = (occur) / cc;
        lambda_0 = 1.0 - lambda_1;

        lambda_0 = prob2scaledprob(lambda_0);
        lambda_1 = prob2scaledprob(lambda_1);
        cc = 0.0;

        for( i = 0; i < sb_temp->num_seq;i++){
                s = sb_temp->sequences[i];
                for(j = 0; j < s->seq_len-w;j++){
                        if(s->u[j]){
                                a = prob2scaledprob(1.0);
                                b = prob2scaledprob(1.0);

                                for(c = 0; c < w;c++){
                                        state = motif->state_sequence[c];
                                        a = a + fhmm->e[state][s->seq[j+c]];
                                        //fprintf(stdout,"%f %f %f %d \n", a,b,fhmm->e[state][s->seq[j+c]],s->seq[j+c]);
                                }
                                a = s->u[j] * (a + lambda_1);

                                for(c = 0; c < w;c++){
                                        b = b + fhmm->background[s->seq[j+c]];
                                }
                                a += (1.0 - s->u[j])* (b+ lambda_0);
                                //fprintf(stdout,"%f %f \n", a,b);
                                cc += a - b;
                        }
                }
        }
        //fprintf(stdout,"%f log", cc);
        // exit(0);
        /*RUN(find_overlapping_hits(sa, motif->start_in_sa,motif->end_in_sa,w));
        occur =  motif->end_in_sa - motif->start_in_sa;
        //fprintf(stdout,"%f occur \n",occur);
        occur = 0.0f;
        for(i = motif->start_in_sa; i < motif->end_in_sa;i++){
                if(sa->lcs[i]->lcp){
                        occur++;
                }
        }
        cc = 0;
        for( i = 0; i < sb->num_seq;i++){
                cc += sb->sequences[i]->seq_len - w;
        }
        //c /= 4;

        lambda_0 = (cc - occur) / cc;
        lambda_1 = (occur) / cc;


        lambda_0 = prob2scaledprob(lambda_0);
        lambda_1 = prob2scaledprob(lambda_1);
        cc = 0.0;
        //fprintf(stdout,"%f occur \n",occur);

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
        */
        motif->log_likelihood  = cc;

        //fprintf(stdout,"occ:%d\tllr:%f\n",motif->end_in_sa - motif->start_in_sa,cc);
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


int write_meme_output(char* filename,struct motif_list* m, struct fhmm*  fhmm)
{
        struct paraclu_cluster* p =  NULL;
        FILE* fptr = NULL;
        int i,j,c;
        float sum;
        ASSERT(filename != NULL, "No file name");
        ASSERT(m != NULL, "No motifs");
        if(fhmm->L == ALPHABET_DNA){
                RUNP(fptr = fopen(filename, "w"));
                fprintf(fptr,"MEME version 4\n");
                fprintf(fptr,"\n");
                fprintf(fptr,"ALPHABET= ACGT\n");
                fprintf(fptr,"\n");
                fprintf(fptr,"strands: +\n");
                fprintf(fptr,"\n");
                fprintf(fptr,"Background letter frequencies\n");
                for(i = 0; i < fhmm->L;i++){
                        fprintf(fptr,"%c %0.3f","ACGT"[i],fhmm->background[i]);
                        if(i != 3){
                                fprintf(fptr," ");
                        }
                }
                fprintf(fptr,"\n");
                fprintf(fptr,"\n");


                for(i = 0 ; i < m->num_items;i++){
                        p = m->plist[i];
                        if(p){
                                fprintf(fptr,"MOTIF M%03d\n", i+1);
                                fprintf(fptr,"letter-probability matrix: alength= %d w= %d nsites= %d E= 4.1e-009\n",fhmm->L, p->len, p->num_hits);
                                for(j = 0 ; j < p->len;j++){
                                        sum = 0.0f;
                                        for(c = 0; c < fhmm->L;c++){
                                                sum += p->count_matrix[j][c];
                                        }
                                        fprintf(fptr,"%0.6f", p->count_matrix[j][0]/sum);
                                        for(c = 1; c < fhmm->L;c++){
                                                fprintf(fptr," %0.6f", p->count_matrix[j][c]/sum);
                                        }
                                        fprintf(fptr,"\n");
                                }
                                fprintf(fptr,"\n");
                        }
                }

                fclose(fptr);

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
                                tmp[j][c] =  p->count_matrix[j][c];//   fhmm->e[p->state_sequence[j]][c];// p->matrix[j][c];//
                        }
                }


                hdf5_data->rank = 2;
                hdf5_data->dim[0] = p->len;
                hdf5_data->dim[1] = fhmm->L;
                hdf5_data->chunk_dim[0] = p->len;
                hdf5_data->chunk_dim[1] = fhmm->L;
                hdf5_data->native_type = H5T_NATIVE_FLOAT;
                snprintf(buffer, BUFFER_LEN, "M%03d;%4.2f;%d", i+1,p->log_likelihood, p->num_present_in_seq );
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

int randomize_freq_matrix(float** m,int len, int L)
{
        int i,j,c;
        float tmp;

        /* column wise  */
        for(i = 0; i < len-1; i++){
                j = i + (int) rk_interval((len-1) - i, &rndstate);
                for(c = 0; c < L;c++){
                        tmp = m[j][c];
                        m[j][c] = m[i][c];
                        m[i][c] = tmp;
                }

        }

        /* row wise; */

        for(i = 0; i < len; i++){
                for(c = 0; c < L-1;c++){
                        j = c + (int) rk_interval((L-1)-c, &rndstate);
                        tmp = m[i][j];
                        m[i][j] = m[i][c];
                        m[i][c] = tmp;
                }
        }
        /*float sum = 0.0f;
        float average = 0.0f;
        float k,l;
        float a[4];
        float b[4];
        for(c = 1; c < 100;c++){
                average= 0.0;
                for(i = 0; i < 100000;i++){
                        sum = 0.0f;
                        while(sum == 0.0f){
                                sum = 0.0f;
                                for(j =0; j < L;j++){
                                        a[j] = rk_gamma(&rndstate, 10.0 / (float) c, 1.0);
                                        sum += a[j];

                                }
                        }

                        for(j =0; j < L;j++){
                                a[j] /= sum;
                                //fprintf(stdout,"%f ",a[j]);

                        }
                        //fprintf(stdout,"\n");
                        ASSERT(sum != 0.0f, "No sum");
                        sum = 0.0f;
                        while(sum == 0.0f){
                                sum = 0.0f;
                                for(j =0; j < L;j++){
                                        b[j] = rk_gamma(&rndstate, 10.0 / (float) c, 1.0);
                                        sum += b[j];

                                }
                        }
                        ASSERT(sum != 0.0f, "No sum");

                        for(j =0; j < L;j++){
                                b[j] /= sum;
                        }
                        k = 0.0f;
                        l = 0.0f;
                        for(j =0; j < L;j++){

                                k += a[j] * log2f(a[j] / b[j]);
                                l += b[j] * log2f(b[j] / a[j]);

                        }
                        if(k+l > 0.0){
                                average+= (k+l) / 2.0f;
                        }

                }
                average /= 100000.0f;
                fprintf(stdout,"beta:%f\t%f %f\n",10.0 / (float) c,average,1 / 10.0 / (float) c );
        }
        exit(0);*/
        /*for(j = 0; j < L;j++){
                fprintf(stdout,"%c  ","ACGT"[j]);
                for(i = 0; i < len;i++){
                        fprintf(stdout," %0.3f",m[i][j]);
                }
                fprintf(stdout,"\n");
        }

        fprintf(stdout,"\n");*/
        return OK;
//ERROR:
//        exit(0);
}

int insert_motif(struct motif_list* m, struct paraclu_cluster* p,struct fhmm* fhmm,struct seq_buffer *sb)
{
        struct paraclu_cluster* a;


        int i;
        int new = 1;
        float score = 0.0f;
        float best_score = 1000000.0f;

        int hit = -1;
        int allowed_overlap;

        if(m->num_items){
                /*rand_freq_matrix = malloc_2d_float(rand_freq_matrix, p->len, fhmm->L, 0.0);
                for(i = 0; i < p->len;i++){
                        for(j = 0; j < fhmm->L;j++){
                                rand_freq_matrix[i][j] = p->freq_matrix[i][j];
                        }
                }
                MMALLOC(background_scores, sizeof(float) * num_random);

                for(j = 0; j < num_random;j++){
                        RUN(randomize_freq_matrix(rand_freq_matrix, p->len,fhmm->L));
                        a =  m->plist[rk_interval(m->num_items-1, &rndstate)];
                        allowed_overlap = MACRO_MIN(a->len, p->len);
                        RUN(motif_dyn_programming(fhmm, a->freq_matrix,rand_freq_matrix, a->len, p->len,allowed_overlap));
                        background_scores[j] = fhmm->F_matrix[a->len][p->len]/allowed_overlap ;
                }

                qsort(background_scores, num_random, sizeof(float),float_cmp);*/

                for(i = 0; i < m->num_items;i++){
                        a = m->plist[i];
                        allowed_overlap = MACRO_MIN(a->len, p->len);
                        RUN(motif_dyn_programming(fhmm, a->freq_matrix, p->freq_matrix, a->len, p->len,allowed_overlap));
                        score = fhmm->F_matrix[a->len][p->len]/allowed_overlap;
                        //RUN(align_motif(fhmm, a , p, &score));

                        /* build background distribution */



                        /*c = 0;
                        for(j = 0; j < num_random;j++){

                                if(background_scores[j] > score ){
                                        break;
                                }
                                c++;
                                if(c > 101){
                                        break;
                                }

                                }*/
                        //if(c < 10){
                                //for(j = 0; j < c+2;j++){
                                //        fprintf(stdout,"%d: %f <%f> %d\n",j, background_scores[j],score,c);
                                //}
                        //}
                        //fprintf(stdout,"there are %d scores better than %f (e.g. %f)\n",c,score,background_scores[c+1]);
                        if(score <= 1.5){//} && c < 100){
                                if(score < best_score){
                                        best_score = score;
                                        hit = i;
                                }
                        }

                        //fprintf(stdout,"%d %f   %d\n",i, score,m->num_items);
                }
                //fprintf(stdout,"Best:%f hit: %d\n",best_score,hit );
                if(hit != -1){

                        a = m->plist[hit];
                        allowed_overlap = MACRO_MIN(a->len, p->len);
                        RUN(motif_dyn_programming(fhmm , a->freq_matrix, p->freq_matrix, a->len, p->len,allowed_overlap));

                        RUN(traceback_and_merge_motifs(fhmm,a,p));
                        RUN(add_motif_count_matrix_based_on_hits(sb,a));
                        new = 0;
                        free_paraclu_cluster(p);
                        }
                //free_2d((void**) rand_freq_matrix);

        }
        //for(i = 0; i < m->num_items;i++){
        //RUN(compare_hit_positions(  m->plist[i]->hits , p->hits, m->plist[i]->num_hits, p->num_hits, MACRO_MAX( m->plist[i]->len, p->len),&jac));
                //if(jac >= 0.0){
                //RUN(align_motif(fhmm, m->plist[i], p, &score));
                //        if(p->len == -1){
                                //merge_motif(m->plist[i],p);
                //                 new = 0;
                //               free_paraclu_cluster(p);
                                //                break;
                                //        }
                        //}
                /* compare motifs - keep longer one... */
                /*RUN(compare_motif(fhmm, m->plist[i], p, &kl));
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
                        }*/

                //}
        //fprintf(stdout,"Insert? %d\n",new);
        if(new){                /* still */
                m->plist[m->num_items] = p;
                m->num_items++;
                if(m->num_items == m->alloc_items){
                        RUN(extend_motif_list(m));
                }
        }
        qsort(m->plist , m->num_items, sizeof(struct paraclu_cluster*), sort_paraclu_cluster_based_on_hits);

        return OK;
ERROR:
        return FAIL;
}

int compare_hit_positions(struct hit** list_a, struct hit** list_b, int len_a,int len_b, int w,float *jac)
{
        int intersection;
        int only_a;
        int only_b;
        int a,b;

        a = 0;
        b = 0;
        intersection =0;
        only_a = 0;
        only_b = 0;

        while(1){
                if(list_a[a]->seq == list_b[b]->seq){
                        if(list_a[a]->pos <= list_b[b]->pos){
                                if(list_a[a]->pos + w >= list_b[b]->pos){
                                        intersection++;
                                        a++;
                                        b++;
                                }else{
                                        only_a++;
                                        a++;
                                }
                        }else{
                                if(list_b[b]->pos + w >= list_a[a]->pos){
                                        intersection++;
                                        a++;
                                        b++;
                                }else{
                                        only_b++;
                                        b++;
                                }
                        }
                }else if(list_a[a]->seq < list_b[b]->seq){
                        only_a++;
                        a++;
                }else{
                        only_b++;
                        b++;
                }
                if(a == len_a){
                        break;
                }
                if(b == len_b){
                        break;
                }
        }

        if(b != len_b){
                only_b += len_b - b;
        }
        if(a != len_a){
                only_a += len_a - a;
        }
        //fprintf(stdout,"i:%d a:%d b:%d    entries: %d %d Jaccard :%f\n", intersection,only_a,only_b, len_a,len_b  ,  (float) intersection / (float)(intersection + only_a + only_b));
        *jac = MACRO_MAX((float) intersection /(float)(len_a),(float) intersection / (float)(len_b));// (float) intersection / (float)(intersection + only_a + only_b);
        return OK;
}

int merge_motif(struct paraclu_cluster* a,struct paraclu_cluster* b)
{
        int i;

        ASSERT(a != NULL, "Nno a motif");
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

int traceback_and_merge_motifs(struct fhmm* fhmm,struct paraclu_cluster* a,struct paraclu_cluster* b)
{
        int i,j,c;
        int new_len;
        int pos_a = 0;
        int pos_b = 0;

        uint8_t* path = NULL;
        float** t;
        int* new_state_sequence = NULL;
        //float** new_count_matrix = NULL;

        ASSERT(fhmm != NULL, "No fhmm model");
        ASSERT(a != NULL,"No A motif");
        ASSERT(b != NULL,"No B motif");

        t = fhmm->B_matrix;

        MMALLOC(path, sizeof(uint8_t) * (a->len + b->len +1));
/*fprintf(stdout,"Motif A\n");
  for(c = 0; c < fhmm->L;c++){
  fprintf(stdout,"%c","ACGT"[c]);
  for(i = 0; i < a->len;i++){
  fprintf(stdout," %0.2f", fhmm->e[a->state_sequence[i]][c]);
  }
  fprintf(stdout,"\n");
  }
  fprintf(stdout,"Motif B\n");

  for(c = 0; c < fhmm->L;c++){
  fprintf(stdout,"%c","ACGT"[c]);
  for(j = 0; j < b->len;j++){
  fprintf(stdout," %0.2f", fhmm->e[b->state_sequence[j]][c]);
  }

  fprintf(stdout,"\n");
  }
  fprintf(stdout,"Dyn matrix\n");
  for(i = 0; i < a->len+1;i++){
  for(j = 0; j < b->len+1;j++){

  fprintf(stdout," %2.2f",m[i][j]);
  }
  fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n");

  fprintf(stdout,"Traceback matrix\n");
  for(i = 0; i < a->len+1;i++){
  for(j = 0; j < b->len+1;j++){

  fprintf(stdout," %d",(int)t[i][j]);
  }
  fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n");
*/
        /* Traceback to get path */

        i = a->len;
        j = b->len;
        new_len = 0;
        while( i + j != 0){
                switch ((int)t[i][j]) {
                case 0: {
                        path[new_len] = 0;
                        i--;
                        j--;
                        break;
                }
                case 1: {
                        path[new_len] = 1;
                        i--;
                        break;
                }
                case 2: {
                        path[new_len] = 2;
                        j--;
                        break;
                }
                default:
                        break;
                }
                new_len++;
        }
        //if(0){
        //new_count_matrix = malloc_2d_float(new_count_matrix, new_len, fhmm->L, 0.0f);
        MMALLOC(new_state_sequence, sizeof(int) * (new_len));
        j = 0;
        for(i = new_len-1; i >=0;i--){
                //fprintf(stdout, "%d\t",path[i]);
                switch (path[i]) {
                case 0: {

                        /*for(c = 0;c < fhmm->L;c++){
                          fprintf(stdout," %0.2f", fhmm->e[a->state_sequence[pos_a]][c]);
                          }
                          for(c = 0;c < fhmm->L;c++){
                          fprintf(stdout," %0.2f", fhmm->e[b->state_sequence[pos_b]][c]);
                          }*/

                        pick_state(fhmm->e[a->state_sequence[pos_a]],fhmm->e[b->state_sequence[pos_b]], fhmm->L, &c);
                        if(c == 1){
                                new_state_sequence[j] = a->state_sequence[pos_a];
                        }else{
                                new_state_sequence[j] = b->state_sequence[pos_b];
                        }
                        for(c = 0;c < fhmm->L;c++){

                                //                      new_count_matrix[j][c] = a->count_matrix[pos_a][c]+b->count_matrix[pos_b][c];
                        }
                        pos_a++;
                        pos_b++;
                        break;
                }
                case 1: {
                        /*        for(c = 0;c < fhmm->L;c++){
                                  fprintf(stdout," %0.2f", fhmm->e[a->state_sequence[pos_a]][c]);
                                  }
                                  for(c = 0;c < fhmm->L;c++){
                                  fprintf(stdout," %0.2f", 0.0);
                                  }*/
                        new_state_sequence[j] = a->state_sequence[pos_a];
                        for(c = 0;c < fhmm->L;c++){

                                //new_count_matrix[j][c] = a->count_matrix[pos_a][c];
                        }
                        pos_a++;
                        break;
                }
                case 2: {
                        /*for(c = 0;c < fhmm->L;c++){
                          fprintf(stdout," %0.2f", 0.0);
                          }
                          for(c = 0;c < fhmm->L;c++){
                          fprintf(stdout," %0.2f", fhmm->e[b->state_sequence[pos_b]][c]);
                          }*/
                        new_state_sequence[j] = b->state_sequence[pos_b];
                        for(c = 0;c < fhmm->L;c++){
                                //new_count_matrix[j][c] = b->count_matrix[pos_b][c];
                        }
                        pos_b++;
                        break;
                }
                default:
                        break;
                }
                //fprintf(stdout, "\n");
                j++;
        }

        //}

        pos_a = 0;
        j = new_len-1;
        while(path[j] == 1){
                pos_a++; /* All hits in b have to pushed back pos_a positions */
                j--;
        }

        if(pos_a){
                for(j = 0; j < b->num_hits;j++){
                        b->hits[j]->pos -= pos_a;
                }
        }

        pos_b= 0;
        j = new_len-1;
        while(path[j] == 2){
                pos_b++; /* All hits in b have to pushed back pos_a positions */
                j--;
        }

        if(pos_b){
                for(j = 0; j < a->num_hits;j++){
                        a->hits[j]->pos -= pos_b;
                }
        }
        //fprintf(stdout,"pa:%d pb:%d\n",pos_a,pos_b);
        //fprintf(stdout,"\n");
        /*for(j = 0; j < new_len;j++){
          fprintf(stdout," %d",new_state_sequence[j]);//fhmm->e[b->state_sequence[j]][i]);
          }
          fprintf(stdout,"\n");
          for(i = 0; i < fhmm->L;i++){
          fprintf(stdout,"%c","ACGT"[i]);
          for(j = 0; j < new_len;j++){
          fprintf(stdout," %0.2f", new_count_matrix[j][i]);//fhmm->e[b->state_sequence[j]][i]);
          }

          fprintf(stdout,"\n");
          }*/
        //free_2d((void**) a->count_matrix);

            //    a->count_matrix = NULL;
            //a->count_matrix = new_count_matrix;
        a->len = new_len;
        MFREE(a->state_sequence);
        a->state_sequence = new_state_sequence;
        MREALLOC(a->hits, sizeof(struct hit*) *(a->num_hits + b->num_hits));
        for(j = 0; j < b->num_hits;j++){
                a->hits[a->num_hits] = b->hits[j];
                a->num_hits++;
                b->hits[j] = NULL;
                //        fprintf(stdout,"%d %d\n", a->hits[j]->seq,a->hits[j]->pos);
        }
        MFREE(b->hits);
        //fprintf(stdout,"\n");
        b->hits = NULL;
        qsort(a->hits, a->num_hits, sizeof(struct hit*), sort_hit_positions);
        for(j = 1; j < a->num_hits;j++){
                if( a->hits[j]->seq ==  a->hits[j-1]->seq){
                        if(a->hits[j]->pos ==  a->hits[j-1]->pos){
                                a->hits[j]->seq = 1000000;
                                a->hits[j]->pos = 1000000;
                        }
                }
        }
        //qsort(a->hits, a->num_hits, sizeof(struct hit*), sort_hit_positions);
        //RUN(remove_overlapping_hits(a));
        qsort(a->hits, a->num_hits, sizeof(struct hit*), sort_hit_positions);
        c = a->num_hits;
        for(j = 0; j < a->num_hits;j++){
                if(a->hits[j]->seq == 1000000){
                        c = j;
                        break;
                }
                //fprintf(stdout,"%d %d\n", a->hits[j]->seq,a->hits[j]->pos);
        }
        //fprintf(stdout,"%d\n",c);

        for(j = c; j < a->num_hits;j++){
                MFREE(a->hits[j]);
                a->hits[j] = NULL;
        }
        a->num_hits = c;

        MREALLOC(a->hits, sizeof(struct hit*) *(a->num_hits));
        b->len = -1;
        a->log_likelihood = MACRO_MAX(a->log_likelihood, b->log_likelihood);

        merge_motif(a,b);

        MFREE(path);
        return OK;
ERROR:
        MFREE(path);
        return FAIL;

}
int align_motif(struct fhmm* fhmm, struct paraclu_cluster* a,struct paraclu_cluster* b, float* score)
{

        int allowed_overlap;
        //int* permuted_states = NULL;
        //unsigned int seed = 42;

        ASSERT(fhmm != NULL, "No fhmm model");
        ASSERT(a != NULL,"No A motif");
        ASSERT(b != NULL,"No B motif");

        //MMALLOC(permuted_states, sizeof(int) * (a->len));


        allowed_overlap = MACRO_MIN(a->len, b->len);
        /* Idea: dynamic programming to align a and b. Allow for
         * leading and trailing gaps; no gaps in the middle and a
         * minimum alignment length *(6 probably) */

        /* Don't need to allocate a matrix - done in the main function... */
        RUN(motif_dyn_programming(fhmm, a->freq_matrix, b->freq_matrix, a->len, b->len,allowed_overlap));

        *score = fhmm->F_matrix[a->len][b->len];
        //fprintf(stdout,"SCORE:%f\n",score);
          /*if(score){
                for(i = 0; i < a->len;i++){
                        permuted_states[i] = a->state_sequence[i];
                }

                c = 0;
                for(i = 0; i < 10000;i++){
                        RUN(shuffle_arr_r(permuted_states,a->len,&seed));
                        RUN(motif_dyn_programming(fhmm, permuted_states, b->state_sequence, a->len, b->len,6));
                        if(m[a->len][b->len] <= score){
                                c++;
                                if(c == 100){
                                        break;
                                }
                        }
                }
        }else{
                c =0;
        }*/
        /*fprintf(stdout,"%d %f \n",c,score);
        for(i = 0; i < a->len;i++){
                fprintf(stdout," %d", a->state_sequence[i]);
        }
        fprintf(stdout,"\n");
        for(i = 0; i < b->len;i++){
                fprintf(stdout," %d", b->state_sequence[i]);
        }
        fprintf(stdout,"\n");*/
        //if(*score <= 3.0){
        //        RUN(motif_dyn_programming(fhmm, a->state_sequence, b->state_sequence, a->len, b->len,allowed_overlap));

        //       RUN(traceback_and_merge_motifs(fhmm,a,b));
        //}

        //MFREE(permuted_states);

        return OK;
ERROR:
        //MFREE(permuted_states);
        return FAIL;
}


int pick_state(float* a, float* b, int len, int* pick)
{
        float ma,mb;
        int i;

        *pick = 0;
        ma = -1.0f;
        mb = -1.0f;
        for(i = 0; i < len;i++){
                if(a[i] > ma){
                        ma = a[i];
                }
                if(b[i] > mb){
                        mb = b[i];
                }


        }
        if(ma < mb){
                *pick = 1;
        }else{
                *pick = 2;
        }


        return OK;
}


int motif_dyn_programming(struct fhmm* fhmm, float** e_a, float** e_b, int len_a, int len_b,int min_aln_len)
{
        //float** e = fhmm->e;
        float** m;
        float** t;
        int i,j,c;
        float x,y;
        m = fhmm->F_matrix;
        t = fhmm->B_matrix;
        /* Step 1: fill matrix */

        m[0][0] = 0.0f;


        for(i = 1; i < len_a+1;i++){
                m[i][0] = 0.0f;
                t[i][0] = 1;
        }
        for(j = 1; j < len_b+1;j++){
                m[0][j] = 0.0f;
                t[0][j] = 2;
        }


        for(i = 1; i < len_a+1;i++){
                for(j = 1; j < len_b+1;j++){
                        x = 0.0f;
                        y = 0.0f;
                        for(c = 0; c < fhmm->L;c++){
                                x += e_a[i-1][c] * log2f(e_a[i-1][c] / e_b[j-1][c]);
                                y += e_b[j-1][c] * log2f(e_b[j-1][c] / e_a[i-1][c]);
                        }
                        m[i][j] = ((x + y)) / 2.0f + m[i-1][j-1];
                        t[i][j] = 0;
                }
        }

        /* Fill transition for last column / row */

        /* Disallow all alignments with fewer than 6 aligned positions*/

        i =  len_a;
        for(j = min_aln_len+1; j < len_b+1;j++){
                if(m[i][j-1] < m[i][j]){
                        m[i][j] = m[i][j-1];
                        t[i][j] = 2;
                }
        }
        j = len_b;
        for(i = min_aln_len+1; i < len_a+1;i++){
                if(m[i-1][j] < m[i][j]){
                        m[i][j] = m[i-1][j];
                        t[i][j] = 1;
                }
        }
        i = len_a;
        j = len_b;
        if(m[i][j-1] < m[i][j]){
                m[i][j] = m[i][j-1];
                t[i][j] = 2;
        }

        if(m[i-1][j] < m[i][j]){
                m[i][j] = m[i-1][j];
                t[i][j] = 1;
        }
        return OK;
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
        p->hits = NULL;
        p->num_hits = 0;
        p->count_matrix = NULL;
        p->freq_matrix = NULL;
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
        int i;
        if(p){
                if(p->hits){
                        for(i = 0; i < p->num_hits;i++){
                                MFREE(p->hits[i]);
                        }
                        MFREE(p->hits);

                }
                if(p->count_matrix){
                        free_2d((void**)p->count_matrix);
                }
                if(p->freq_matrix){
                        free_2d((void**)p->freq_matrix);
                }
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


int shuffle_arr_r(int* arr,int n,unsigned int* seed)
{
        int i,j;
        int tmp;
        for (i = 0; i < n - 1; i++) {
                j = i +  (int) (rand_r(seed) % (int) (n-i));
                tmp = arr[j];
                arr[j] = arr[i];
                arr[i] = tmp;
        }
        return OK;
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

int sort_paraclu_cluster_based_on_hits(const void *a, const void *b)
{
    struct paraclu_cluster* const *one = a;
    struct paraclu_cluster* const *two = b;

    if((*one)->num_hits <  (*two)->num_hits ){
            return 1;
    }else{
            return -1;
    }
}

int sort_hit_positions(const void *a, const void *b)
{
        struct hit* const *one = a;
        struct hit* const *two = b;
        if((*one)->seq > (*two)->seq){
                return 1;
        }else if( (*one)->seq == (*two)->seq){
                if((*one)->pos >  (*two)->pos){
                        return 1;
                }else{
                        return -1;
                }
        } else{
                return -1;
        }
}


int float_cmp(const void *a, const void *b)
{
    const float *ia = (const float *)a;
    const float *ib = (const float *)b;
    if(*ia > *ib ){
            return 1;
    }else{
            return -1;
    }
}