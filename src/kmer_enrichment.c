#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>

#include "tldevel.h"
#include "hash_table.h"
#include "ihmm_seq.h"

#include "finite_hmm.h"
#include "emit_random.h"
#include "run_score.h"

#include "model.h"

#include "kalign.h"

#include "motif_refinement.h"

int esl_stats_GTest(int ca, int na, int cb, int nb, double *ret_G, double *ret_P);
int esl_stats_IncompleteGamma(double a, double x, double *ret_pax, double *ret_qax);
int esl_stats_LogGamma(double x, double *ret_answer);

struct parameters{
        char* in_model;
        unsigned long seed;
        int kmer_len;
        int N;
};



struct motif{
        double** freq_matrix;
        double** count_matrix;
        double log_likelihood;
        double dust_score;
        int L;
        int len;
        int count;
};

struct alignment{
        struct motif* motif;
        char** sequences;
        int* seq_counts;
        int num_seq;
        int alloc_num_seq;

};

struct motif_list{
        struct motif** mlist;
        int num_items;
        int alloc_items;


};

rk_state rndstate;


struct alignment* init_alignment(int e_num_seq);
int resize_alignment(struct alignment* a);
void free_alignment(struct alignment* a);
int create_motif_based_on_alignment(struct alignment* a,double* background, int L);



int dust_sequence( uint8_t* seq, int len,double* score);
struct motif_list* init_motif_list(int num);
int extend_motif_list(struct motif_list* ml);
void free_motif_list(struct motif_list* ml);
int insert_motif(struct motif_list* ml, struct motif* m);


int hierarchal_merge_motif(struct motif_list* ml,struct fhmm* fhmm,int best_i,int best_j);
int set_overlap(int len_a, int len_b);
int motif_dyn_programming(struct fhmm* fhmm, double** e_a, double** e_b, int len_a, int len_b,int min_aln_len);
int randomize_freq_matrix(double** m,int len, int L);
int traceback_and_merge_motifs(struct fhmm* fhmm,struct motif* a,struct motif* b,struct motif* new_motif);

int write_meme_output(char* filename,struct motif_list* m, struct fhmm*  fhmm);

static int run_kmer_counting(struct parameters* param);
static int per_state_kmer_counting(struct parameters* param);

struct motif* create_motif_based_on_kmer_count(char* seq, int len,int count, int L, double* background, int total_counts);
void free_motif(struct motif* m);

static int double_cmp(const void *a, const void *b);
static int sort_motif_list_based_on_likelihood(const void *a, const void *b);
static int compare_kmer_based_on_count(const void *a, const void *b);

static int print_help(char **argv);
static int free_parameters(struct parameters* param);

/* Initialize hash table functions  */
HT_GLOBAL_INIT(KMERHASH, char*);

int fill_background_hash(HT_TYPE(KMERHASH)* ht,int num_samples, int len,int L, double* background);

static int insert_kmers_from_sampled_sequences(HT_TYPE(KMERHASH)* ht, struct seq_buffer* sb, int num_sampled, int kmer_len);

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        print_program_header(argv, "Build HDPHMM model(s).");

        MMALLOC(param, sizeof(struct parameters));
        param->in_model = NULL;
        param->kmer_len = 10;
        param->N = 1000000;

        param->seed = 0;

        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"len",required_argument,0,'l'},
                        {"seed",required_argument,0,'s'},
                        {"num",required_argument,0,'n'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"m:tl:n:s:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 'm':
                        param->in_model = optarg;
                        break;
                case 'l':
                        param->kmer_len = atoi(optarg);
                        break;
                case 's':
                        param->seed = atoi(optarg);
                        break;
                case 'n':
                        param->N = atoi(optarg);
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

        if(!param->in_model){
                RUN(print_help(argv));
                ERROR_MSG("No model file! use -m  <blah.h5>");
        }else{
                if(!my_file_exists(param->in_model)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->in_model);
                }
        }

        if(!param->kmer_len){
                 RUN(print_help(argv));
                 ERROR_MSG("kmer len is 0! Use -l [4-16]");
        }


        if(param->kmer_len >= 50){
                 RUN(print_help(argv));
                 ERROR_MSG("kmer len is very long! Use -l [4-16]");
        }
        if(param->seed){
                rk_seed(param->seed, &rndstate);
        }else{
                rk_randomseed(&rndstate);
        }
        //RUN(score_sequences_for_command_line_reporting(param));

        //RUN(run_kmer_counting(param));

        per_state_kmer_counting(param);
        exit(0);
        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}



int per_state_kmer_counting(struct parameters* param)
{
        char buffer[BUFFER_LEN];
        HT_TYPE(KMERHASH)* ht_back = NULL;
        HT_TYPE(KMERHASH)* ht_state = NULL;
        struct model_bag* model_bag = NULL;
        struct alignment* a = NULL;
        struct seq_buffer* sb = NULL;
        struct seq_buffer* sb_org = NULL;
        struct ihmm_sequence* s = NULL;
        struct fhmm* fhmm = NULL;

        struct motif_list* ml = NULL;

        char* kmer = NULL;
        int best;
        int i,j,g,l;
        int num_samples = 100000;
        int kmer_len_start = 10;
        int kmer_len_end = 10;
        int sampled = 0;
        hash_table_node_KMERHASH_t* n = NULL;



        LOG_MSG("Reading models.");
        RUNP(model_bag = read_model_bag_hdf5(param->in_model));
        LOG_MSG("Done.");

        //ht_back = HT_INIT(KMERHASH,1000000);

        RUNP(fhmm = read_best_fmodel(param->in_model, &best));
        RUN(alloc_dyn_matrices(fhmm));
        RUN(realloc_dyn_matrices(fhmm, 128));
        RUN(convert_fhmm_log_to_prob_for_sampling(fhmm));
        RUNP(ml = init_motif_list(fhmm->K));

        RUNP(sb_org = get_sequences_from_hdf5_model(param->in_model, IHMM_SEQ_READ_ONLY_SEQ));
        for(i = 2; i < fhmm->K;i++){
                //LOG_MSG("Working on state: %d",i);

                ht_state = HT_INIT(KMERHASH,1000000);
                ht_back = HT_INIT(KMERHASH,1000000);
                for(l = kmer_len_start; l <= kmer_len_end;l++){
                        sampled = 0;
                        sb = emit_kmers_from_state(fhmm, i, num_samples, l,rk_ulong(&rndstate));
                        for(j = 0; j < num_samples;j++){
                                s = sb->sequences[j];
                                kmer= NULL;
                                RUNP(kmer = galloc(kmer, sizeof(char) * (l+1)));
                                if(l == s->seq_len){
                                        sampled++;
                                        for(g = 0; g < l;g++){
                                                kmer[g] = "ACGT"[s->seq[g]];

                                        }
                                        kmer[l] = 0;
                                        //fprintf(stdout ,"%s\n",  kmer);
                                        HT_INSERT(KMERHASH, ht_state, kmer, NULL);
                                }

                        }
                        free_ihmm_sequences(sb);
                        RUN(fill_background_hash(ht_back,sampled, l,fhmm->L,  fhmm->background));
                }
                //LOG_MSG("hash size: %d", ht_state->total_count );



                //LOG_MSG("hash size: %d", ht_back->total_count );

                HT_FLATTEN(KMERHASH,ht_state);
                qsort(ht_state->flat,ht_state->num_items, sizeof(hash_table_node_KMERHASH_t*),compare_kmer_based_on_count);
                double sum_g = 0.0;

                RUNP(a = init_alignment(512));

                for(j = 0; j < ht_state->num_items;j++){

                        //fprintf(stdout,"%s\t%d\t",ht_state->flat[j]->key,ht_state->flat[j]->count);
                        n =  HT_SEARCH(KMERHASH,ht_back,ht_state->flat[j]->key);
                        double ret_G,ret_p;
                        l = 0;
                        if(n){
                                l = n->count;
                                //fprintf(stdout,"COUNT:%d\n", l);

                        }

                        RUN(esl_stats_GTest(ht_state->flat[j]->count,
                                        ht_state->total_count,
                                        l,
                                        ht_back->total_count,
                                        &ret_G,
                                            &ret_p));
                        if(ret_p * (double)ht_state->total_count *(double) (fhmm->K -2) < 1e-6){
                                l = strlen(ht_state->flat[j]->key)+2;
                                MMALLOC( a->sequences[a->num_seq], sizeof(char) * l);

                                strncpy(a->sequences[a->num_seq], ht_state->flat[j]->key, l-2);
                                a->sequences[a->num_seq][l-2] =0;
                                a->seq_counts[a->num_seq] = ht_state->flat[j]->count;
                                a->num_seq++;


                                if(a->num_seq == a->alloc_num_seq){
                                        RUN(resize_alignment(a));


                                }

                                fprintf(stdout,"%s\t%d\t%d\t",ht_state->flat[j]->key,ht_state->flat[j]->count,l);
                                fprintf(stdout,"%f\t%e\n",ret_G, ret_p * (double)ht_state->total_count * (double) (fhmm->K -2)) ;
                                sum_g += ret_G;
                                if(a->num_seq == 100){
                                        break;
                                }
                        }


                }
                if(a->num_seq){

                        a->sequences = kalign_align(a->sequences,a->num_seq);


                        RUN(create_motif_based_on_alignment(a,fhmm->background, fhmm->L));
                        //em_algorithm(a->motif->count_matrix, a->motif->len, fhmm->L, sb_org);
                        //exit(0);
                        //for(j =0; j <  a->num_seq;j++){
                        //        fprintf(stdout,"%d\t%d\t%s\n",j, a->seq_counts[j], a->sequences[j]);
                        // }
                        fprintf(stdout,"%d\t%f\tSumG \n",i,sum_g);
                        a->motif->count = (int)((double) a->motif->count / (double) num_samples * 1000.0);
                        RUN(insert_motif(ml,a->motif));
                        a->motif = NULL;

                }
                free_alignment(a);
                a = NULL;
                HT_FREE(KMERHASH,ht_state);
                HT_FREE(KMERHASH,ht_back);
        }

         for(i = fhmm->L-1; i >= 1;i--){
                fhmm->background[i] -= fhmm->background[i-1];
        }



        qsort(ml->mlist , ml->num_items, sizeof(struct motif*), sort_motif_list_based_on_likelihood);

        hierarchal_merge_motif(ml,fhmm,-1,-1);
        snprintf(buffer, BUFFER_LEN, "model_%d.meme",1 );




        RUN(write_meme_output(buffer,ml,fhmm));
        free_motif_list(ml);

        free_fhmm(fhmm);
        free_ihmm_sequences(sb_org);

        return OK;
ERROR:
        return FAIL;
}



int fill_background_hash(HT_TYPE(KMERHASH)* ht,int num_samples, int len,int L, double* background)
{

        int i,j,c;
        char* kmer = NULL;
        double r;


        for(i = 0; i < num_samples;i++){
                kmer = NULL;
                RUNP(kmer = galloc(kmer, sizeof(char) * (len+1)));
                for(j = 0; j < len;j++){
                        r = rk_double(&rndstate);
                        //fprintf(stdout,"%f\n",r);
                        for(c = 0 ; c < L;c++){
                                if(r <= background[c]){
                                        kmer[j] = "ACGT"[c];
                                        break;
                                }
                        }
                }
                kmer[len] = 0;
                RUN(HT_INSERT(KMERHASH,ht,kmer,NULL));
        }
        return OK;
ERROR:
        return FAIL;

}
int run_kmer_counting(struct parameters* param)
{
        char buffer[BUFFER_LEN];
        struct model_bag* model_bag = NULL;

        struct seq_buffer* sb = NULL;
        struct fhmm* fhmm = NULL;
        struct motif* m = NULL;

        struct motif_list* ml = NULL;
        int best;
        int i;


        //char buffer[BUFFER_LEN];



        LOG_MSG("Reading models.");
        RUNP(model_bag = read_model_bag_hdf5(param->in_model));
        LOG_MSG("Done.");

        HT_TYPE(KMERHASH)* ht = NULL;
        //RUNP(all_scores = galloc(all_scores, limit,model_bag->num_models, 0.0));
        RUNP(fhmm = read_best_fmodel(param->in_model, &best));
        //for(c = 0; c < model_bag->num_models;c++){
        //        fhmm = model_bag->finite_models[c];
        RUN(alloc_dyn_matrices(fhmm));
        RUN(realloc_dyn_matrices(fhmm, 128));
        RUN(convert_fhmm_log_to_prob_for_sampling(fhmm));
        /* Horrible hack - we ASSUME that fhmm->r_score is 0.0  this will force emit to convert log probs first time it get's called  */

        RUNP(ml = init_motif_list(1024));
        ht = HT_INIT(KMERHASH,10000);
        while(ht->total_count  != param->N){

                RUNP(sb = emit_sequences_from_fhmm_model(fhmm,1000, rk_ulong(&rndstate)));



                RUN(insert_kmers_from_sampled_sequences(ht, sb, param->N, param->kmer_len));

                free_ihmm_sequences(sb);
                sb = NULL;
                LOG_MSG("%d items sampled.",ht->total_count);
        }
        HT_FLATTEN(KMERHASH,ht);

        qsort(ht->flat,ht->num_items, sizeof(hash_table_node_KMERHASH_t*),compare_kmer_based_on_count);

        /* convert background into a normal probability vector */
        for(i = fhmm->L-1; i >= 1;i--){
                fhmm->background[i] -= fhmm->background[i-1];
        }


        for(i = 0; i < ht->num_items;i++){
                m = create_motif_based_on_kmer_count( ht->flat[i]->key,param->kmer_len ,ht->flat[i]->count ,  fhmm->L, fhmm->background, ht->total_count);

                if(m->log_likelihood >= 10 && i <= 250){
                        RUN(insert_motif(ml,m));
                }else{
                        free_motif(m);
                }


        }
        HT_FREE(KMERHASH,ht);
        ht = NULL;

        qsort(ml->mlist , ml->num_items, sizeof(struct motif*), sort_motif_list_based_on_likelihood);
        for(i = 0; i < ml->num_items;i++){
                m = ml->mlist[i];
                fprintf(stdout,"Count:%d\t\tLL:%f\t Dust:%f\n",m->count,m->log_likelihood,m->dust_score );
        }
        hierarchal_merge_motif(ml,fhmm,-1,-1);
        snprintf(buffer, BUFFER_LEN, "model_%d.meme",1 );
        RUN(write_meme_output(buffer,ml,fhmm));
        free_motif_list(ml);
        free_fhmm(fhmm);
        //exit(0);


                //}
        //LOG_MSG("Got past seq free");
                // free_model_bag(model_bag);
        //LOG_MSG("Got past model bag free");


        return OK;
ERROR:
        free_model_bag(model_bag);
        return FAIL;
}

struct motif* create_motif_based_on_kmer_count(char* seq, int len,int count, int L, double* background, int total_counts)
{
        struct motif* m = NULL;
        int i,j,c;
        double sum;
        double lambda_0;
        double lambda_1;
        double cc;
        double a,b;
        uint8_t dustseq[128];
        double dust_score;
        MMALLOC(m, sizeof(struct motif));
        m->count_matrix = NULL;
        m->freq_matrix = NULL;
        m->len = len;
        m->log_likelihood = 0.0;
        m->count = count;
        m->L = L;
        RUNP(m->count_matrix = galloc(m->count_matrix,len,L,0.0));
        RUNP(m->freq_matrix = galloc(m->freq_matrix,len,L,0.0));


        lambda_1 = (double) m->count / (double) total_counts;
        lambda_0 = 1.0 - lambda_1;
        lambda_0 = prob2scaledprob(lambda_0);
        lambda_1 = prob2scaledprob(lambda_1);

        cc = 0.0;
        a = prob2scaledprob(1.0);
        b = prob2scaledprob(1.0);
        for(i = 0; i < m->len;i++){
                for(j = 0; j < m->L;j++){
                        m->count_matrix[i][j] = 0.0;//background[j];
                }
                c = seq[i] -65;
                dustseq[i] = c;

                m->count_matrix[i][c] += m->count;

                sum = 0.0;
                for(j = 0; j < m->L;j++){
                        sum += m->count_matrix[i][j] + background[j];
                }
                for(j = 0; j < m->L;j++){
                        m->freq_matrix[i][j] = (m->count_matrix[i][j]+background[j]) / sum;
                        //        fprintf(stdout,"%f ", m->freq_matrix[i][j]);
                }
                //fprintf(stdout,"\n");
                a = a + prob2scaledprob( m->freq_matrix[i][c]);
                //fprintf(stdout,"%f %f %f %d \n", a,b,fhmm->e[state][s->seq[j+c]],s->seq[j+c]);

                b = b + prob2scaledprob(background[c]);
        }

        dust_sequence( dustseq, m->len, &dust_score );
        a = 1.0 * (a + lambda_1);

        //a += (1.0 - s->u[j])* (b+ lambda_0); - evaluates to zero
//fprintf(stdout,"%f %f \n", a,b);
        //fprintf(stdout,"Lambda:%f\n", lambda_1);
        cc = (m->count)* (a - b);
        m->log_likelihood = cc;
        m->dust_score = dust_score;
        //fprintf(stdout,"Count:%d\tLambda:%f\tLL:%f  DUST:%f %s \n",m->count,lambda_1,cc,m->dust_score, seq );
        return m;
ERROR:
        return NULL;

}

void free_motif(struct motif* m)
{
        if(m){
                if(m->count_matrix){
                        gfree(m->count_matrix);
                }
                if(m->freq_matrix){
                        gfree(m->freq_matrix);
                }
                MFREE(m);
        }

}

int compare_kmer_based_on_count(const void *a, const void *b)
{
        hash_table_node_KMERHASH_t* const *one = a;



        hash_table_node_KMERHASH_t* const *two = b;
        if((*one)->count  < (*two)->count){
                return 1;
        }
        if((*one)->count  ==  (*two)->count){
                return 0;
        }
        return -1;

/* integer comparison: returns negative if b > a
   and positive if a > b */
}



int insert_kmers_from_sampled_sequences(HT_TYPE(KMERHASH)* ht, struct seq_buffer* sb, int num_sampled, int kmer_len)
{
        int i,j,g;
        struct ihmm_sequence* s = NULL;
        char* kmer = NULL;
        ASSERT(sb != NULL, "No sequences");
        for(i = 0 ; i < sb->num_seq;i++){

                s = sb->sequences[i];
                /* if(s->seq_len > kmer_len ){

                   RUNP(kmer = galloc(kmer, sizeof(char) * (kmer_len+1)));
                   strncpy(kmer,(char*) s->seq+0, kmer_len);
                   kmer[kmer_len] = 0;
                   RUN(HT_INSERT(KMERHASH,ht,kmer,NULL));
                   kmer= NULL;
                   if(ht->total_count == num_sampled){
                   return OK;
                   }


                   if(strncmp((char*) s->seq +1, (char*) s->seq +1-1 ,kmer_len)){
                   RUNP(kmer = galloc(kmer, sizeof(char) * (kmer_len+1)));
                   strncpy(kmer,(char*) s->seq+1, kmer_len);
                   kmer[kmer_len] = 0;
                   RUN(HT_INSERT(KMERHASH,ht,kmer,NULL));
                   kmer= NULL;
                   if(ht->total_count == num_sampled){
                   return OK;
                   }
                   }*/
                //c = rk_interval(kmer_len-1, &rndstate);

                for(j = 0;j < s->seq_len-kmer_len;j+=1){



                        //if(strncmp((char*) s->seq +j, (char*) s->seq +j-2 ,kmer_len)){
                        //        if(strncmp((char*) s->seq +j, (char*) s->seq +j-1 ,kmer_len)){
                        RUNP(kmer = galloc(kmer, sizeof(char) * (kmer_len+1)));
                        //strncpy(kmer,(char*) s->seq+j, kmer_len);
                        for(g = 0; g < kmer_len;g++){
                                kmer[g] = s->seq[j+g]+65;

                        }
                        kmer[kmer_len] = 0; /* make sure string is NULL terminated */
                        RUN(HT_INSERT(KMERHASH,ht,kmer,NULL));
                        kmer= NULL;
                        if(ht->total_count == num_sampled){
                                return OK;
                        }
                        //               }
                        //}



                }
                //}
        }


        return OK;
ERROR:
        return FAIL;

}



struct motif_list* init_motif_list(int num)
{
        struct motif_list* ml = NULL;
        int i;
        ASSERT(num > 1, "List needs at least two entry");
        MMALLOC(ml, sizeof(struct motif_list));
        ml->alloc_items = num;
        ml->num_items = 0;
        ml->mlist = NULL;

        MMALLOC(ml->mlist, sizeof(struct motif*) * ml->alloc_items);
        for(i = 0; i < ml->alloc_items;i++){
                ml->mlist[i] = NULL;
        }

        return ml;
ERROR:
        return NULL;
}

int extend_motif_list(struct motif_list* ml )
{
        int i,old_len;

        ASSERT(ml != NULL, "No list");
        old_len = ml->alloc_items;
        ml->alloc_items = ml->alloc_items << 1;
        MREALLOC(ml->mlist, sizeof(struct paraclu_cluster*)* ml->alloc_items);
        for(i = old_len; i < ml->alloc_items;i++){
                ml->mlist[i] = NULL;
        }
        return OK;
ERROR:

        return FAIL;
}

void free_motif_list(struct motif_list* ml)
{
        int i;
        if(ml){
                for(i = 0; i < ml->num_items;i++){
                        if(ml->mlist[i]){
                                free_motif(ml->mlist[i]);
                        }
                }
                MFREE(ml->mlist);
                MFREE(ml);
        }

}

int insert_motif(struct motif_list* ml, struct motif* m)
{

        ml->mlist[ml->num_items] = m;
        ml->num_items++;
        if(ml->num_items == ml->alloc_items){
                RUN(extend_motif_list(ml));
        }

        return OK;
ERROR:
        return FAIL;
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
        const char usage[] = " -m <model.h5> -[l,b] ";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--model","Input model." ,"[]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--len","kmer length." ,"[10]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--num","Number of Samples." ,"[1M]"  );
        return OK;
}

int sort_motif_list_based_on_likelihood(const void *a, const void *b)
{
        struct motif* const *one = a;
        struct motif* const *two = b;

        if((*one)->log_likelihood <  (*two)->log_likelihood){
                return 1;
        }else{
                return -1;
        }
}

int hierarchal_merge_motif(struct motif_list* ml,struct fhmm* fhmm,int best_i,int best_j)
{
        struct motif* a;
        struct motif* b;
        struct motif* new_motif = NULL;
        double** rand_freq_matrix = NULL;

        double* background_scores = NULL;
        int i,j;
        double score, best_score;
        int allowed_overlap;
        int num_random = 10000;
        int limit;
        int local_best_i;
        int local_best_j;
        /* all pairwise comparisons... */

        limit = ml->num_items;
        local_best_i = best_i;
        local_best_j = best_j;

        best_score = 1000000;
        for(i = 0; i < limit-1;i++){
                if(ml->mlist[i] && i != best_i){
                        a = ml->mlist[i];
                        for(j = i+1; j <limit;j++){
                                if(ml->mlist[j] && j != best_j){
                                        b = ml->mlist[j];


                                        allowed_overlap = set_overlap(a->len, b->len);// MACRO_MAX(6, MACRO_MIN(a->len, b->len)-1);
                                        fhmm->F_matrix[a->len][b->len] = 1000000;
                                        RUN(motif_dyn_programming(fhmm, a->freq_matrix, b->freq_matrix, a->len, b->len,allowed_overlap));
                                        score = fhmm->F_matrix[a->len][b->len];

                                        if(score < best_score){

                                                best_score = score;
                                                best_i = i;
                                                best_j = j;
                                        }
                                }
                        }
                }
        }


        if(best_i != local_best_i  || best_j != local_best_j){
                a = ml->mlist[best_i];
                b = ml->mlist[best_j];
                if(best_score != 0.0){
                        allowed_overlap = set_overlap(a->len, b->len);
                        RUNP(rand_freq_matrix = galloc(rand_freq_matrix, a->len, fhmm->L, 0.0));
                        for(i = 0; i < a->len;i++){
                                for(j = 0; j < fhmm->L;j++){
                                        rand_freq_matrix[i][j] = a->freq_matrix[i][j];
                                }
                        }
                        MMALLOC(background_scores, sizeof(double) * num_random);

                        for(j = 0; j < num_random;j++){
                                RUN(randomize_freq_matrix(rand_freq_matrix, a->len,fhmm->L));


                                allowed_overlap = set_overlap(a->len, b->len);
                                fhmm->F_matrix[a->len][b->len] = 1000000;
                                RUN(motif_dyn_programming(fhmm,rand_freq_matrix,b->freq_matrix, a->len, b->len,allowed_overlap));
                                background_scores[j] = fhmm->F_matrix[a->len][b->len];
                                //fprintf(stdout,"%f\n",background_scores[j] );
                        }
                        //fprintf(stdout,"\n");
                        //exit(0);
                        gfree(rand_freq_matrix);
                        qsort(background_scores, num_random, sizeof(double),double_cmp);

                        i = 0;

                        for(j = 0; j < num_random;j++){
                                if(background_scores[j] > best_score ){
                                        break;
                                }
                                i++;
                                if(i > 1001){
                                        break;
                                }

                        }
                        //if(c < 10){

                        //}
                        //fprintf(stdout,"there are %d scores better than %f (e.g. %f)\n",i,best_score,background_scores[i+1]);

                        MFREE(background_scores);
                }else{
                        i = 0;
                }

                if(i <= 1000){
                        a = ml->mlist[best_i];
                        b = ml->mlist[best_j];

                        allowed_overlap = set_overlap(a->len, b->len);

                        LOG_MSG("Merging %d %d  len:%d %d allowed over:%d score:%f",best_i, best_j, a->len, b->len,   allowed_overlap,best_score);
                        RUN(motif_dyn_programming(fhmm , a->freq_matrix, b->freq_matrix, a->len, b->len,allowed_overlap));


                        MMALLOC(new_motif, sizeof(struct motif));
                        new_motif->L = a->L;
                        new_motif->dust_score = 0.0;
                        new_motif->len = 0;
                        new_motif->log_likelihood = 0.0;
                        new_motif->count = 0;
                        new_motif->count_matrix = NULL;
                        new_motif->freq_matrix = NULL;

                        //RUNP(new_motif = init_paraclu_cluster());
                        RUN(traceback_and_merge_motifs(fhmm,a,b,new_motif));
                        //RUN(add_motif_count_matrix_based_on_hits(sb,new_motif));
                        //RUN(calculate_information_content( new_motif->freq_matrix , fhmm->background, new_motif->len, fhmm->L, &new_motif->total));
                        //  RUN(calculate_log_likelihood_based_on_hits(fhmm_log,sb,p,sa,sb_temp));
                        //fprintf(stdout,"merging? %f\t%f\t%f ", a->total,p->total,new_motif->total);
                        //if(new_motif->total >= MACRO_MIN(a->total , p->total)){//}  , b) (a->total + p->total) / 2.0){
                        free_motif(a);
                        free_motif(b);
                        ml->mlist[best_i] = new_motif;
                        ml->mlist[best_j] = NULL;
                        RUN(hierarchal_merge_motif( ml,fhmm,best_i,best_j  ));
                }
        }

        return OK;
ERROR:
        return FAIL;
}

int set_overlap(int len_a, int len_b)
{
        int c;
        c = MACRO_MAX(len_a, len_b);
        if(c >= 128){
                return MACRO_MIN(len_a, len_b);
        }
        return MACRO_MAX(6, MACRO_MIN(len_a, len_b)-1);
}

int motif_dyn_programming(struct fhmm* fhmm, double** e_a, double** e_b, int len_a, int len_b,int min_aln_len)
{
        //double** e = fhmm->e;
        double** m;
        double** t;
        int i,j,c;
        double x,y;
        m = fhmm->F_matrix;
        t = fhmm->B_matrix;
        /* Step 1: fill matrix */

        m[0][0] = 0.0f;




        j = len_a - min_aln_len;

        for(i = 1; i < j+1;i++){
                m[i][0] = 0.0f;
                t[i][0] = 1;
        }

        for(i = j+1; i < len_a+1;i++){
                m[i][0] = 100000.0f;
                t[i][0] = 1;
        }

        i = len_b - min_aln_len;
        for(j = 1; j < i+1;j++){
                m[0][j] = 0.0f;
                t[0][j] = 2;
        }

        for(j = i+1; j < len_b+1;j++){
                m[0][j] = 100000.0f;
                t[0][j] = 2;
        }

        for(i = 1; i < len_a+1;i++){
                for(j = 1; j < len_b+1;j++){
                        x = 0.0f;
                        y = 0.0f;
                        /* Euclidian Distance  */
                        /*for(c = 0; c < fhmm->L;c++){
                                y =  (e_a[i-1][c] -  e_b[j-1][c]) ;
                                x += y*y;
                                }*/

                        /* KL divergence   */
                        for(c = 0; c < fhmm->L;c++){
                                x += e_a[i-1][c] * log2f(e_a[i-1][c] / e_b[j-1][c]);
                                y += e_b[j-1][c] * log2f(e_b[j-1][c] / e_a[i-1][c]);
                        }
                        x = ((x + y)) / 2.0f;
                        /* logistic function & re-scaled to 0 - 1.0  */
                        x = (expf(x) / ( 1.0f + expf(x)) - 0.5f) *2.0;
                        m[i][j] = x  + m[i-1][j-1];
                        //m[i][j] =  sqrtf(x);
                        t[i][j] = 0;
                }
        }

        /* Fill transition for last column / row */

        /* Disallow all alignments with fewer than 6 aligned positions*/

        i =  len_a;
        for(j = 1; j < len_b+1;j++){
                if(m[i][j-1] < m[i][j]){
                        m[i][j] = m[i][j-1];
                        t[i][j] = 2;
                }
        }
        j = len_b;
        for(i = 1; i < len_a+1;i++){
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

int randomize_freq_matrix(double** m,int len, int L)
{
        int i,j,c;
        double tmp;

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

        /*for(j = 0; j < L;j++){
                fprintf(stdout,"%c  ","ACGT"[j]);
                for(i = 0; i < len;i++){
                        fprintf(stdout," %0.3f",m[i][j]);
                }
                fprintf(stdout,"\n");
        }

        fprintf(stdout,"\n");*/
        /*double sum = 0.0f;
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
        return OK;
//ERROR:
//        exit(0);
}


int traceback_and_merge_motifs(struct fhmm* fhmm,struct motif* a,struct motif* b,struct motif* new_motif)
{
        int i,j,c;
        int new_len;
        int pos_a = 0;
        int pos_b = 0;

        uint8_t* path = NULL;
        double** t;
        double sum;


        //float** new_count_matrix = NULL;

        ASSERT(fhmm != NULL, "No fhmm model");
        ASSERT(a != NULL,"No A motif");
        ASSERT(b != NULL,"No B motif");

        t = fhmm->B_matrix;

        MMALLOC(path, sizeof(uint8_t) * (a->len + b->len +1));

        /*fprintf(stdout,"Dyn matrix\n");
        for(i = 0; i < a->len+1;i++){
                for(j = 0; j < b->len+1;j++){

                        fprintf(stdout," %2.2f",fhmm->F_matrix[i][j]);
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
        fprintf(stdout,"\n");*/

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
        new_motif->len = new_len;
        new_motif->count_matrix = galloc(new_motif->count_matrix,new_motif->len, new_motif->L,0.0);
        new_motif->freq_matrix = galloc(new_motif->freq_matrix,new_motif->len, new_motif->L,0.0);
        new_motif->count = a->count + b->count;
        //new_count_matrix = malloc_2d_float(new_count_matrix, new_len, fhmm->L, 0.0f);
        //MMALLOC(new_state_sequence, sizeof(int) * (new_len));
        j = 0;
        for(i = new_len-1; i >=0;i--){
                //fprintf(stdout, "%d\t",path[i]);
                switch (path[i]) {
                case 0: {

                        //for(c = 0;c < fhmm->L;c++){
                        //        fprintf(stdout," %*.2f",3,a->count_matrix[pos_a][c]);// fhmm->e[a->state_sequence[pos_a]][c]);
                        //}
                        //for(c = 0;c < fhmm->L;c++){
                        //        fprintf(stdout," %*.2f",3,b->count_matrix[pos_b][c]);// fhmm->e[b->state_sequence[pos_b]][c]);
                        //}


                        for(c = 0;c < fhmm->L;c++){
                                new_motif->count_matrix[j][c] = a->count_matrix[pos_a][c] + b->count_matrix[pos_b][c];
                        }
                        pos_a++;
                        pos_b++;
                        break;
                }
                case 1: {
                        //for(c = 0;c < fhmm->L;c++){
                        //       fprintf(stdout," %*.2f",3,a->count_matrix[pos_a][c]);// fhmm->e[a->state_sequence[pos_a]][c]);
                        //}
                        //for(c = 0;c < fhmm->L;c++){
                        //        fprintf(stdout," %*.2f",3,0.0);
                        //}
                        //new_state_sequence[j] = a->state_sequence[pos_a];
                        for(c = 0;c < fhmm->L;c++){
                                new_motif->count_matrix[j][c] = a->count_matrix[pos_a][c];
                                //new_freq_matrix[j][c] = a->freq_matrix[pos_a][c];
                        }
                        pos_a++;
                        break;
                }
                case 2: {
                        //for(c = 0;c < fhmm->L;c++){
                        //        fprintf(stdout," %*.2f",3,0.0);
                        //}
                        //for(c = 0;c < fhmm->L;c++){
                        //        fprintf(stdout," %*.2f",3,b->count_matrix[pos_b][c]);// fhmm->e[b->state_sequence[pos_b]][c]);
                        //        }
                        //new_state_sequence[j] = b->state_sequence[pos_b];
                        for(c = 0;c < fhmm->L;c++){
                                new_motif->count_matrix[j][c] = b->count_matrix[pos_b][c];
                                //new_freq_matrix[j][c] = b->freq_matrix[pos_b][c];
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
        //fprintf(stdout, "\n");
        for(i = 0; i < new_motif->len;i++){
                sum = 0.0;
                for(j = 0; j < new_motif->L;j++){
                        sum += new_motif->count_matrix[i][j] + fhmm->background[j];
                }

                for(j = 0; j < new_motif->L;j++){
                        new_motif->freq_matrix[i][j] = (new_motif->count_matrix[i][j]+ fhmm->background[j]) / sum;
                }



        }


        //}

        /*new_motif->len = new_len;
        new_motif->state_sequence = new_state_sequence;
        new_motif->num_hits = a->num_hits + b->num_hits;
        MMALLOC(new_motif->hits, sizeof(struct hit*) *(a->num_hits + b->num_hits));

        for(i= 0; i < new_motif->num_hits;i++){
                new_motif->hits[i] = NULL;
                MMALLOC(new_motif->hits[i], sizeof(struct hit));
        }

        new_motif->num_hits = 0;


        */






        //done.
                //        fprintf(stdout,"%d %d\n", a->hits[j]->seq,a->hits[j]->pos);

        /* fill with hits from a and b */
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
        //exit(0);
        MFREE(path);
        return OK;
ERROR:
        MFREE(path);
        return FAIL;
}


int double_cmp(const void *a, const void *b)
{
    const double *ia = (const double *)a;
    const double *ib = (const double *)b;
    if(*ia > *ib ){
            return 1;
    }else{
            return -1;
    }
}



int dust_sequence( uint8_t* seq, int len,double* score)
{
        int i,c;
        int key = 0;
        double triplet[64];
        double s = 0.0;
        for(i = 0;i < 64;i++){
                triplet[i] = 0.0;
        }
        c = 0;
        key = ((seq[c]& 0x3 ) << 2 )|  (seq[c+1]& 0x3);
        c+= 2;

        for(i = c;i < len;i++){

                key = key << 2 | (seq[i]& 0x3);
                triplet[key & 0x3F]++; /* 6 bits = 3nucleotides */
                c++;
        }
        s = 0.0;
        for(i = 0;i < 64;i++){

                s+= (triplet[i] * (triplet[i] -1.0)) / 2.0;
                //triplet[i] = 0.0;
        }
        s = s / (double)(len-3);
        *score = s;
        return OK;
}


int write_meme_output(char* filename,struct motif_list* m, struct fhmm*  fhmm)
{
        struct motif* p =  NULL;
        FILE* fptr = NULL;
        int i,j,c;
        double sum;
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
                        p = m->mlist[i];
                        if(p){
                                fprintf(fptr,"MOTIF M%03d\n", i+1);
                                fprintf(stdout,"MOTIF M%03d\n", i+1);
                                fprintf(fptr,"letter-probability matrix: alength= %d w= %d nsites= %d E= 4.1e-009\n",fhmm->L, p->len, p->count);
                                for(j = 0 ; j < p->len;j++){
                                        sum = 0.0f;
                                        for(c = 0; c < fhmm->L;c++){
                                                sum += p->count_matrix[j][c];
                                        }
                                        fprintf(fptr,"%0.6f", p->count_matrix[j][0]/sum);
                                        fprintf(stdout,"%0.6f", p->count_matrix[j][0]);
                                        for(c = 1; c < fhmm->L;c++){
                                                fprintf(fptr," %0.6f", p->count_matrix[j][c]/sum);
                                                fprintf(stdout," %0.6f", p->count_matrix[j][c]);
                                        }
                                        fprintf(fptr,"\n");
                                        fprintf(stdout,"\tsum:%f\n",sum);
                                }
                                fprintf(fptr,"\n");
                                fprintf(stdout,"\n");
                        }
                }

                fclose(fptr);

        }

        return OK;
ERROR:
        return FAIL;
}




/*****************************************************************
 * 3. Standard statistical tests.
 *****************************************************************/

/* Function:  esl_stats_GTest()
 * Synopsis:  Calculates a G-test on 2 vs. 1 binomials.
 *
 * Purpose:   In experiment a, we've drawn <ca> successes in <na> total
 *            trials; in experiment b, we've drawn <cb> successes in
 *            <nb> total trials. Are the counts different enough to
 *            conclude that the two experiments are different? The
 *            null hypothesis is that the successes in both experiments
 *            were drawn from the same binomial distribution with
 *            per-trial probability $p$. The tested hypothesis is that
 *            experiments a,b have different binomial probabilities
 *            $p_a,p_b$. The G-test is a log-likelihood-ratio statistic,
 *            assuming maximum likelihood values for $p,p_a,p_b$.
 *            $2G$ is distributed approximately as $X^2(1)$,
 *            %"X" is "Chi"
 *            which we use to calculate a P-value for the G statistic.
 *
 * Args:      ca    - number of positives in experiment a
 *            na    - total number in experiment a
 *            cb    - number of positives in experiment b
 *            nb    - total number in experiment b
 *            ret_G - RETURN: G statistic, a log likelihood ratio, in nats
 *            ret_P - RETURN: P-value for the G-statistic
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      Archive1999/0906-sagescore/sagescore.c
 */
int esl_stats_GTest(int ca, int na, int cb, int nb, double *ret_G, double *ret_P)
{
        double a,b,c,d,n;
        double G = 0.;

        a = (double) ca;
        b = (double) (na - ca);
        c = (double) cb;
        d = (double) (nb - cb);
        n = (double) na+nb;

        /* Yes, the calculation here is correct; algebraic
         * rearrangement of the log-likelihood-ratio with
         * p_a = ca/na, p_b = cb/nb, and p = (ca+cb)/(na+nb).
         * Guard against 0 probabilities; assume 0 log 0 => 0.
         */
        if (a   > 0.) G  = a * log(a);
        if (b   > 0.) G += b * log(b);
        if (c   > 0.) G += c * log(c);
        if (d   > 0.) G += d * log(d);
        if (n   > 0.) G += n * log(n);
        if (a+b > 0.) G -= (a+b) * log(a+b);
        if (c+d > 0.) G -= (c+d) * log(c+d);
        if (a+c > 0.) G -= (a+c) * log(a+c);
        if (b+d > 0.) G -= (b+d) * log(b+d);

        *ret_G = G;
        if(G > 0.0){
                RUN(esl_stats_IncompleteGamma( 0.5, G, NULL, ret_P));
        }else{
                *ret_P = 1.0;
        }

        return OK;
ERROR:
        return FAIL;
}


/* Function: esl_stats_IncompleteGamma()
 * Synopsis: Calculates the incomplete Gamma function.
 *
 * Purpose:  Returns $P(a,x)$ and $Q(a,x)$ where:
 *
 *           \begin{eqnarray*}
 *             P(a,x) & = & \frac{1}{\Gamma(a)} \int_{0}^{x} t^{a-1} e^{-t} dt \\
 *                    & = & \frac{\gamma(a,x)}{\Gamma(a)} \\
 *             Q(a,x) & = & \frac{1}{\Gamma(a)} \int_{x}^{\infty} t^{a-1} e^{-t} dt\\
 *                    & = & 1 - P(a,x) \\
 *           \end{eqnarray*}
 *
 *           $P(a,x)$ is the CDF of a gamma density with $\lambda = 1$,
 *           and $Q(a,x)$ is the survival function.
 *
 *           For $x \simeq 0$, $P(a,x) \simeq 0$ and $Q(a,x) \simeq 1$; and
 *           $P(a,x)$ is less prone to roundoff error.
 *
 *           The opposite is the case for large $x >> a$, where
 *           $P(a,x) \simeq 1$ and $Q(a,x) \simeq 0$; there, $Q(a,x)$ is
 *           less prone to roundoff error.
 *
 * Method:   Based on ideas from Numerical Recipes in C, Press et al.,
 *           Cambridge University Press, 1988.
 *
 * Args:     a          - for instance, degrees of freedom / 2     [a > 0]
 *           x          - for instance, chi-squared statistic / 2  [x >= 0]
 *           ret_pax    - RETURN: P(a,x)
 *           ret_qax    - RETURN: Q(a,x)
 *
 * Return:   <eslOK> on success.
 *
 * Throws:   <eslERANGE> if <a> or <x> is out of accepted range.
 *           <eslENOHALT> if approximation fails to converge.
 */
int
esl_stats_IncompleteGamma(double a, double x, double *ret_pax, double *ret_qax)
{
        int    iter;			/* iteration counter */
        double pax;			/* P(a,x) */
        double qax;			/* Q(a,x) */
        ASSERT(a > 0.0, "a must be > 0");
        ASSERT(x >= 0.0, "x must be >= 0");
        //if (a <= 0.) ESL_EXCEPTION(eslERANGE, "esl_stats_IncompleteGamma(): a must be > 0");
        //if (x <  0.) ESL_EXCEPTION(eslERANGE, "esl_stats_IncompleteGamma(): x must be >= 0");

        /* For x > a + 1 the following gives rapid convergence;
         * calculate Q(a,x) = \frac{\Gamma(a,x)}{\Gamma(a)},
         * using a continued fraction development for \Gamma(a,x).
         */
        if (x > a+1)
        {
                double oldp;		/* previous value of p    */
                double nu0, nu1;		/* numerators for continued fraction calc   */
                double de0, de1;		/* denominators for continued fraction calc */

                nu0 = 0.;			/* A_0 = 0       */
                de0 = 1.;			/* B_0 = 1       */
                nu1 = 1.;			/* A_1 = 1       */
                de1 = x;			/* B_1 = x       */

                oldp = nu1;
                for (iter = 1; iter < 100; iter++)
                {
                        /* Continued fraction development:
                         * set A_j = b_j A_j-1 + a_j A_j-2
                         *     B_j = b_j B_j-1 + a_j B_j-2
                         * We start with A_2, B_2.
                         */
                        /* j = even: a_j = iter-a, b_j = 1 */
                        /* A,B_j-2 are in nu0, de0; A,B_j-1 are in nu1,de1 */
                        nu0 = nu1 + ((double)iter - a) * nu0;
                        de0 = de1 + ((double)iter - a) * de0;
                        /* j = odd: a_j = iter, b_j = x */
                        /* A,B_j-2 are in nu1, de1; A,B_j-1 in nu0,de0 */
                        nu1 = x * nu0 + (double) iter * nu1;
                        de1 = x * de0 + (double) iter * de1;
                        /* rescale */
                        if (de1 != 0.)
                        {
                                nu0 /= de1;
                                de0 /= de1;
                                nu1 /= de1;
                                de1 =  1.;
                        }
                        /* check for convergence */
                        if (fabs((nu1-oldp)/nu1) < 1.e-7)
                        {
                                RUN(esl_stats_LogGamma(a, &qax));
                                //if ((status = esl_stats_LogGamma(a, &qax)) != eslOK) return status;
                                qax = nu1 * exp(a * log(x) - x - qax);

                                if (ret_pax != NULL) *ret_pax = 1 - qax;
                                if (ret_qax != NULL) *ret_qax = qax;
                                return OK;
                        }

                        oldp = nu1;
                }
                ERROR_MSG("esl_stats_IncompleteGamma(): series failed to converge");
        }
        else /* x <= a+1 */
        {
                double p;			/* current sum               */
                double val;		/* current value used in sum */

                /* For x <= a+1 we use a convergent series instead:
                 *   P(a,x) = \frac{\gamma(a,x)}{\Gamma(a)},
                 * where
                 *   \gamma(a,x) = e^{-x}x^a \sum_{n=0}{\infty} \frac{\Gamma{a}}{\Gamma{a+1+n}} x^n
                 * which looks appalling but the sum is in fact rearrangeable to
                 * a simple series without the \Gamma functions:
                 *   = \frac{1}{a} + \frac{x}{a(a+1)} + \frac{x^2}{a(a+1)(a+2)} ...
                 * and it's obvious that this should converge nicely for x <= a+1.
                 */
                p = val = 1. / a;
                for (iter = 1; iter < 10000; iter++)
                {
                        val *= x / (a+(double)iter);
                        p   += val;

                        if (fabs(val/p) < 1.e-7)
                        {

                                RUN(esl_stats_LogGamma(a, &pax));
                                pax = p * exp(a * log(x) - x - pax);

                                if (ret_pax != NULL) *ret_pax = pax;
                                if (ret_qax != NULL) *ret_qax = 1. - pax;
                                return OK;
                        }
                }
                ERROR_MSG("esl_stats_IncompleteGamma(): series failed to converge");
        }
        /*NOTREACHED*/
        return OK;
ERROR:
        return FAIL;
}


struct alignment* init_alignment(int e_num_seq)
{
        struct alignment* a = NULL;
        int i;

        MMALLOC(a, sizeof(struct alignment));
        a->alloc_num_seq = e_num_seq;
        a->seq_counts = NULL;
        a->sequences = NULL;
        a->num_seq = 0;
        a->motif = NULL;

        MMALLOC(a->seq_counts, sizeof(int) * a->alloc_num_seq);
        MMALLOC(a->sequences, sizeof(char*) * a->alloc_num_seq);

        for(i = 0; i < a->alloc_num_seq;i++){
                a->seq_counts[i] = 0;
                a->sequences[i] = NULL;
        }
        return a;
ERROR:
        return NULL;
}

int resize_alignment(struct alignment* a)
{
        int i,old;
        ASSERT(a!= NULL, "no alignment");
        old =  a->alloc_num_seq;
        a->alloc_num_seq =  a->alloc_num_seq << 1;

        MREALLOC(a->seq_counts, sizeof(int) * a->alloc_num_seq);
        MREALLOC(a->sequences, sizeof(char*) * a->alloc_num_seq);
        for(i = old; i < a->alloc_num_seq;i++){
                a->seq_counts[i] = 0;
                a->sequences[i] = NULL;
        }
        return OK;
ERROR:
        return FAIL;

}

void free_alignment(struct alignment* a)
{
        int i;

        if(a){
                if(a->sequences){
                        for(i = 0; i < a->num_seq;i++){
                                MFREE(a->sequences[i]);
                        }
                        MFREE(a->sequences);

                }
                if(a->seq_counts){
                        MFREE(a->seq_counts);
                }

                if(a->motif){
                        free_motif(a->motif);
                }
                MFREE(a);
        }
}

int create_motif_based_on_alignment(struct alignment* a,double* background, int L)
{
        struct motif* m = NULL;
        char* tmp = NULL;

        double sum;
        int motif_len = 0;
        int i,j;
        ASSERT(a != NULL, "No alignment");
        ASSERT(a->num_seq != 0, " Weird - alignment does not contain sequences");
        motif_len = 0;
        for(i = 0; i < a->num_seq;i++){
                motif_len = MACRO_MAX(motif_len, strlen(a->sequences[i]));
                //LOG_MSG("%d", strlen(a->sequences[i]));
        }
        /* Check if sequences are actually aligned */
        for(i = 0; i < a->num_seq;i++){
                if(motif_len != strlen(a->sequences[i])){
                        ERROR_MSG("Sequences are not aligned!");
                }
        }

        //LOG_MSG("Motif_len: %d:", motif_len);
        MMALLOC(m, sizeof(struct motif));
        m->count_matrix = NULL;
        m->freq_matrix = NULL;
        m->len = motif_len;
        m->log_likelihood = 0.0;
        m->count = 0;
        m->L = L;
        RUNP(m->count_matrix = galloc(m->count_matrix,m->len,L,0.0));
        RUNP(m->freq_matrix = galloc(m->freq_matrix,m->len,L,0.0));

        for(i = 0; i < a->num_seq;i++){
                tmp= a->sequences[i];
                m->count += a->seq_counts[i];
                for(j = 0; j < m->len;j++){
                        switch (tmp[j]) {
                        case 'A':
                        case 'a':
                                m->count_matrix[j][0] += a->seq_counts[i];

                                break;

                        case 'C':
                        case 'c':

                                m->count_matrix[j][1] += a->seq_counts[i];
                                break;
                        case 'G':
                        case 'g':

                                m->count_matrix[j][2] += a->seq_counts[i];
                                break;
                        case 'T':
                        case 't':

                                m->count_matrix[j][3] += a->seq_counts[i];
                                break;

                        default:
                                break;
                        }
                }
        }
        for(i = 0;i <  m->len;i++){

                sum = 0.0;
                for(j = 0; j < m->L;j++){
                        sum += m->count_matrix[i][j] + background[j];
                }
                for(j = 0; j < m->L;j++){
                        m->freq_matrix[i][j] = (m->count_matrix[i][j]  + background[j])/ sum;
                        //        fprintf(stdout,"%f ", m->freq_matrix[i][j]);
                }
        }
        a->motif = m;

        return OK;
ERROR:
        return FAIL;
}
