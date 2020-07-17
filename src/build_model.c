
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>

#include "tldevel.h"

#include "tllogsum.h"
#include "tlmisc.h"
#include "tlrng.h"

//#include "ihmm_seq.h"
#include "beam_sample.h"

#include "finite_hmm.h"

#include "distributions.h"

#include "run_score.h"

#include "sequence_alloc.h"
#include "sequence_io.h"
#include "sequence_struct.h"
#include "sequence_prep.h"

#include "pst_build.h"


#include "model_core.h"

#include "model_struct.h"
#include "model_io.h"
#include "fast_hmm_param.h"
//#include "init_seq_label.h"
#include "hmm_conversion.h"

#include "esl_stopwatch.h"

#include "thread_data_io.h"

#define OPT_SEED 1
#define OPT_NUM_MODELS 2
#define OPT_COMPETITIVE 3
#define OPT_SEQDB 4

struct parameters{
        char* input;
        char* in_model;
        char* seq_db;
        char* cmd_line;
        double alpha;
        double gamma;
        unsigned long seed;
        rk_state rndstate;
        struct rng_state* rng;
        int active_file;
        int competitive;
        int num_iter;
        int inner_iter;
        int local;
        int num_models;
        int num_threads;
        int num_start_states;
        int num_max_states;
        int rev;
};

static int run_build_ihmm(struct parameters* param);
static int score_all_vs_all(struct model_bag* mb, struct seq_buffer* sb, struct seqer_thread_data** td);
static int analyzescores(struct seq_buffer* sb,  struct model_bag* model_bag );
static int reset_sequence_weights(struct seq_buffer* sb, int num_models);
static int set_sequence_weights(struct seq_buffer* sb, int num_models, double temperature);
static int write_program_state(char* filename, struct model_bag* mb, struct seq_buffer* sb, struct seqer_thread_data** td, int n_thread);


static int random_score_sequences(struct seq_buffer* sb,double* background );
static int score_sequences_for_command_line_reporting(struct parameters* param);

static int init_num_state_array(int* num_state_array, int maxK, int len, struct parameters* param, double mean);
static int print_help(char **argv);
static int free_parameters(struct parameters* param);

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        //print_program_header(argv, "Build HDPHMM model(s).");

        MMALLOC(param, sizeof(struct parameters));
        param->input = NULL;
        param->in_model = NULL;
        param->seq_db = NULL;
        param->cmd_line = NULL;
        //param->tmp_file_A = NULL;
        //param->tmp_file_B = NULL;
        param->num_threads = 1;
        param->num_start_states = 0;
        param->local = 0;
        param->rev = 0;
        param->active_file = 0;
        param->num_iter = 2000;
        param->inner_iter = 100;
        param->alpha = IHMM_PARAM_PLACEHOLDER;
        param->gamma = IHMM_PARAM_PLACEHOLDER;
        param->seed = 0;
        param->num_models = 1;
        param->competitive = 0;
        param->num_max_states = 1000;
        param->rng = NULL;
        while (1){
                static struct option long_options[] ={
                        {"in",required_argument,0,'i'},
                        {"seqdb",required_argument,0,OPT_SEQDB},
                        {"out",required_argument,0,'o'},
                        {"states",required_argument,0,'s'},
                        {"local",no_argument,0,'l'},
                        {"nthreads",required_argument,0,'t'},
                        {"nmodels",required_argument,0, OPT_NUM_MODELS},
                        {"niter",required_argument,0,'n'},
                        {"model",required_argument,0,'m'},
                        {"alpha",required_argument,0,'a'},
                        {"gamma",required_argument,0,'g'},
                        {"seed",required_argument,0,OPT_SEED},
                        {"competitive",no_argument,0,OPT_COMPETITIVE},
                        {"rev",0,0,'r'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"rhi:o:t:n:m:s:la:g:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case OPT_SEQDB:
                        param->seq_db = optarg;
                        break;
                case OPT_COMPETITIVE:
                        param->competitive =1;
                        break;
                case OPT_SEED:
                        param->seed = atoi(optarg);
                        break;
                case OPT_NUM_MODELS:
                        param->num_models = atoi(optarg);
                        break;
                case 'a':
                        param->alpha = atof(optarg);
                        break;
                case 'g':
                        param->gamma = atof(optarg);
                        break;
                case 'r':
                        param->rev = 1;
                        break;
                case 'l':
                        param->local = 1;
                        break;
                case 'i':
                        param->input = optarg;
                        break;
                case 'o':
                        //param->output = optarg;
                        break;
                case 'm':
                        param->in_model = optarg;
                        break;
                case 's':
                        param->num_start_states = atoi(optarg);
                        break;
                case 'n':
                        param->num_iter = atoi(optarg);
                        break;
                case 't':
                        param->num_threads = atoi(optarg);
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
                ERROR_MSG("No input file! use --in <blah.fa>");
        }else{
                if(!my_file_exists(param->input)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->input);
                }
        }

        if(!param->in_model){
                RUN(print_help(argv));
                ERROR_MSG("No output file! use --in <blah.fa>");
        }else{
                if(my_file_exists(param->in_model)){

                        LOG_MSG("Will continue training models in %s.",param->in_model);
                        if(param->input){
                                WARNING_MSG("Sequences in %s will be ignored.", param->input);
                                WARNING_MSG("Instead models will continue to be trained on the sequences stored alongside the model in %s", param->in_model );
                        }
                }else{
                        LOG_MSG("Will construct new model(s).");
                }
        }

        if(!param->seq_db){
                RUN(print_help(argv));
                ERROR_MSG("No seqDB use --seqdb <blah.fa>");

        }else{
                if(!my_file_exists(param->seq_db)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->seq_db);
                }
        }

        if(param->num_models < 1){
                RUN(print_help(argv));
                ERROR_MSG("To few models! use -nmodels <1+>");
        }

        if(param->seed){
                RUNP(param->rng = init_rng(param->seed));
                rk_seed(param->seed, &param->rndstate);
        }else{
                RUNP(param->rng = init_rng(0));
                rk_randomseed(&param->rndstate);
        }

        RUN(make_cmd_line(&param->cmd_line,argc,argv));

        RUN(run_build_ihmm(param));

        RUN(score_sequences_for_command_line_reporting(param));

        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}

int run_build_ihmm(struct parameters* param)
{
        struct fast_param_bag* ft_bag = NULL;
        struct model_bag* model_bag = NULL;
        struct seq_buffer* sb = NULL;
        struct seqer_thread_data** td = NULL;

        int i;

        ASSERT(param!= NULL, "No parameters found.");
        init_logsum();

        /* If we have saved a model continue from there */
        if(my_file_exists(param->in_model)){
                /* PROBABLY need to re-alloc num_state_array */
                //MMALLOC(num_state_array, sizeof(int)* param->num_models);
                RUNP(model_bag = read_model_bag_hdf5(param->in_model));
                RUNP(sb = get_sequences_from_hdf5_model(param->in_model,IHMM_SEQ_READ_ALL));

                MMALLOC(sb->num_state_arr, sizeof(int) * model_bag->num_models);
                for(i = 0; i < model_bag->num_models;i++){
                        sb->num_state_arr[i] = model_bag->models[i]->num_states;
                }

                RUNP(td = read_thread_data_to_hdf5(param->in_model));

                RUN(convert_ihmm_to_fhmm_models(model_bag));
                exit(0);
        }else{                  /* New run - start from the beginning */
                /* create pst search model  */
                RUN(create_pst_model(param->rng, param->input,param->seq_db, param->in_model,0.000001, 0.01, 20.0));
/* Step one read in sequences */
                LOG_MSG("Loading sequences.");
                //RUNP(rng = init_rng(0));

                RUN(read_sequences_file(&sb, param->input));


                RUN(prep_sequences(sb,param->rng, param->num_models,param->num_start_states,0.0));

                //RUNP(sb = load_sequences(param->input,&param->rndstate));

                //sb = concatenate_sequences(sb);
                //RUN(label_seq_based_on_random_fhmm(sb, initial_states, 0.3));
                //LOG_MSG("%d models ",param->num_models);
                /* This function sets the number of starting states to max seq len /2 +- 25  */

                //LOG_MSG("MAXK: %d", param->num_max_states);
                RUNP(model_bag = alloc_model_bag(sb->num_state_arr, sb->L, param->num_models,MAX_NUM_STATES , param->seed));

                //label_ihmm_sequences_based_on_guess_hmm(struct seq_buffer *sb, int k, float alpha)
                /* New label sequences  */
                /* We need to label sequences somehow so that we can
                 * set the u arrays properly based on the initial
                 * model guess */
                //model_bag->max_num_states = 0;
                LOG_MSG("%d ",sb->num_seq);
                for(i = 0; i < model_bag->num_models;i++){

                        RUN(fill_counts(model_bag->models[i], sb,i));

                        //model_bag->max_num_states = MACRO_MAX(model_bag->max_num_states,model_bag->models[i]->num_states);
                }

                //RUN(check_labels(sb,model_bag->num_models ));
                /* Set seed in sequence buffer */
                //sb->seed = rk_ulong(&model_bag->rndstate);
                //rk_seed(sb->seed, &sb->rndstate);

                RUN(set_model_hyper_parameters(model_bag, param->alpha, param->gamma));


                /* Allocating thread structure. */
                RUNP(td = create_seqer_thread_data(&param->num_threads, (sb->max_len+2)  ,model_bag->max_num_states, &model_bag->rndstate, THREAD_DATA_BEAM));
        }

        LOG_MSG("Will use %d threads.", param->num_threads);
        //if((pool = thr_pool_create(param->num_threads,param->num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");

        RUNP(ft_bag = alloc_fast_param_bag(model_bag->num_models, sb->num_state_arr, sb->L));
        //RUNP(ft = alloc_fast_hmm_param(initial_states,sb->L));
        /* fill background of first fast hmm param struct  */
        //RUN(fill_background_emission(ft_bag->fast_params[0], sb));
        /* Now copy remaining 1:N first fast hmm param structs */
        /*for(i = 1; i < ft_bag->num_models;i++){
                for(j = 0; j < sb->L;j++){
                        ft_bag->fast_params[i]->background_emission[j] = ft_bag->fast_params[0]->background_emission[j];
                }
                }*/
        //LOG_MSG("lasty max state: %d",ft_bag->max_last_state);
        /* Do a random score */
        RUN(random_score_sequences(sb, model_bag->models[0]->background  ));

        /* Main function */
        int outer_iter = param->num_iter / param->inner_iter;

        DECLARE_TIMER(n);
        //RUN(convert_ihmm_to_fhmm_models(model_bag));
        for(i = 0; i < outer_iter;i++){
                //LOG_MSG("Outer: %d %d", i, param->inner_iter);
/* run inner iter beam sampling iterations */
                LOG_MSG("Start beam");
                START_TIMER(n);
                RUN(run_beam_sampling(model_bag,ft_bag, sb,td, param->inner_iter, param->num_threads));
                STOP_TIMER(n);
                GET_TIMING(n);

                LOG_MSG("Start scoring");
                START_TIMER(n);
                /* convert to fhmm */
                RUN(convert_ihmm_to_fhmm_models(model_bag));
                /* score */
                RUN(score_all_vs_all(model_bag,sb,td));
                /* analyzescores */
                RUN(analyzescores(sb, model_bag));
                /* need to reset weights before writing models to disk!  */
                if(param->competitive){ /* competitive training */
                        set_sequence_weights(sb,  model_bag->num_models, 2.0 / log10f( (float) (i+1) + 1.0));
                }else{          /* reset sequence weights.  */
                        reset_sequence_weights(sb, model_bag->num_models);
                }
                STOP_TIMER(n);
                GET_TIMING(n);

                /* write temporary results */
                LOG_MSG("Writing model");
                START_TIMER(n);
                write_program_state(param->in_model, model_bag, sb, td, param->num_threads);
                STOP_TIMER(n);
                GET_TIMING(n);

        }

        DESTROY_TIMER(n);

        /* Write results */
        RUN(convert_ihmm_to_fhmm_models(model_bag));
        //RUN(score_all_vs_all(model_bag,sb,td));
        RUN(write_model_bag_hdf5(model_bag,param->in_model));
        //RUN(add_annotation(param->in_model, "seqwise_model_cmd", param->cmd_line));
        RUN(add_sequences_to_hdf5_model(param->in_model, sb,  model_bag->num_models));
        RUN(write_thread_data_to_hdf5(param->in_model, td, param->num_threads, sb->max_len+2, model_bag->max_num_states));
        //RUN(write_thread_data_to_)
        //RUN(write_model(model, param->output));
        /*for(i = 0; i < model_bag->num_models;i++){
                char buffer[BUFFER_LEN];
                snprintf(buffer, BUFFER_LEN, "%s_%d.h5", param->output,i);
                RUN(write_model_hdf5(model_bag->models[i],buffer));
                RUN(add_annotation(buffer, "wims_model_cmd", param->cmd_line));
                RUN(add_background_emission(buffer,ft_bag->fast_params[0]->background_emission,ft_bag->fast_params[0]->L));

                RUN(add_sequences_to_hdf5_model(buffer, sb, i));

                RUN(run_build_fhmm_file(buffer,/models/m10));
                }*/
//        RUN(add_annotation)

        //RUN(add_sequences_to_hdf5_model(param->output, sb));
        //RUN(print_states_per_sequence(sb));
        //RUN(write_ihmm_sequences(sb,"test.lfs","testing"));
        //sb, num thread, guess for aplha and gamma.. iterations.

        free_ihmm_sequences(sb);
        free_model_bag(model_bag);
        free_fast_param_bag(ft_bag);
        free_seqer_thread_data(td);
        //thr_pool_destroy(pool);
        //MFREE(num_state_array);
        return OK;
ERROR:
        free_ihmm_sequences(sb);
        free_model_bag(model_bag);
        free_fast_param_bag(ft_bag);
        free_seqer_thread_data(td);
        //thr_pool_destroy(pool);
        //MFREE(num_state_array);
        return FAIL;
}

int write_program_state(char* filename, struct model_bag* mb, struct seq_buffer* sb, struct seqer_thread_data** td, int n_thread)
{

        RUN(convert_ihmm_to_fhmm_models(mb));
        //RUN(score_all_vs_all(model_bag,sb,td));
        RUN(write_model_bag_hdf5(mb,filename));
        //RUN(add_annotation(param->in_model, "seqwise_model_cmd", param->cmd_line));
        RUN(add_sequences_to_hdf5_model(filename, sb, mb->num_models));
        RUN(write_thread_data_to_hdf5(filename, td, n_thread, sb->max_len+2, mb->max_num_states));
        return OK;
ERROR:
        return FAIL;

}

int set_sequence_weights(struct seq_buffer* sb, int num_models, double temperature)
{
        int i,j;
        double sum;
        double prior;
        /* assume model prior is 1 / num model */
/* re-set!!!  */
        prior = prob2scaledprob(1.0 / (double) num_models);
        for(i = 0; i < sb->num_seq;i++){
                //fprintf(stdout,"Seq%d\t",i);
                sum = prob2scaledprob(0.0);
                for(j = 0; j < num_models;j++){
                        sb->sequences[i]->score_arr[j] = sb->sequences[i]->score_arr[j] + prior;
                        sum = logsum(sum, sb->sequences[i]->score_arr[j]);


                }
                //fprintf(stdout,"Seq%d\t",i);
                for(j = 0; j < num_models;j++){
                       sb->sequences[i]->score_arr[j] = sb->sequences[i]->score_arr[j] - sum;
                       //fprintf(stdout,"%f ", scaledprob2prob(sb->sequences[i]->score_arr[j]));
                }
                //fprintf(stdout,"\n");
                sum = prob2scaledprob(0.0);
                for(j = 0; j < num_models;j++){
                        sb->sequences[i]->score_arr[j] = sb->sequences[i]->score_arr[j] * (1.0 / temperature);
                        sum = logsum(sum, sb->sequences[i]->score_arr[j]);

                }
                //fprintf(stdout,"Seq%d\t",i);
                for(j = 0; j < num_models;j++){
                        sb->sequences[i]->score_arr[j] = sb->sequences[i]->score_arr[j] - sum;
                        //fprintf(stdout,"%f ", scaledprob2prob(sb->sequences[i]->score_arr[j]));
                        sb->sequences[i]->score_arr[j] = 1.0 +  (double) num_models *  scaledprob2prob(sb->sequences[i]->score_arr[j]);
                }
                //fprintf(stdout,"\n");



        }
        return OK;
}


int reset_sequence_weights(struct seq_buffer* sb, int num_models)
{
        int i,j;
        /* re-set!!!  */
        for(i = 0; i < sb->num_seq;i++){
                //fprintf(stdout,"Seq%d\t",i);
                for(j = 0; j < num_models;j++){
                        sb->sequences[i]->score_arr[j] = 1.0;

                }
        }

        return OK;
}

int analyzescores(struct seq_buffer* sb, struct model_bag* model_bag)
{
        double s0,s1,s2;
        int i,j;
        //int max_print;
        int num_models = model_bag->num_models;
        ASSERT(sb!= NULL, "No sequences");
        /*max_print = MACRO_MIN(5, sb->num_seq);
        for(i = 0; i < max_print;i++){
                fprintf(stdout,"Seq%d\t",i);
                for(j = 0; j < num_models;j++){
                        fprintf(stdout,"%f ", sb->sequences[i]->score_arr[j]);

                }
                fprintf(stdout,"\n");

                }*/

        for(j = 0; j < num_models;j++){
                s0 = 0.0;
                s1 = 0.0;
                s2 = 0.0;
                for(i = 0; i < sb->num_seq;i++){
                        s0 += 1.0;
                        s1 += sb->sequences[i]->score_arr[j];
                        s2 += sb->sequences[i]->score_arr[j] * sb->sequences[i]->score_arr[j];

                }

                s2 = sqrt((s0 * s2 - s1 * s1)/ (s0 * (s0 -1.0)));
                s1 = s1 / s0 ;
                fprintf(stdout,"Model %d:\t%f\t%f\t(%d states)  alpha = %f, gamma = %f\n",j, s1,s2  ,model_bag->models[j]->num_states, model_bag->models[j]->alpha ,model_bag->models[j]->gamma);
        }
        return OK;
ERROR:
        return FAIL;
}

/* Score all sequences using all models */
int score_all_vs_all(struct model_bag* mb, struct seq_buffer* sb, struct seqer_thread_data** td)
{
        int i,j,c;
        int num_threads = td[0]->num_threads;
        int run;

        ASSERT(sb != NULL, "No sequences");
        ASSERT(mb != NULL, "No models");
        c = 0;

        for(run = 0; run < mb->num_models;run+= num_threads){
                j = 0;
                for(i = 0; i < num_threads;i++){
                        td[i]->thread_ID = c;
                        td[i]->fhmm = mb->finite_models[c];
                        td[i]->sb = sb;
                        j++;
                        c++;
                        if(c == mb->num_models){
                                break;
                        }
                }
#ifdef HAVE_OPENMP
                omp_set_num_threads( MACRO_MIN(num_threads,j));
#pragma omp parallel shared(td) private(i)
                {
#pragma omp for schedule(dynamic) nowait
#endif
                        for(i = 0; i < j;i++){
                                do_score_sequences_per_model(td[i]);
                        }
#ifdef HAVE_OPENMP
                }
#endif
        }
        return OK;
ERROR:
        return FAIL;
}

int random_score_sequences(struct seq_buffer* sb,double* background )
{
        struct ihmm_sequence* s;
        double* back = NULL;
        int expected_len;
        int i;

        ASSERT(sb!=NULL, "No sequences");
        ASSERT(background != NULL, "No background");

        MMALLOC(back, sizeof(double) * sb->L);
        for(i= 0; i < sb->L;i++){
                back[i] = prob2scaledprob(background[i]);
        }
        expected_len = 0;
        for(i = 0; i < sb->num_seq;i++){
                expected_len += sb->sequences[i]->seq_len;
        }
        expected_len = expected_len / sb->num_seq;
        LOG_MSG("Average sequence length: %d",expected_len);

        for(i =0; i < sb->num_seq;i++){
                s = sb->sequences[i];
                RUN(random_model_score(s->seq_len, &s->r_score));// , s->seq, s->seq_len,expected_len));
        }
        MFREE(back);
        //exit(0);
        return OK;

ERROR:
        MFREE(back);
        return FAIL;
}


/* After training a model it would be nice to know the logodds scores
 * of a number of sequences to see if something was found */

int score_sequences_for_command_line_reporting(struct parameters* param)
{
        struct model_bag* model_bag = NULL;
        struct seqer_thread_data** td = NULL;
        //struct thr_pool* pool = NULL;
        struct seq_buffer* sb = NULL;
        struct fhmm* fhmm = NULL;


        //double** all_scores = NULL;
        //double s1,s2;
        double best_score;
        //double sum;
        double max_likelihood;

        double BIC;
        double total_seq_len;
        int max = 1000;
        int c;
        int i;
        int limit;

        //char buffer[BUFFER_LEN];



        LOG_MSG("Reading models.");
        RUNP(model_bag = read_model_bag_hdf5(param->in_model));
        LOG_MSG("Done.");

        LOG_MSG("Loading sequences.");

        RUNP(sb = get_sequences_from_hdf5_model(param->in_model , IHMM_SEQ_READ_ONLY_SEQ  ));

        LOG_MSG("Read %d sequences.",sb->num_seq);

        LOG_MSG("Starting thread pool.");
       /* start threadpool  */
        //if((pool = thr_pool_create(param->num_threads , param->num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");


        /* allocate data for threads; */
        RUNP(td = create_seqer_thread_data(&param->num_threads,(sb->max_len+2)  , model_bag->max_num_states , NULL, THREAD_DATA_FULL));


        LOG_MSG("Done.");

        model_bag->best_model = -1;
        best_score = -100.0;
        max_likelihood = prob2scaledprob(1.0);

        limit = MACRO_MIN(max, sb->num_seq);

        total_seq_len = 0.0;
        for(i = 0; i < limit;i++){
                total_seq_len += sb->sequences[i]->seq_len;
        }
        //RUNP(all_scores = galloc(all_scores, limit,model_bag->num_models, 0.0));

        for(c = 0; c < model_bag->num_models;c++){
                fhmm = model_bag->finite_models[c];
                //RUN(alloc_dyn_matrices(fhmm));
                RUN(run_score_sequences(fhmm,sb,td));

                max_likelihood = prob2scaledprob(1.0);
                for(i = 0; i < limit;i++){
                        max_likelihood = logsum(max_likelihood,sb->sequences[i]->score);
                }

                RUN(calculate_BIC(fhmm, max_likelihood, total_seq_len, &BIC));
                LOG_MSG("Model: %d (%d states) ML: %.2f BIC %f (based on first %d seqs)",c,fhmm->K, max_likelihood,BIC,limit);
                if(BIC > best_score){
                        best_score = BIC;
                        model_bag->best_model = c;
                }
        }

        LOG_MSG("Best Model: %d",model_bag->best_model);
        RUN(write_best_model(param->in_model, model_bag->best_model));

        /*for(i = 0; i < limit;i++){
                s1 = 0.0;
                s2 = 0.0;

                for(c = 0; c < model_bag->num_models;c++){
                        s1 += all_scores[i][c];
                        s2 += all_scores[i][c] * all_scores[i][c];
                        fprintf(stdout,"%5.2f ",all_scores[i][c]);
                }

                s2 = sqrt(((double) model_bag->num_models * s2 - s1 * s1)/ ((double) model_bag->num_models * ((double) model_bag->num_models -1.0)));
                s1 = s1 / (double) model_bag->num_models;

                fprintf(stdout," %f stdev: %f \n",s1,s2);


        }


        //LOG_MSG("Mean KL divergence: %f stdev: %f (based on first %d seqs)",s1,s2,limit);
        gfree(all_scores);*/
        //LOG_MSG("Got past writing");
        free_seqer_thread_data(td);
        //LOG_MSG("Got past free thread data ");
        //thr_pool_destroy(pool);
        //LOG_MSG("Got past poolfree");
        free_ihmm_sequences(sb);
        //LOG_MSG("Got past seq free");
        free_model_bag(model_bag);
        //LOG_MSG("Got past model bag free");
        return OK;
ERROR:
        free_seqer_thread_data(td);
        //thr_pool_destroy(pool);
        free_ihmm_sequences(sb);
        free_fhmm(fhmm);
        return FAIL;
}


int free_parameters(struct parameters* param)
{
        ASSERT(param != NULL, " No param found - free'd already???");
        if(param->cmd_line){
                gfree(param->cmd_line);
        }
        free_rng(param->rng);
        /*if(param->tmp_file_A){
                MFREE(param->tmp_file_A);
        }
        if(param->tmp_file_B){
                MFREE(param->tmp_file_B);
                }*/
        MFREE(param);
        return OK;
ERROR:
        return FAIL;
}

int print_help(char **argv)
{
        const char usage[] = " -in <fasta> -out <outfile>";
        char* tmp = NULL;

        RUN(tlfilename(argv[0], &tmp));
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",tmp,usage);
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthreads","Number of threads." ,"[8]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--states","Number of starting states." ,"[10]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nmodels","Number of models to train." ,"[10]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--model","Continue training model <>." ,"[]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--niter","Number of iterations." ,"[1000]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--alpha","Alpha hyper parameter." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gamma","Gamma hyper oparameter." ,"[NA]"  );
        MFREE(tmp);
        return OK;
ERROR:
        MFREE(tmp);
        return FAIL;
}
