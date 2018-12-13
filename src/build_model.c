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

#include "ihmm_seq.h"

#include "beam_sample.h"

#include "finite_hmm.h"

#include "distributions.h"

#include "run_score.h"

#include "init_seq_label.h"


#define OPT_SEED 1
#define OPT_NUM_MODELS 2

struct parameters{
        char* input;
        char* output;
        char* in_model;
        char* cmd_line;
        double alpha;
        double gamma;
        unsigned long seed;
        int num_iter;
        int local;
        int num_models;
        int num_threads;
        int num_start_states;
        int rev;
};

static int run_build_ihmm(struct parameters* param);

static int random_score_sequences(struct seq_buffer* sb,double* background );
static int score_sequences_for_command_line_reporting(struct parameters* param);

static int print_help(char **argv);
static int free_parameters(struct parameters* param);

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        print_program_header(argv, "Build HDPHMM model(s).");

        MMALLOC(param, sizeof(struct parameters));
        param->input = NULL;
        param->output = NULL;
        param->in_model = NULL;
        param->cmd_line = NULL;
        param->num_threads = 8;
        param->num_start_states = 10;
        param->local = 0;
        param->rev = 0;
        param->num_iter = 1000;
        param->alpha = IHMM_PARAM_PLACEHOLDER;
        param->gamma = IHMM_PARAM_PLACEHOLDER;
        param->seed = 0;
        param->num_models = 10;
        while (1){
                static struct option long_options[] ={
                        {"in",required_argument,0,'i'},
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
                        param->output = optarg;
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

        if(!param->output){
                RUN(print_help(argv));
                ERROR_MSG("No output file! use --in <blah.fa>");
        }else{
                if(my_file_exists(param->output)){
                        WARNING_MSG("Will overwrite: %s.",param->output);
                }
        }

        if(param->in_model){
                if(!my_file_exists(param->in_model)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->in_model);
                }
        }

        if(param->num_models < 1){
                RUN(print_help(argv));
                ERROR_MSG("To few models! use -nmodels <1+>");
        }

        RUNP(param->cmd_line = make_cmd_line(argc,argv));

        //rk_save_testing();

        //return EXIT_SUCCESS;
        RUN(run_build_ihmm(param));
        /* 1 means allow transitions that are not seen in the training
         * data */
        //

        RUN(score_sequences_for_command_line_reporting(param));

        /* calibrate model parameters */
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

        struct spotseq_thread_data** td = NULL;

        struct thr_pool* pool = NULL;
        int* num_state_array = NULL;

        //int initial_states;
        int i,j;

        ASSERT(param!= NULL, "No parameters found.");
        init_logsum();


        MMALLOC(num_state_array, sizeof(int)* param->num_models);
        if(param->in_model){
                /* PROBABLY need to re-alloc num_state_array */
                //MMALLOC(num_state_array, sizeof(int)* param->num_models);
                RUNP(model_bag = read_model_bag_hdf5(param->in_model));
                RUNP(sb = get_sequences_from_hdf5_model(param->in_model,IHMM_SEQ_READ_ALL));

                for(i = 0; i < model_bag->num_models;i++){
                        num_state_array[i] = model_bag->models[i]->num_states;
                }



                RUNP(td = read_thread_data_to_hdf5(param->in_model));
                //MFREE(num_state_array);
                //free_ihmm_sequences(sb);
                //free_model_bag(model_bag);
                //return EXIT_SUCCESS;
                /*
                if(param->alpha != IHMM_PARAM_PLACEHOLDER){
                        model->alpha = param->alpha;
                        model->alpha_a = IHMM_PARAM_PLACEHOLDER;
                        model->alpha_b = IHMM_PARAM_PLACEHOLDER;
                }
                if(param->gamma != IHMM_PARAM_PLACEHOLDER){
                        model->gamma = param->gamma;
                        model->gamma_a = IHMM_PARAM_PLACEHOLDER;
                        model->gamma_b = IHMM_PARAM_PLACEHOLDER;
                }
                */
                /*RUN(write_model_bag_hdf5(model_bag,param->output));
                RUN(add_annotation(param->output, "spotseq_model_cmd", param->cmd_line));
                RUN(add_background_emission(param->output,  sb->background, sb->L));

                RUN(add_sequences_to_hdf5_model(param->output, sb,  model_bag->num_models));
                RUN(write_thread_data_to_hdf5(param->output, td, param->num_threads, sb->max_len,  model_bag->max_num_states));
                return EXIT_SUCCESS;*/
        }else{
                /* Step one read in sequences */
                LOG_MSG("Loading sequences.");

                RUNP(sb = load_sequences(param->input));

                //sb = concatenate_sequences(sb);
                //RUN(label_seq_based_on_random_fhmm(sb, initial_states, 0.3));
                LOG_MSG("%d models ",param->num_models);

                if(param->rev && sb->L == ALPHABET_DNA){
                        LOG_MSG("Add revcomp sequences.");
                        RUN(add_reverse_complement_sequences_to_buffer(sb));
                }

                num_state_array[0] = 10;
                for(i = 1; i < param->num_models;i++){
                        num_state_array[i] = num_state_array[i-1] + 10;//( (sb->max_len / 2) / param->num_models);
                }
                RUNP(model_bag = alloc_model_bag(num_state_array, sb->L, param->num_models,   param->seed));


                RUN(add_multi_model_label_and_u(sb, model_bag->num_models));
                //label_ihmm_sequences_based_on_guess_hmm(struct seq_buffer *sb, int k, float alpha)
                /* New label sequences  */
                /* We need to label sequences somehow so that we can
                 * set the u arrays properly based on the initial
                 * model guess */
                model_bag->max_num_states = 0;
                for(i = 0; i < model_bag->num_models;i++){
                        RUN(random_label_based_on_multiple_models(sb, model_bag->models[i]->num_states,i,&model_bag->rndstate));

                        RUN(fill_counts(model_bag->models[i], sb,i));


                        model_bag->max_num_states = MACRO_MAX(model_bag->max_num_states,model_bag->models[i]->num_states);
                        //print_counts(model_bag->models[i]);


                }

                RUN(check_labels(sb,model_bag->num_models ));
                /* Set seed in sequence buffer */
                sb->seed = rk_ulong(&model_bag->rndstate);
                rk_seed(sb->seed, &sb->rndstate);

                RUN(set_model_hyper_parameters(model_bag, param->alpha, param->gamma));

                /* Allocating thread structure. */
                RUNP(td = create_spotseq_thread_data(&param->num_threads,(sb->max_len+2)  ,model_bag->max_num_states, &model_bag->rndstate));


        }

        LOG_MSG("Will use %d threads.", param->num_threads);
        if((pool = thr_pool_create(param->num_threads,param->num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");

        RUNP(ft_bag = alloc_fast_param_bag(model_bag->num_models, num_state_array, sb->L));
        //RUNP(ft = alloc_fast_hmm_param(initial_states,sb->L));
        /* fill background of first fast hmm param struct  */
        RUN(fill_background_emission(ft_bag->fast_params[0], sb));
        /* Now copy remaining 1:N first fast hmm param structs */
        for(i = 1; i < ft_bag->num_models;i++){
                for(j = 0; j < sb->L;j++){
                        ft_bag->fast_params[i]->background_emission[j] = ft_bag->fast_params[0]->background_emission[j];
                }

        }
        /*for(i = 0; i < 1000000;i++){

                long a = rk_random(&model_bag->rndstate) >> 9, b = rk_random(&model_bag->rndstate) >> 4;
                float test = a / 8388609.0f;
                if(test < 0.00001 || test > 0.9999){
                        fprintf(stdout,"%ld %ld  %f\n ",a,b, test);
                }

        }
        exit(0);*/

        /* Do a random score */

        RUN(random_score_sequences(sb, ft_bag->fast_params[0]->background_emission  ));

        /* Main function */
        RUN(run_beam_sampling(model_bag,ft_bag, sb,td, pool,  param->num_iter,  param->num_threads));


        /* Write results */
        RUN(convert_ihmm_to_fhmm_models(model_bag));
        RUN(write_model_bag_hdf5(model_bag,param->output));
        RUN(add_annotation(param->output, "spotseq_model_cmd", param->cmd_line));
        RUN(add_sequences_to_hdf5_model(param->output, sb,  model_bag->num_models));
        RUN(write_thread_data_to_hdf5(param->output, td, param->num_threads, sb->max_len, model_bag->max_num_states));
        //RUN(write_thread_data_to_)
        //RUN(write_model(model, param->output));
        /*for(i = 0; i < model_bag->num_models;i++){
                char buffer[BUFFER_LEN];
                snprintf(buffer, BUFFER_LEN, "%s_%d.h5", param->output,i);
                RUN(write_model_hdf5(model_bag->models[i],buffer));
                RUN(add_annotation(buffer, "spotseq_model_cmd", param->cmd_line));
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
        free_spotseq_thread_data(td, param->num_threads);
        thr_pool_destroy(pool);
        MFREE(num_state_array);
        return OK;
ERROR:
        free_ihmm_sequences(sb);
        free_model_bag(model_bag);
        free_fast_param_bag(ft_bag);
        free_spotseq_thread_data(td, param->num_threads);
        thr_pool_destroy(pool);
        MFREE(num_state_array);
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
                RUN(random_model_score(back, &s->r_score , s->seq, s->seq_len,expected_len));
        }
        MFREE(back);
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

        struct spotseq_thread_data** td = NULL;

        struct thr_pool* pool = NULL;

        struct seq_buffer* sb = NULL;
        struct fhmm* fhmm = NULL;

        double s1,s2;
        double best_score;
        int max = 1000;
        int c;
        int i;
        int limit;

        //char buffer[BUFFER_LEN];



        LOG_MSG("Reading models.");
        RUNP(model_bag = read_model_bag_hdf5(param->output));
        LOG_MSG("Done.");

        LOG_MSG("Loading sequences.");

        RUNP(sb =get_sequences_from_hdf5_model(param->output, IHMM_SEQ_READ_ONLY_SEQ  ));

        LOG_MSG("Read %d sequences.",sb->num_seq);

        LOG_MSG("Starting thread pool.");
       /* start threadpool  */
        if((pool = thr_pool_create(param->num_threads , param->num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");


        /* allocate data for threads; */
        RUNP(td = create_spotseq_thread_data(&param->num_threads,(sb->max_len+2)  , model_bag->max_num_states , NULL));


        LOG_MSG("Done.");

        model_bag->best_model = -1;
        best_score = -100.0;
        for(c = 0; c < model_bag->num_models;c++){
                fhmm = model_bag->finite_models[c];
                RUN(alloc_dyn_matrices(fhmm));
                RUN(run_score_sequences(fhmm,sb,td, pool));

                limit = MACRO_MIN(max, sb->num_seq);

                s1 = 0.0;
                s2 = 0.0;
                for(i = 0; i < limit;i++){
                        //sb_in->sequences[i]->senq_len = 10 + (int)(rk_double(&rndstate)*10.0) - 5.0;
                        s1 += sb->sequences[i]->score;
                        s2 += (sb->sequences[i]->score * sb->sequences[i]->score);
                }

                s2 = sqrt(((double) limit * s2 - s1 * s1)/ ((double) limit * ((double) limit -1.0)));
                s1 = s1 / (double) limit;



                LOG_MSG("Model: %d Mean log-odds ratio: %f stdev: %f (based on first %d seqs)",c,s1,s2,limit);
                if(s1 > best_score){
                        best_score = s1;
                        model_bag->best_model = c;
                }
        }

        LOG_MSG("Best Model: %d",model_bag->best_model);
        RUN(write_best_model(param->output, model_bag->best_model));
        //LOG_MSG("Got past writing");
        free_spotseq_thread_data(td,param->num_threads);
        //LOG_MSG("Got past free thread data ");
        thr_pool_destroy(pool);
        //LOG_MSG("Got past poolfree");
        free_ihmm_sequences(sb);
        //LOG_MSG("Got past seq free");
        free_model_bag(model_bag);
        //LOG_MSG("Got past model bag free");


        return OK;
ERROR:
        free_spotseq_thread_data(td,param->num_threads);
        thr_pool_destroy(pool);
        free_ihmm_sequences(sb);
        free_fhmm(fhmm);
        return FAIL;
}


int free_parameters(struct parameters* param)
{
        ASSERT(param != NULL, " No param found - free'd already???");
        if(param->cmd_line){
                MFREE(param->cmd_line);
        }
        MFREE(param);
        return OK;
ERROR:
        return FAIL;
}

int print_help(char **argv)
{
        const char usage[] = " -in <fasta> -out <outfile>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthreads","Number of threads." ,"[8]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--states","Number of starting states." ,"[10]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nmodels","Number of models to train." ,"[10]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--model","Continue training model <>." ,"[]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--niter","Number of iterations." ,"[1000]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--alpha","Alpha hyper parameter." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gamma","Gamma hyper oparameter." ,"[NA]"  );
        return OK;
}
