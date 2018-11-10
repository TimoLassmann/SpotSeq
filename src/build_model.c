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

struct parameters{
        char* input;
        char* output;
        char* in_model;
        char* cmd_line;
        float alpha;
        float gamma;
        unsigned long seed;
        int num_iter;
        int local;
        int num_threads;
        int num_start_states;
        int rev;
};

static int run_build_ihmm(struct parameters* param);

static int random_score_sequences(struct seq_buffer* sb,float* background );
static int score_sequences_for_command_line_reporting(struct parameters* param);

static int print_help(char **argv);
static int free_parameters(struct parameters* param);

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        tlog.echo_build_config();

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
        while (1){
                static struct option long_options[] ={
                        {"in",required_argument,0,'i'},
                        {"out",required_argument,0,'o'},
                        {"states",required_argument,0,'s'},
                        {"local",no_argument,0,'l'},
                        {"nthreads",required_argument,0,'t'},
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

        RUNP(param->cmd_line = make_cmd_line(argc,argv));

        RUN(run_build_ihmm(param));
        /* 1 means allow transitions that are not seen in the training
         * data */
        RUN(run_build_fhmm_file(param->output,0));

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
        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;
        struct seq_buffer* sb = NULL;

        int initial_states;
        int i;

        ASSERT(param!= NULL, "No parameters found.");
        init_logsum();
        initial_states = param->num_start_states;

        if(param->in_model){
                RUNP(model = read_model_hdf5(param->in_model));
                // print_model_parameters(model);
                // print_counts(model);
                RUNP(sb = get_sequences_from_hdf5_model(param->in_model));

                /*model->alpha_a = 6.0f;
                model->alpha_b = 15.0f;
                model->alpha = rk_gamma(&model->rndstate, model->alpha_a,1.0 / model->alpha_b);
                model->gamma_a = 16.0f;
                model->gamma_b = 4.0f;
                model->gamma = rk_gamma(&model->rndstate, model->gamma_a,1.0 / model->gamma_b);*/
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


        }else{
                /* Step one read in sequences */
                LOG_MSG("Loading sequences.");


                RUNP(sb = load_sequences(param->input));
                //sb = concatenate_sequences(sb);
                //RUN(label_seq_based_on_random_fhmm(sb, initial_states, 0.3));



                if(param->rev && sb->L == ALPHABET_DNA){
                        LOG_MSG("Add revcomp sequences.");
                        RUN(add_reverse_complement_sequences_to_buffer(sb));
                }
                RUNP(model = alloc_ihmm_model(initial_states+2, sb->L, param->seed));




                if(param->alpha == IHMM_PARAM_PLACEHOLDER){
                        model->alpha_a = 6.0f;
                        model->alpha_b = 15.0f;
                        model->alpha = rk_gamma(&model->rndstate, model->alpha_a,1.0 / model->alpha_b);
                }else{
                        model->alpha = param->alpha;
                }

                if(param->gamma == IHMM_PARAM_PLACEHOLDER){
                        model->gamma_a = 16.0f;
                        model->gamma_b = 4.0f;
                        model->gamma = rk_gamma(&model->rndstate, model->gamma_a,1.0 / model->gamma_b);
                }else{
                        model->gamma = param->gamma;
                }
                model->gamma_limit = 50.0f;
                model->alpha_limit = 2.0f;

                //RUN(fill_counts(model,sb));
                /* I am doing this as a pre-caution. I don't want the inital model
                 * contain states that are not visited.. */
                //RUN(remove_unused_states_labels(model, sb));
                //RUN(fill_counts(model,sb));
                //RUN(inititalize_model(model, sb,initial_states));// initial_states) );

                /* Note this also initializes the last (to infinity state) */

                for(i = 0; i < model->num_states;i++){
                        model->beta[i] = (float)(model->num_states);
                }

                for(i = 0;i < 10;i++){
                        RUN(iHmmHyperSample(model, 20));
                }
        }
        /* Set seed in sequence buffer */
        if(param->seed){
                rk_seed(param->seed + 2, &sb->rndstate);
        }else{
                rk_randomseed(&sb->rndstate);
        }

        /*LOG_MSG("testing");
        int a,b,r;
        float test[4];
        float sum;
        float min, max;

        for(a = 1;a < 100;a++){

                for(r = 0; r < 10;r++){
                        sum = 0.0f;
                        for(b = 0; b < 4;b++){
                                test[b] = rk_gamma(&model->rndstate,10.0 + (float) a, 1.0);
                                sum += test[b];
                        }
                        min = 2.0;
                        max = -2.0;
                        fprintf(stdout,"%d\t",a);
                        for(b = 0; b < 4;b++){
                                test[b] = test[b] / sum;
                                fprintf(stdout," %0.2f",test[b]);
                                if(max <  test[b]){
                                        max= test[b];
                                }
                                if(min > test[b]){
                                        min= test[b];
                                }
                        }
                        fprintf(stdout,"\t%f-%f range:%f\n",min,max,max-min);
                }
        }

        exit(0);*/
        LOG_MSG("Read %d sequences.",sb->num_seq);

        RUNP(ft = alloc_fast_hmm_param(initial_states,sb->L));
        RUN(fill_background_emission(ft, sb));

        /* Do a random score */

        RUN(random_score_sequences(sb, ft->background_emission  ));

        RUN(run_beam_sampling( model, sb, ft,NULL,  param->num_iter,  param->num_threads));

        //RUN(write_model(model, param->output));
        RUN(write_model_hdf5(model, param->output));
//        RUN(add_annotation)
        RUN(add_annotation(param->output, "spotseq_model_cmd", param->cmd_line));
        //RUN(add_background_emission(param->output,ft->background_emission,ft->L));
        RUN(add_sequences_to_hdf5_model(param->output, sb));
        //RUN(print_states_per_sequence(sb));
        //RUN(write_ihmm_sequences(sb,"test.lfs","testing"));
        //sb, num thread, guess for aplha and gamma.. iterations.

        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        free_ihmm_sequences(sb);
        return OK;
ERROR:
        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        free_ihmm_sequences(sb);
        return FAIL;
}

int random_score_sequences(struct seq_buffer* sb,float* background )
{
        struct ihmm_sequence* s;
        float* back = NULL;
        int expected_len;
        int i;

        ASSERT(sb!=NULL, "No sequences");
        ASSERT(background != NULL, "No background");


        MMALLOC(back, sizeof(float) * sb->L);
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
        return OK;
ERROR:
        return FAIL;
}


/* After training a model it would be nice to know the logodds scores
 * of a number of sequences to see if something was found */

int score_sequences_for_command_line_reporting(struct parameters* param)
{
        struct seq_buffer* sb = NULL;
        struct fhmm* fhmm;
        double s1,s2;
        int max = 100;
        int i;
        int limit;
        /* Maybe runscoring on top X sequences to report if something was found... */
        LOG_MSG("Loading model.");
        RUNP(fhmm = init_fhmm(param->output));

        /* load sequences.  */
        LOG_MSG("Loading sequences.");
        RUNP(sb =get_sequences_from_hdf5_model(param->output));

        LOG_MSG("Read %d sequences.",sb->num_seq);

        RUN(run_score_sequences(fhmm,sb, param->num_threads));

        limit = MACRO_MIN(max, sb->num_seq);

        s1 = 0.0;
        s2 = 0.0;
        for(i = 0; i < limit;i++){
                //sb_in->sequences[i]->seq_len = 10 + (int)(rk_double(&rndstate)*10.0) - 5.0;
                s1 += sb->sequences[i]->score;
                s2 += (sb->sequences[i]->score * sb->sequences[i]->score);
        }

        s2 = sqrt(((double) limit * s2 - s1 * s1)/ ((double) limit * ((double) limit -1.0)));
        s1 = s1 / (double) limit;



        LOG_MSG("Mean log-odds ratio: %f stdev: %f (based on first %d seqs)",s1,s2,limit);

        free_ihmm_sequences(sb);
        free_fhmm(fhmm);

        return OK;
ERROR:
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
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--model","Continue training model <>." ,"[]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--niter","Number of iterations." ,"[1000]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--alpha","Alpha hyper parameter." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gamma","Gamma hyper oparameter." ,"[NA]"  );
        return OK;
}
