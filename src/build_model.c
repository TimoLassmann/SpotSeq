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

struct parameters{
        char* input;
        char* output;
        char* in_model;
        char* cmd_line;
        int num_iter;
        int local;
        int num_threads;
        int num_start_states;
};


static int run_build_ihmm(struct parameters* param);
static int print_help(char **argv);
static int free_parameters(struct parameters* param);

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        /*rk_state rndstate;
        float r;
        rk_randomseed(&rndstate);

        int i;
        for(i = 0; i < 10000;i++){

                r = rk_beta(&rndstate, 1.0, 11) * 0.0001;
                if(r == 0.0){
                        fprintf(stdout,"%d Got a zero: %f\n",i, r);
                }
                if(fabs(r-0.0) < FLT_EPSILON){
                        fprintf(stdout,"%d Got a zero: %f\n",i, r);
                }

        }
        exit(0);*/

        tlog.echo_build_config();

        MMALLOC(param, sizeof(struct parameters));
        param->input = NULL;
        param->output = NULL;
        param->in_model = NULL;
        param->cmd_line = NULL;
        param->num_threads = 8;
        param->num_start_states = 10;
        param->local = 0;
        param->num_iter = 1000;

        while (1){
                static struct option long_options[] ={
                        {"in",required_argument,0,'i'},
                        {"out",required_argument,0,'o'},
                        {"states",required_argument,0,'s'},
                        {"local",no_argument,0,'l'},
                        {"nthreads",required_argument,0,'t'},
                        {"niter",required_argument,0,'n'},
                        {"model",required_argument,0,'m'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"hi:o:t:n:m:s:l",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
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

        RUN(run_build_fhmm_file(param->output));
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

        int initial_states = 400;
        int i;

        ASSERT(param!= NULL, "No parameters found.");

        initial_states = param->num_start_states;

        if(param->in_model){
                RUNP(model = read_model_hdf5(param->in_model));
                // print_model_parameters(model);
                // print_counts(model);
                RUNP(sb = get_sequences_from_hdf5_model(param->in_model));

                //RUN(fill_counts(model, sb));

        }else{
                /* Step one read in sequences */
                LOG_MSG("Loading sequences.");


                RUNP(sb = load_sequences(param->input));

                RUNP(model = alloc_ihmm_model(initial_states, sb->L));
                /* Initial guess... */
                model->alpha0_a = 6.0f;
                model->alpha0_b = 15.0f;
                model->gamma_a = 16.0f;
                model->gamma_b = 4.0f;
                /* vague */
                /*model->alpha0_a = 1.0f;
                model->alpha0_b = 1.0f;
                model->gamma_a = 2.0f;
                model->gamma_b = 1.0f;*/

                model->alpha = IHMM_PARAM_PLACEHOLDER;
                model->gamma = IHMM_PARAM_PLACEHOLDER;
                RUN(inititalize_model(model, sb,initial_states));// initial_states) );

                for(i = 0;i < 10;i++){
                        RUN(iHmmHyperSample(model, 20));
                }
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
        return OK;
}

