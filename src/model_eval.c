#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <libgen.h>

#include "tldevel.h"
#include "tlmisc.h"
#include "tllogsum.h"
#include "tlrng.h"

#include "tlalphabet.h"

//#include "model.h"
#include "sequence_struct.h"
#include "sequence_alloc.h"
#include "sequence_io.h"
#include "sequence_prep.h"
#include "sequences_sim.h"

#include "thread_data.h"

#include "null_model_emission.h"
//#include "ihmm_seq.h"

#include "finite_hmm.h"
#include "finite_hmm_io.h"
#include "finite_hmm_alloc.h"

#include "run_score.h"

struct parameters{
        char* in_model;
        char* output;
        int num_threads;
        struct rng_state* rng;
};



static int print_help(char **argv);
static int free_parameters(struct parameters* param);

int main (int argc, char *argv[])
{
        FILE* fptr = NULL;

        struct parameters* param = NULL;
        struct fhmm* fhmm = NULL;
        struct seq_buffer* sb = NULL;


        struct seqer_thread_data** td = NULL;

        double p;
        int i,c;
        int num_test_seq = 10000;
        //print_program_header(argv, "Scores sequences.");

        MMALLOC(param, sizeof(struct parameters));
        param->in_model = NULL;
        param->output = NULL;
        param->num_threads = 8;
        param->rng = NULL;

        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"out",required_argument,0,'o'},
                        {"nthreads",required_argument,0,'t'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"o:t:m:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 'o':
                        param->output = optarg;
                        break;
                case 't':
                        param->num_threads = atoi(optarg);
                        break;
                case 'm':
                        param->in_model = optarg;
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

        if(1){
                //rk_seed(42, &param->rndstate);
                RUNP(param->rng = init_rng(0));
        }else{
                //rk_randomseed(&param->rndstate);
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

        if(!param->output){
                RUN(print_help(argv));
                ERROR_MSG("No output file! use -o   <blah.csv>");
        }else{
                if(my_file_exists(param->output)){
                       WARNING_MSG("The file %s will be over-written.",param->output);
                }
        }


        LOG_MSG("Loading model.");


        init_logsum();

        /* struct tl_seq_buffer* hits = NULL; */


        /* RUN(convert_tl_seq_buf_into_ihmm_seq_buf(hits, &sb)); */
        /* free_tl_seq_buffer(hits); */

        /* LOG_MSG("Found %d putative hits", sb->num_seq); */

        LOG_MSG("Read search fhmm");
        //sb->sequences[0]->score_arr
        RUN(read_searchfhmm(param->in_model, &fhmm));

        RUN(sim_sequences(num_test_seq, fhmm->L, 100,&sb, param->rng));
        RUN(create_seqer_thread_data(&td,param->num_threads,(sb->max_len+2)  , fhmm->K+1, NULL));
        RUN(run_score_sequences(fhmm,sb, td));
        int g = 0;
        for(i = 0; i < sb->num_seq;i++){
                p = esl_exp_surv(sb->sequences[i]->score, fhmm->tau, fhmm->lambda);
                if(p < 0.02){
                        fprintf(stdout, "%f %f %s\n",sb->sequences[i]->score ,p ,   sb->sequences[i]->name);
                        g++;
                }


        }
        LOG_MSG("%d %f ", g , (double)g / (double)sb->num_seq);
        RUN(sim_sequences(num_test_seq, fhmm->L, 400,&sb, param->rng));
        RUN(resize_seqer_thread_data(td  ,(sb->max_len+2)  , fhmm->K+1));
        RUN(run_score_sequences(fhmm,sb, td));
        g = 0;
        for(i = 0; i < sb->num_seq;i++){
                p = esl_exp_surv(sb->sequences[i]->score, fhmm->tau, fhmm->lambda);
                if(p < 0.02){
                        fprintf(stdout, "%f %f %s\n",sb->sequences[i]->score ,p ,   sb->sequences[i]->name);
                        g++;
                }



        }
        LOG_MSG("%d %f ", g , (double)g / (double)sb->num_seq);

        RUN(sim_sequences(num_test_seq , fhmm->L, 1600,&sb, param->rng));
        RUN(resize_seqer_thread_data(td  ,(sb->max_len+2)  , fhmm->K+1));
        RUN(run_score_sequences(fhmm,sb, td));
        g = 0;
        for(i = 0; i < sb->num_seq;i++){
                p = esl_exp_surv(sb->sequences[i]->score, fhmm->tau, fhmm->lambda);
                if(p < 0.02){
                        fprintf(stdout, "%f %f %s\n",sb->sequences[i]->score ,p ,   sb->sequences[i]->name);
                        g++;
                }


        }
        LOG_MSG("%d %f ", g , (double)g / (double)sb->num_seq);

        free_ihmm_sequences(sb);
        free_seqer_thread_data(td);
        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        if(fptr){
                fclose(fptr);
        }

        free_seqer_thread_data(td);
        //
        //free_ihmm_sequences(sb);
        free_fhmm(fhmm);
        free_parameters(param);
        return EXIT_FAILURE;
}


int free_parameters(struct parameters* param)
{
        ASSERT(param != NULL, " No param found - free'd already???");
        if(param->rng){
                free_rng(param->rng);
        }
        MFREE(param);
        return OK;
ERROR:
        return FAIL;
}

int print_help(char **argv)
{
        const char usage[] = " -m <model.h5> -o <outfile> ";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthreads","Number of threads." ,"[8]"  );

        return OK;
}
