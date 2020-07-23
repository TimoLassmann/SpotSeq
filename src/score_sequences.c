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
#include "tlseqio.h"
#include "tlalphabet.h"

//#include "model.h"
#include "sequence_struct.h"

#include "sequence_io.h"
#include "sequence_prep.h"

#include "pst.h"
#include "pst_io.h"
#include "pst_search.h"


#include "thread_data.h"

//#include "ihmm_seq.h"

#include "finite_hmm.h"
#include "finite_hmm_io.h"
#include "finite_hmm_alloc.h"

#include "run_score.h"

struct parameters{
        char* in_model;
        char* in_sequences;
        char* background_sequences;
        char* output;
        char* summary_file;
        double threshold;
        int num_threads;
        rk_state rndstate;
        struct rng_state* rng;

};

static int scan_sequences_pst(struct parameters* param,struct tl_seq_buffer** hits);

static int print_help(char **argv);
static int free_parameters(struct parameters* param);

int main (int argc, char *argv[])
{
        FILE* fptr = NULL;

        struct parameters* param = NULL;
        struct fhmm* fhmm = NULL;
        struct seq_buffer* sb = NULL;


        struct seqer_thread_data** td = NULL;


        int i,c;

        //print_program_header(argv, "Scores sequences.");

        MMALLOC(param, sizeof(struct parameters));
        param->in_model = NULL;
        param->in_sequences = NULL;
        param->background_sequences = NULL;
        param->output = NULL;
        param->num_threads = 8;
        param->summary_file = NULL;
        param->threshold = 5.0;   /* z_score cutoff for pst model scores  */
        param->rng = NULL;

        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"in",required_argument,0,'i'},
                        {"out",required_argument,0,'o'},
                        {"nthreads",required_argument,0,'t'},
                        {"background",required_argument,0,'b'},
                        {"summary",required_argument,0,'s'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"hm:i:s:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 'i':
                        param->in_sequences = optarg;
                        break;
                case 'b':
                        param->background_sequences = optarg;
                        break;
                case 'o':
                        param->output = optarg;
                        break;
                case 't':
                        param->num_threads = atoi(optarg);
                        break;
                case 'm':
                        param->in_model = optarg;
                        break;
                case 's':
                        param->summary_file = optarg;
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

        if(42){
                rk_seed(42, &param->rndstate);
                RUNP(param->rng = init_rng(42));
        }else{
                rk_randomseed(&param->rndstate);
        }
        LOG_MSG("Starting run");

        if(!param->in_sequences){
                RUN(print_help(argv));
                ERROR_MSG("No input sequences! use -i <blah.fa>");

        }else{
                if(!my_file_exists(param->in_sequences)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->in_sequences);
                }
        }

        if(!param->in_model){
                RUN(print_help(argv));
                ERROR_MSG("No model file! use -m  <blah.h5>");
        }else{
                if(!my_file_exists(param->in_model)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->in_model);
                }
        }

        if(param->background_sequences ){
                if(!my_file_exists(param->background_sequences)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->background_sequences);
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

        struct tl_seq_buffer* hits = NULL;
        RUN(scan_sequences_pst(param, &hits));

        RUN(convert_tl_seq_buf_into_ihmm_seq_buf(hits, &sb));
        free_tl_seq_buffer(hits);

        LOG_MSG("Read search fhmm");
        //sb->sequences[0]->score_arr
        RUN(read_searchfhmm(param->in_model, &fhmm));
        //RUN(alloc_dyn_matrices(fhmm));
        RUN(create_seqer_thread_data(&td,param->num_threads,(sb->max_len+2)  , fhmm->K+1, NULL));
        LOG_MSG("Run scoring");
        RUN(run_score_sequences(fhmm,sb, td));
         /* Print scores.. */
        //RUNP(fptr = fopen(param->output, "w"));
        //fprintf(fptr, "Name,Score_%s\n",  param->in_model);
        for(i = 0; i < sb->num_seq;i++){
                fprintf(stdout, "%f %f %s\n",sb->sequences[i]->score ,esl_exp_logsurv(sb->sequences[i]->score, fhmm->tau, fhmm->lambda)  ,   sb->sequences[i]->name);
        }
        //fclose(fptr);

        free_seqer_thread_data(td);
        //thr_pool_destroy(pool);
        //free_ihmm_sequences(sb);
        free_fhmm(fhmm);
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
        const char usage[] = " -m <model.h5> -i <input sequences> ";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthreads","Number of threads." ,"[8]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--background","Background sequences - residue counts from these will be ADDED to the background model. " ,"[8]"  );
        return OK;
}


int scan_sequences_pst(struct parameters* param,struct tl_seq_buffer** hits)
{
        struct tl_seq_buffer* h = NULL;
        struct pst* p = NULL;
        int i;
        init_logsum();
        LOG_MSG("Load PST model");
        RUN(read_pst_hdf5(&p, param->in_model));
        RUN(search_db(p, param->in_sequences, param->threshold,&h));
        /*for(i = 0; i < h->num_seq;i++){
                fprintf(stdout,"%d: %s\n", i, h->sequences[i]->name);
                }*/
        free_pst(p);

        //free_tl_seq_buffer(hits);

        *hits = h;
        return OK;
ERROR:
        return FAIL;
}
