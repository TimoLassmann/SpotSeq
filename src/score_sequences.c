#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <libgen.h>


#include "tldevel.h"
#include "tlmisc.h"
#include "tllogsum.h"
#include "tlrng.h"
#include "tlseqio.h"
#include "tlalphabet.h"
#include "tlseqbuffer.h"
//#include "model.h"
#include "sequence_struct.h"



//#include "sequence_io.h"
//#include "sequence_prep.h"

#include "pst.h"
#include "pst_io.h"
#include "pst_search.h"


#include "thread_data.h"

//#include "ihmm_seq.h"
#include "bias_model.h"

#include "finite_hmm.h"
#include "finite_hmm_io.h"
#include "finite_hmm_alloc.h"

#include "finite_hmm_score.h"

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



static int scan_sequences_pst(struct parameters* param,struct tl_seq_buffer** hits,uint64_t* db_size);
static int print_help(char **argv);
static int free_parameters(struct parameters* param);

int main (int argc, char *argv[])
{
        FILE* fptr = NULL;

        struct parameters* param = NULL;
        struct fhmm** fhmm = NULL;
        //struct fhmm* bias = NULL;
        //struct tl_seq_buffer* sb = NULL;
        double* s = NULL;
        uint64_t db_size = 0;
        struct seqer_thread_data** td = NULL;


        int i;
        int c;

        //print_program_header(argv, "Scores sequences.");

        MMALLOC(param, sizeof(struct parameters));
        param->in_model = NULL;
        param->in_sequences = NULL;
        param->background_sequences = NULL;
        param->output = NULL;
        param->num_threads = 8;
        param->summary_file = NULL;
        param->threshold = 2.0;   /* z_score cutoff for pst model scores  */
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

        struct tl_seq_buffer* sb = NULL;
        RUN(scan_sequences_pst(param, &sb,&db_size ));

        LOG_MSG("Found %d putative hits", sb->num_seq);
        LOG_MSG("Read search fhmm");
        //sb->sequences[0]->score_arr

        MMALLOC(fhmm, sizeof(struct fhmm) * 2);
        RUN(read_searchfhmm(param->in_model, &fhmm[0]));

        RUN(read_biasfhmm(param->in_model, &fhmm[1]));
        //RUN(alloc_dyn_matrices(fhmm));
        RUN(create_seqer_thread_data(&td,param->num_threads,(sb->max_len+2)  , fhmm[0]->K+1, NULL));

        LOG_MSG("Run scoring");
        for(i = 0; i < sb->num_seq;i++){
                s = NULL;
                MMALLOC(s, sizeof(double)* 4);
                sb->sequences[i]->data = s;
        }
        RUN(run_score_sequences(fhmm,sb, td, 2, FHMM_SCORE_FULL ));
         /* Print scores.. */
        RUNP(fptr = fopen(param->output, "w"));
        fprintf(fptr, "Name,score,score_bias,p_score,p_score_bias,e,e_bias\n");
        
        for(i = 0; i < sb->num_seq;i++){
                s = sb->sequences[i]->data;
                sb->sequences[i]->name[strcspn(sb->sequences[i]->name, " ")] = 0;
                fprintf(fptr,"%s,%f,%f,%e,%e,%f,%f\n", sb->sequences[i]->name, s[0],s[1],s[2],s[3],s[2]* (double) db_size, s[3] * (double) db_size);
        }
        fclose(fptr);

        free_seqer_thread_data(td);
        for(i = 0; i < sb->num_seq;i++){
                 s = sb->sequences[i]->data;
                 MFREE(s);
                 sb->sequences[i]->data = NULL;
        }
        free_tl_seq_buffer(sb);
        free_fhmm(fhmm[0]);
        free_fhmm(fhmm[1]);
        MFREE(fhmm);
        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        if(fptr){
                fclose(fptr);
        }

        if(fhmm){
                free_fhmm(fhmm[0]);
                free_fhmm(fhmm[1]);
                MFREE(fhmm);
        }
        free_seqer_thread_data(td);
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


int scan_sequences_pst(struct parameters* param,struct tl_seq_buffer** hits,uint64_t* db_size)
{
        struct tl_seq_buffer* h = NULL;
        struct pst* p = NULL;

        init_logsum();
        LOG_MSG("Load PST model");
        RUN(read_pst_hdf5(&p, param->in_model));
        RUN(search_db(p, param->in_sequences, param->threshold,&h,db_size));
        free_pst(p);

        *hits = h;
        return OK;
ERROR:
        return FAIL;
}
