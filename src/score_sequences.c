#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include "tldevel.h"
#include "model.h"
#include "ihmm_seq.h"

#include "finite_hmm.h"

#include "run_score.h"

struct parameters{
        char* in_model;
        char* in_sequences;
        char* background_sequences;
        char* output;
        char* summary_file;
        int num_threads;
        rk_state rndstate;
};


static int print_help(char **argv);
static int free_parameters(struct parameters* param);

int main (int argc, char *argv[])
{
        FILE* fptr = NULL;

        struct parameters* param = NULL;
        struct fhmm* fhmm = NULL;
        struct seq_buffer* sb = NULL;
        struct seq_buffer* sb_back = NULL;

        struct wims_thread_data** td = NULL;

        struct thr_pool* pool = NULL;


        int i,c;

        print_program_header(argv, "Scores sequences.");

        MMALLOC(param, sizeof(struct parameters));
        param->in_model = NULL;
        param->in_sequences = NULL;
        param->background_sequences = NULL;
        param->output = NULL;
        param->num_threads = 8;
        param->summary_file = NULL;

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

        int best = -1;
        init_logsum();
        RUNP(fhmm = read_best_fmodel(param->in_model, &best));
        RUN(alloc_dyn_matrices(fhmm));
        /* load sequences.  */
        LOG_MSG("Loading sequences.");
        RUNP(sb = load_sequences(param->in_sequences, &param->rndstate));
        LOG_MSG("Read %d sequences.",sb->num_seq);
        /* we need to use the background residue distribution in the sequences to test for our random model! Somehow I did not do this before... */

        for(i =0; i < sb->L;i++){
                fhmm->background[i] = 0.0;
        }

        RUN(get_res_counts(sb, fhmm->background));
        if(param->background_sequences){
                LOG_MSG("Loading background sequences.");
                RUNP(sb_back = load_sequences(param->background_sequences, &param->rndstate));
                LOG_MSG("Read %d sequences.",sb->num_seq);
                RUN(get_res_counts(sb_back, fhmm->background));
        }

        double sum = 0;
        for(i =0; i < sb->L;i++){
                sum += fhmm->background[i];
                //fprintf(stdout,"%d: %f\n", i,scaledprob2prob( fhmm->background[i]));
        }

        //for(i =0; i < sb->L;i++){
        //        fhmm->background[i]/=sum;
        //        fprintf(stdout,"%d: %f\n", i,fhmm->background[i]);
        //}


        for(i =0; i < sb->L;i++){
                fhmm->background[i] = prob2scaledprob(fhmm->background[i] / sum);
        }

       /* start threadpool  */
        if((pool = thr_pool_create(param->num_threads , param->num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");

        /* allocate data for threads; */
        RUNP(td = create_wims_thread_data(&param->num_threads,(sb->max_len+2)  , fhmm->K+1, NULL));

        RUN(run_score_sequences(fhmm,sb, td, pool));
         /* Print scores.. */
        RUNP(fptr = fopen(param->output, "w"));
        fprintf(fptr, "Name,Score_%s\n",  param->in_model);
        for(i = 0; i < sb->num_seq;i++){
                fprintf(fptr, "%s,%f\n",sb->sequences[i]->name, sb->sequences[i]->score);// /  (1.0 + ex
        }
        fclose(fptr);

        if(param->summary_file){
                double s1,s2;
                s1 = 0.0;
                s2 = 0.0;

                for(i = 0; i < sb->num_seq;i++){
                        s1 += sb->sequences[i]->score;
                        s2 += sb->sequences[i]->score * sb->sequences[i]->score;
                }

                s2 = sqrt(((double) sb->num_seq * s2 - s1 * s1)/ ((double) sb->num_seq * ((double)sb->num_seq -1.0)));
                s1 = s1 / (double) sb->num_seq;

                fptr = NULL;
                if(!my_file_exists(param->summary_file)){
                        RUNP(fptr = fopen(param->summary_file, "w"));
                        fprintf(fptr,"Model,SequenceFile,Mean,Stdev\n");
                }else{
                        RUNP(fptr = fopen(param->summary_file, "a"));
                }
                fprintf(fptr,"%s,%s,%f,%f\n",basename(param->in_model), basename(param->in_sequences),s1,s2);

                fclose(fptr);
                fprintf(stdout," %f stdev: %f \n",s1,s2);
                /*if(!my_file_exists(param->in_sequences)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->in_sequences);
                        }*/
        }


        free_wims_thread_data(td);
        thr_pool_destroy(pool);
        free_ihmm_sequences(sb);
        free_fhmm(fhmm);
        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        if(fptr){
                fclose(fptr);
        }

        free_wims_thread_data(td);
        if(pool){
                thr_pool_destroy(pool);
        }
        free_ihmm_sequences(sb);
        free_fhmm(fhmm);
        free_parameters(param);
        return EXIT_FAILURE;
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
        const char usage[] = " -m <model.h5> -i <input sequences> ";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthreads","Number of threads." ,"[8]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--background","Background sequences - residue counts from these will be ADDED to the background model. " ,"[8]"  );
        return OK;
}
