#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <libgen.h>

#include "tldevel.h"

#include "ihmm_seq.h"

#include "finite_hmm.h"

#include "run_score.h"

#include "model.h"

#include "tlmisc.h"
#include "tllogsum.h"


#define SELECT_MODE_BIC 1
#define SELECT_MODE_ML 2



struct parameters{
        char* in_model;
        int num_threads;
        int mode;
};


static int run_model_selection(struct parameters* param);

static int print_help(char **argv);
static int free_parameters(struct parameters* param);

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        //print_program_header(argv, "Build HDPHMM model(s).");

        MMALLOC(param, sizeof(struct parameters));
        param->in_model = NULL;
        param->mode = 0;
        param->num_threads = 4;

        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"nthreads",required_argument,0,'t'},
                        {"ml",0,0,'l'},
                        {"bic",0,0,'b'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"m:t:lbh",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 'm':
                        param->in_model = optarg;
                        break;
                case 'l':
                        param->mode = SELECT_MODE_ML;
                        break;
                case 'b':
                        param->mode = SELECT_MODE_BIC;
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

        if(!param->in_model){
                RUN(print_help(argv));
                ERROR_MSG("No model file! use -m  <blah.h5>");
        }else{
                if(!my_file_exists(param->in_model)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->in_model);
                }
        }

        if(!param->mode){
                 RUN(print_help(argv));
                 ERROR_MSG("No strategy selected. Select either -l (ML) or -b (BIC) to select the best model");
        }


        //RUN(score_sequences_for_command_line_reporting(param));

        RUN(run_model_selection(param));
        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}

int run_model_selection(struct parameters* param)
{
        struct model_bag* model_bag = NULL;

        struct wims_thread_data** td = NULL;

        struct thr_pool* pool = NULL;

        struct seq_buffer* sb = NULL;
        struct fhmm* fhmm = NULL;


        //double** all_scores = NULL;
        //double s1,s2;
        double best_score;
        //double sum;
        double max_likelihood;

        double BIC;
        double total_seq_len;
        //int max = 1000;
        int c;
        int i;
        int limit;

        //char buffer[BUFFER_LEN];



        LOG_MSG("Reading models.");
        RUNP(model_bag = read_model_bag_hdf5(param->in_model));
        LOG_MSG("Done.");

        LOG_MSG("Loading sequences.");

        RUNP(sb =get_sequences_from_hdf5_model(param->in_model, IHMM_SEQ_READ_ONLY_SEQ  ));

        LOG_MSG("Read %d sequences.",sb->num_seq);

        LOG_MSG("Starting thread pool.");
       /* start threadpool  */
        //if((pool = thr_pool_create(param->num_threads , param->num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");


        /* allocate data for threads; */
        RUNP(td = create_wims_thread_data(&param->num_threads,(sb->max_len+2)  , model_bag->max_num_states , NULL));


        LOG_MSG("Done.");

        model_bag->best_model = -1;
        best_score = -100.0;
        max_likelihood = prob2scaledprob(1.0);

        limit = sb->num_seq;

        total_seq_len = 0.0;
        for(i = 0; i < limit;i++){
                total_seq_len += sb->sequences[i]->seq_len;
        }



        //RUNP(all_scores = galloc(all_scores, limit,model_bag->num_models, 0.0));

        for(c = 0; c < model_bag->num_models;c++){
                fhmm = model_bag->finite_models[c];
                RUN(alloc_dyn_matrices(fhmm));
                RUN(run_score_sequences(fhmm,sb,td));

                max_likelihood = prob2scaledprob(1.0);
                for(i = 0; i < limit;i++){
                        max_likelihood = logsum(max_likelihood,sb->sequences[i]->score);
                }

                RUN(calculate_BIC(fhmm, max_likelihood, total_seq_len, &BIC));
                LOG_MSG("Model: %d ML: %.2f BIC %f (based on first %d seqs)",c,max_likelihood,BIC,limit);
                if(param->mode == SELECT_MODE_BIC){
                        if(BIC > best_score){
                                best_score = BIC;
                                model_bag->best_model = c;
                        }
                }else{
                        if(max_likelihood > best_score){
                                best_score = max_likelihood;
                                model_bag->best_model = c;
                        }
                }
        }

        LOG_MSG("Best Model: %d",model_bag->best_model);
        //RUN(write_best_model(param->output, model_bag->best_model));

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
        free_wims_thread_data(td);
        //LOG_MSG("Got past free thread data ");
        //thr_pool_destroy(pool);
        //LOG_MSG("Got past poolfree");
        free_ihmm_sequences(sb);
        //LOG_MSG("Got past seq free");
        free_model_bag(model_bag);
        //LOG_MSG("Got past model bag free");


        return OK;
ERROR:
        free_wims_thread_data(td);
        //thr_pool_destroy(pool);
        free_ihmm_sequences(sb);
        free_fhmm(fhmm);
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
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--ml","Select model based on maximum likelihood." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--bic","Select model based on the Bayesian information criterion (BIC)." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthreads","Number of threads." ,"[4]"  );
        return OK;
}
