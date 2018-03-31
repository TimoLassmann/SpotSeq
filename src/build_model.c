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


struct parameters{       
        char* input;
        char* output;
        char* in_model;
        int num_threads;
};


static int run_build_ihmm(struct parameters* param);

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
        param->num_threads = 8;
                
        while (1){	
                static struct option long_options[] ={
                        {"in",required_argument,0,'i'},
                        {"out",required_argument,0,'o'},
                        {"nthreads",required_argument,0,'n'},
                        {"model",required_argument,0,'m'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"hi:o:n:m:",long_options, &option_index);
		
                if (c == -1){
                        break;
                }
                switch(c) {
                case 'i':
                        param->input = optarg;
                        break;
                                 
                case 'o':
                        param->output = optarg;
                        break;
                case 'm':
                        param->in_model = optarg;
                        break;
                case 'n':
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
        
        RUN(run_build_ihmm(param));

        //RUN(seed_controller_thread(param));
        
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
        int i;
        int initial_states = 10;
        ASSERT(param!= NULL, "No parameters found.");
        
      
        
        /* Step one read in sequences */
        LOG_MSG("Loading sequences.");


        RUNP(sb = load_sequences(param->input));

        
        LOG_MSG("Read %d sequences.",sb->num_seq);

        if(param->in_model){
                RUNP(model = read_model(param->in_model));
        }else{
                RUNP(model = alloc_ihmm_model(initial_states, sb->L));
                /* Initial guess... */
                model->alpha0_a = 4.0f;
                model->alpha0_b = 2.0f;
                model->gamma_a = 3.0f;
                model->gamma_b = 6.0f;
                model->alpha = IHMM_PARAM_PLACEHOLDER;
                model->gamma = IHMM_PARAM_PLACEHOLDER;        
                RUN(inititalize_model(model, sb, 0) );
        }
        RUNP(ft = alloc_fast_hmm_param(initial_states,sb->L));
        RUN(fill_background_emission(ft, sb));
        RUN(run_beam_sampling( model, sb, ft,NULL, 1000, 10));

        RUN(write_model(model, param->output));
        
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
        //    	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--seed-step","Distance between seeds." ,"[8]"  );
//      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthread","Number of threads." ,"[8]"  );
//      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--output","Output file name." ,"[?]"  );
        return OK;
}

