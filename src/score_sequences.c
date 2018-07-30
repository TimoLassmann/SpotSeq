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

#include "hmm_score.h"

struct parameters{       
        char* in_model;
        char* in_sequences;
        char* output;
};


static int print_help(char **argv);
static int free_parameters(struct parameters* param);

static int run_score_sequences(struct parameters* param);

int main (int argc, char *argv[]) 
{		
        struct parameters* param = NULL;
        int c;
        
        tlog.echo_build_config();
        
        MMALLOC(param, sizeof(struct parameters));
        param->in_model = NULL;
        param->in_sequences = NULL;
        param->output = NULL;
        
        while (1){	
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"in",required_argument,0,'i'},
                        {"out",required_argument,0,'o'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"hm:i:",long_options, &option_index);
		
                if (c == -1){
                        break;
                }
                switch(c) {
                case 'i':
                        param->in_sequences = optarg;
                        break;
                case 'o':
                        param->output = optarg;
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

        if(!param->output){
                RUN(print_help(argv));
                ERROR_MSG("No model file! use -m  <blah.h5>");
        }else{
                if(my_file_exists(param->output)){
                        //RUN(print_help(argv));
                       WARNING_MSG("The file %s will be over-written.",param->output);                
                }   
        }

        RUN(run_score_sequences(param));
        
        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}

int run_score_sequences(struct parameters* param)
{
        struct fhmm* fhmm = NULL;
        struct seq_buffer* sb = NULL;
        int i;
        int expected_len;
        FILE* fptr = NULL;
        
        ASSERT(param != NULL, "no parameters");
        
        /* get model set up  */
        RUNP(fhmm = init_fhmm(param->in_model));

        /* load sequences.  */
        LOG_MSG("Loading sequences.");
        RUNP(sb = load_sequences(param->in_sequences));

        LOG_MSG("Read %d sequences.",sb->num_seq);
        
        /* allocate dyn programming matrices.  */
        RUN(realloc_dyn_matrices(fhmm, sb->max_len+1));

        /* calculate average length */
        expected_len = 0;
        for(i = 0; i < sb->num_seq;i++){
                expected_len += sb->sequences[i]->seq_len;
        }
        expected_len = expected_len / sb->num_seq;
        LOG_MSG("Average sequence length: %d",expected_len);
        /* score sequences  */
        RUNP(fptr = fopen(param->output, "w"));
        fprintf(fptr, "Name,Score_%s\n",  param->in_model);
        for(i = 0; i < sb->num_seq;i++){
                fprintf(stdout,"Running %d (len: %d) %d%d%d\n",i,sb->sequences[i]->seq_len,sb->sequences[i]->seq[0],sb->sequences[i]->seq[1],sb->sequences[i]->seq[2]);
                RUN(forward(fhmm, sb->sequences[i]->seq, sb->sequences[i]->seq_len));
                RUN(random_model_score(fhmm, sb->sequences[i]->seq, sb->sequences[i]->seq_len, expected_len));
                fprintf(stdout,"%d f:%f  r:%f\tlog-odds:%f\tP(M):%f\n",  i, fhmm->f_score, fhmm->r_score, fhmm->f_score - fhmm->r_score, expf(fhmm->f_score - fhmm->r_score ) /  (1.0 + expf(fhmm->f_score - fhmm->r_score ) ));
                fprintf(fptr, "%s,%f\n",sb->sequences[i]->name,  expf(fhmm->f_score - fhmm->r_score ) /  (1.0 + expf(fhmm->f_score - fhmm->r_score ) ));
        }

        fclose(fptr);
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

        //    	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--seed-step","Distance between seeds." ,"[8]"  );
        return OK;
}

