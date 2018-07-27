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


#include "distributions.h"

struct parameters{       
        char* input;
        char* output;
        char* in_model;
        char* cmd_line;
        int local;
        int num_threads;
        int num_start_states;
};


static int run_build_ihmm(struct parameters* param);
static int run_build_fhmm(struct parameters* param);


static int print_help(char **argv);
static int free_parameters(struct parameters* param);

static double loggam(double x);
double loggam(double x)
{
    double x0, x2, xp, gl, gl0;
    long k, n;

    static double a[10] = {8.333333333333333e-02,-2.777777777777778e-03,
         7.936507936507937e-04,-5.952380952380952e-04,
         8.417508417508418e-04,-1.917526917526918e-03,
         6.410256410256410e-03,-2.955065359477124e-02,
         1.796443723688307e-01,-1.39243221690590e+00};
    x0 = x;
    n = 0;
    if ((x == 1.0) || (x == 2.0))
    {
        return 0.0;
    }
    else if (x <= 7.0)
    {
        n = (long)(7 - x);
        x0 = x + n;
    }
    x2 = 1.0/(x0*x0);
    xp = 2*M_PI;
    gl0 = a[9];
    for (k=8; k>=0; k--)
    {
        gl0 *= x2;
        gl0 += a[k];
    }
    gl = gl0/x0 + 0.5*log(xp) + (x0-0.5)*log(x0) - x0;
    if (x <= 7.0)
    {
        for (k=1; k<=n; k++)
        {
            gl -= log(x0-1.0);
            x0 -= 1.0;
        }
    }
    return gl;
}

int main (int argc, char *argv[]) 
{		
        struct parameters* param = NULL;
        int c;
        
        tlog.echo_build_config();
        /*int i;
        for(i = 0 ; i < 100;i++){
                fprintf(stdout,"%d %f\n",i,loggam((double) i / 10.0f));
                }
        exit(0);
        */
        MMALLOC(param, sizeof(struct parameters));
        param->input = NULL;
        param->output = NULL;
        param->in_model = NULL;
        param->cmd_line = NULL;
        param->num_threads = 8;
        param->num_start_states = 10;
        param->local = 0;

                
        while (1){	
                static struct option long_options[] ={
                        {"in",required_argument,0,'i'},
                        {"out",required_argument,0,'o'},
                        {"states",required_argument,0,'s'},
                        {"local",no_argument,0,'l'},
                        {"nthreads",required_argument,0,'n'},
                        {"model",required_argument,0,'m'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"hi:o:n:m:s:l",long_options, &option_index);
		
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

        RUN(run_build_fhmm(param));
        //RUN(seed_controller_thread(param));
        
        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}

int run_build_fhmm(struct parameters* param)
{
        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;
        
        int initial_states = 10;
        int iter;
        int i,j,c;
        float** s1_e = NULL;
        float** s2_e = NULL;

        float** s1_t = NULL;
        float** s2_t = NULL;
        float sum;
        int iterations = 1000;
        
        ASSERT(param!= NULL, "No parameters found.");
        
   
        RUNP(model = read_model_hdf5(param->output));
        RUNP(ft = alloc_fast_hmm_param(initial_states,model->L));

        /* first index is state * letter ; second is sample (max = 100) */

        s1_e = malloc_2d_float(s1_e, model->num_states, model->L, 0.0);
        s2_e = malloc_2d_float(s2_e, model->num_states, model->L, 0.0);

        s1_t = malloc_2d_float(s1_t, model->num_states, model->num_states, 0.0);
        s2_t = malloc_2d_float(s2_t, model->num_states, model->num_states, 0.0);

        
        for( iter=  0;iter < iterations;iter++){
                RUN(fill_fast_transitions_only_matrices(model,ft));
                for(i = 0;i < model->num_states;i++){
                        for(c = 0; c < model->L;c++){
                                s1_e[i][c] += ft->emission[c][i];
                                s2_e[i][c] += (ft->emission[c][i] * ft->emission[c][i]);
                        }
                }
                for(i = 0;i < model->num_states;i++){
                        for(c = 0;c < model->num_states;c++){
                                s1_t[i][c] += ft->transition[i][c];
                                s2_t[i][c] += (ft->transition[i][c] * ft->transition[i][c]);
                        }
                }
                
        }
        
        for(i = 0; i < model->num_states;i++){
                sum = 0.0;
                for(j = 0; j < model->L;j++){
                        s2_e[i][j] = sqrt(  ((double) iterations * s2_e[i][j] - s1_e[i][j] * s1_e[i][j])/ ((double) iterations * ((double) iterations -1.0)));
                        s1_e[i][j] = s1_e[i][j] / (double) iterations;
                        sum+= s1_e[i][j];
                        //fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_e[i][j], s2_e[i][j]);
                }
                if(sum){
                        fprintf(stdout,"Emission:%d\n",i);
                        for(j = 0; j < model->L;j++){
                                s1_e[i][j] /= sum;
                                fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_e[i][j], s2_e[i][j]);
                        }
                }

                sum = 0;
                for(j = 0;j < model->num_states;j++){
                        s2_t[i][j] = sqrt(  ((double) iterations * s2_t[i][j] - s1_t[i][j] * s1_t[i][j])/ ((double) iterations * ((double) iterations -1.0)));
                        s1_t[i][j] = s1_t[i][j] / (double) iterations;
                        sum+= s1_t[i][j];
                        //fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_t[i][j], s2_t[i][j]);
                }
                if(sum){
                        fprintf(stdout,"transition:%d\n",i);
                        for(j = 0;j < model->num_states;j++){
                                s1_t[i][j] /= sum;
                                fprintf(stdout,"%d %d : %f stdev:%f\n",i,j,s1_t[i][j], s2_t[i][j]);
                        }
                }
                
        }

        RUN(add_fhmm(param->output,s1_t,s1_e, model->num_states, model->L  ));
        
        free_2d((void**) s1_e);
        free_2d((void**) s2_e);

        free_2d((void**) s1_t);
        free_2d((void**) s2_t);
   
        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        return OK;
ERROR:
        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        return FAIL;
}


int run_build_ihmm(struct parameters* param)
{
        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;
        struct seq_buffer* sb = NULL;
        
        int initial_states = 400;
        ASSERT(param!= NULL, "No parameters found.");
        
        initial_states = param->num_start_states;
        
        /* Step one read in sequences */
        LOG_MSG("Loading sequences.");


        RUNP(sb = load_sequences(param->input));

        
        LOG_MSG("Read %d sequences.",sb->num_seq);
       
        if(param->in_model){
                RUNP(model = read_model(param->in_model));
        }else{
                RUNP(model = alloc_ihmm_model(initial_states, sb->L));
                /* Initial guess... */
                model->alpha0_a = 6.0f;
                model->alpha0_b = 15.0f;
                model->gamma_a = 16.0f;
                model->gamma_b = 4.0f;
                model->alpha = IHMM_PARAM_PLACEHOLDER;
                model->gamma = IHMM_PARAM_PLACEHOLDER;
                model->target_len = param->local;
                RUN(inititalize_model(model, sb,initial_states));// initial_states) );
                
        }
        RUNP(ft = alloc_fast_hmm_param(initial_states,sb->L));
        RUN(fill_background_emission(ft, sb));
        RUN(run_beam_sampling( model, sb, ft,NULL, 10000, 10));

        //RUN(write_model(model, param->output));
        RUN(write_model_hdf5(model, param->output));
//        RUN(add_annotation)
        RUN(add_annotation(param->output, "spotseq_model_cmd", param->cmd_line));

        RUN(print_states_per_sequence(sb));
        RUN(write_ihmm_sequences(sb,"test.lfs","testing"));
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
        //    	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--seed-step","Distance between seeds." ,"[8]"  );
//      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthread","Number of threads." ,"[8]"  );
//      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--output","Output file name." ,"[?]"  );
        return OK;
}

