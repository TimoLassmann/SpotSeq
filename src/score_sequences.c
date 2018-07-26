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

#include "hdf5_glue.h"

struct parameters{       
        char* in_model;
        char* in_sequences;
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
                
        while (1){	
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"in",required_argument,0,'i'},
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
        struct hdf5_data* hdf5_data = NULL;
        float** e = NULL;
        float** t = NULL;
        int num_states;
        int l;
        int i,j;
        ASSERT(param != NULL, "no parameters");

        /* read in hdf5 file and get emission and transition matrix */
        ASSERT(param->in_model != NULL, "No filename");
        ASSERT(my_file_exists(param->in_model) != 0,"File %s does not exist.",param->in_model);


        hdf5_data = hdf5_create();
	
        hdf5_open_file(param->in_model,hdf5_data);

        hdf5_read_attributes(hdf5_data,hdf5_data->file);
        ASSERT(hdf5_data->num_attr != 0 , "Could not find attributes");
        print_attributes(hdf5_data);
        get_group_names(hdf5_data);
        fprintf(stdout,"Groups:\n");
        for(i = 0; i < hdf5_data->grp_names->num_names;i++){
                fprintf(stdout,"%d %s\n",i,hdf5_data->grp_names->names[i]);
        }
        
        hdf5_open_group("imodel",hdf5_data);
        hdf5_read_attributes(hdf5_data,hdf5_data->group);
        ASSERT(hdf5_data->num_attr != 0 , "Could not find attributes");
        print_attributes(hdf5_data);
        num_states = 0;
        l = 0;
        for(i = 0; i < hdf5_data->num_attr;i++){
                if(!strncmp("Number of states", hdf5_data->attr[i]->attr_name, 16)){
                        num_states = hdf5_data->attr[i]->int_val;
                }
                if(!strncmp("Number of letters", hdf5_data->attr[i]->attr_name, 17)){
                        l = hdf5_data->attr[i]->int_val;
                }
                
        }
        hdf5_close_group(hdf5_data);
        hdf5_open_group("fmodel",hdf5_data);

        hdf5_read_dataset("emission",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read transition_counts");
        e = (float**) hdf5_data->data;

        
        hdf5_read_dataset("transition",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read transition_counts");
        t = (float**) hdf5_data->data;

        hdf5_close_group(hdf5_data);

        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);


        for(i = 0; i < num_states;i++){
                for(j =0; j < l;j++){
                        fprintf(stdout,"%f ",e[i][j]);
                }
                fprintf(stdout,"\n");
        }
        
        return OK;
ERROR:
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
