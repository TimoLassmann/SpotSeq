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
#include "tlmisc.h"




#include "model.h"

#include "adjusted_rand_index.h"



struct parameters{
        char* in_model;
};



int calculate_pair_wise_ari(struct parameters* param);


int free_parameters(struct parameters* param);
int print_help(char **argv);

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;


        //print_program_header(argv, "Build HDPHMM model(s).");

        MMALLOC(param, sizeof(struct parameters));
        param->in_model = NULL;

        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"m:h",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
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

        if(!param->in_model){
                RUN(print_help(argv));
                ERROR_MSG("No model file! use -m  <blah.h5>");
        }else{
                if(!my_file_exists(param->in_model)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->in_model);
                }
        }




        calculate_pair_wise_ari(param);

        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}

int calculate_pair_wise_ari(struct parameters* param)
{
        struct seq_buffer* sb = NULL;

        struct model_bag* mb = NULL;

        long long** count_matrix = NULL;
        double adjRandIndex = 0.0;
        int** tmp_l;

        int* number_of_states = NULL;
        int num_models = -1;
        int max_K = -1;
        int i,j,c,g;
        LOG_MSG("Reading model bag.");
        RUNP(mb = read_model_bag_hdf5(param->in_model));
        LOG_MSG("Done.");

        /* want to store number of models and number of states in each  */
        num_models = mb->num_models;

        ASSERT(num_models != -1, "No models found!");

        max_K = mb->max_num_states;
        ASSERT(max_K > 0, "No state in any model.");


        MMALLOC(number_of_states, sizeof(int) * num_models);
        for(i = 0; i < num_models;i++){
                number_of_states[i] = mb->models[i]->num_states;
        }

        free_model_bag(mb);


        LOG_MSG("Reading labelled sequences.");
        RUNP(sb = get_sequences_from_hdf5_model(param->in_model, IHMM_SEQ_READ_ALL));
        LOG_MSG("Done.");

        MMALLOC(count_matrix, sizeof(long long*) * max_K);
        for(i =0;i < max_K;i++){
                count_matrix[i] = NULL;
                MMALLOC(count_matrix[i], sizeof(long long) * max_K);
        }
        for(i = 0; i < num_models;i++){
                for(j = i;j< num_models;j++){
                        for(c = 0; c < max_K;c++){
                                for(g = 0; g < max_K;g++){
                                        count_matrix[c][g] = 0;
                                }
                        }

                        for(c = 0; c < sb->num_seq;c++){
                                tmp_l = sb->sequences[c]->label_arr;
                                for(g = 0;g < sb->sequences[c]->seq_len;g++){
                                        count_matrix[tmp_l[i][g]-2][tmp_l[j][g]-2]++;
                                }
                        }
                        for(c = 0; c < number_of_states[i];c++){
                                for(g = 0;g < number_of_states[j];g++){
                                        //             fprintf(stdout,"%lld ", count_matrix[c][g]);
                                }
                                //fprintf(stdout,"\n");
                        }

                        //  fprintf(stdout,"\n");
                        RUN(ari(count_matrix,   number_of_states[i], number_of_states[j] , &adjRandIndex));
                        LOG_MSG("%d-%d: %f",i,j, adjRandIndex);
                }

        }
        for(i =0;i < max_K;i++){
                MFREE(count_matrix[i]);
        }
        MFREE(count_matrix);


        MFREE(number_of_states);

        free_ihmm_sequences(sb);

        return OK;
ERROR:
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
        const char usage[] = " -m <model.h5>  ";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--model","Input model." ,"[]"  );
        return OK;
}
