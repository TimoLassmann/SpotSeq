#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <float.h>
#include <libgen.h>

#include "tldevel.h"


#include "randomkit.h"
#include "outdir.h"

#include "matrix_io.h"

#include "benchmark_seq.h"

#define BUFFER_LEN 128



double background_CpG[4] = {0.155,0.341,0.35,0.154};

double background_norm[4] = {0.262,0.246,0.239,0.253};
double background_mix[4] = {0.0,0.0,0.0,0.0};


struct parameters{
        char* outdir;
        char* run_name;
        double max_mix_ratio;
        int number_negative;
        unsigned long seed;
        rk_state rndstate;
};



int standard_MOTIF_challenge(struct parameters* param);
int create_n_random_sequences_norm_dist(struct seq_buffer* sb,struct parameters* param, int N);
int create_n_random_sequences_mix_dist(struct seq_buffer* sb,struct parameters* param, int N);


static int print_help(char **argv);
static int free_parameters(struct parameters* param);

int mix_challenge(struct parameters* param);

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        //print_program_header(argv, "Generates standard challenge motif benchmarks.");

        MMALLOC(param, sizeof(struct parameters));

        param->outdir = NULL;
        param->run_name = NULL;
        param->number_negative = 1000;
        param->max_mix_ratio = 1.0;
        param->seed = 0;


        while (1){
                static struct option long_options[] ={
                        {"max_mix_ratio",required_argument,0,'m'},
                        {"seed",required_argument,0,'s'},
                        {"negative",required_argument,0,'n'},
                        {"out",required_argument,0,'o'},
                        {"run",required_argument,0,'r'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"he:o:l:s:n:r:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 's':
                        param->seed = atoi(optarg);
                        break;
                case 'n':
                        param->number_negative = atoi(optarg);
                        break;
                case 'o':
                        param->outdir = optarg;
                        break;
                case 'r':
                        param->run_name = optarg;
                        break;
                case 'm':
                        param->max_mix_ratio = atof(optarg);
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

        if(!param->outdir){
                RUN(print_help(argv));
                ERROR_MSG("No output file! use -o <blah.fa>");

        }

        if(!param->run_name){
                RUN(print_help(argv));
                ERROR_MSG("No run name! use -r <blah");

        }
        if(param->max_mix_ratio  > 1.0 || param->max_mix_ratio < 0.0){
                RUN(print_help(argv));
                ERROR_MSG("Ratio has to be between 0 and 1 ");
        }


        /* initialise random number generator  */
        if(param->seed){
                rk_seed(param->seed, &param->rndstate);
        }else{
                rk_randomseed(&param->rndstate);
        }



        mix_challenge(param);

        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}


/*

  Here I model sequences switching from a random to a CpG biased nucleotide distribution.

*/


int mix_challenge(struct parameters* param)
{
        struct seq_buffer* sb = NULL;

        char buffer[BUFFER_LEN];

        int i;

        double r;
        ASSERT(param != NULL,"No parameters");


        for(i = 0; i < 4; i++){
                background_mix[i] = background_norm[i] * (1.0 - param->max_mix_ratio) + background_CpG[i] * param->max_mix_ratio;

        }
        /* make sure everything sums to zero  */
        r = 0.0;
        for(i = 0; i < 4;i++){
                r += background_mix[i];
        }
        for(i = 0; i < 4;i++){
                background_mix[i] /= r;
        }



        r = 0.0;
        for(i = 0; i < 4;i++){
                r += background_CpG[i];
        }
        for(i = 0; i < 4;i++){
                background_CpG[i] /= r;
        }

        r = 0.0;
        for(i = 0; i < 4;i++){
                r += background_norm[i];
        }

        for(i = 0; i < 4;i++){
                background_norm[i] /= r;
        }

        //for(i = 0; i < 4;i++){
        //        fprintf(stdout,"%f\t%f\t%f\n",background_CpG[i],background_norm[i],background_mix[i]);
        //}
        /* Step zero cumsum */

        for(i = 1; i < 4;i++){
                background_CpG[i] += background_CpG[i-1];
                background_norm[i] += background_norm[i-1];
                background_mix[i] += background_mix[i-1];

        }

        /* step one create output directories  */
        RUN(create_dir(param->outdir,1));


        RUNP(sb = alloc_seq_buffer( MACRO_MAX(20,param->number_negative)));



        /* create output filename */
        snprintf(buffer, BUFFER_LEN, "%s/train_%s_k%d.fa",param->outdir, param->run_name  ,(int)(param->max_mix_ratio*100));





        /* create training sequence */
        RUN(create_n_random_sequences_mix_dist(sb,param, 20));



        RUN(write_sequences_to_file(sb,buffer));
        reset_sb(sb); // clear out stuff from before...

        /* create testing sequence */
        RUN(create_n_random_sequences_mix_dist(sb,param, param->number_negative));

        snprintf(buffer, BUFFER_LEN, "%s/test_%s_k%d.fa", param->outdir,param->run_name,(int)(param->max_mix_ratio*100));
        RUN(write_sequences_to_file(sb,buffer));



        reset_sb(sb);// clear out stuff from before...


        /* create negative sequences */
        RUN( create_n_random_sequences_norm_dist(sb,param,param->number_negative));
        snprintf(buffer, BUFFER_LEN, "%s/neg_%s_k%d.fa",param->outdir,param->run_name,(int)(param->max_mix_ratio*100));
        RUN(write_sequences_to_file(sb,buffer));




        //snprintf(buffer, BUFFER_LEN, "%s/%s_%s_mis_%d_%d.fa", outdir,"Standard_Challenge", motif,n_mismatches,lf*10);
        //LOG_MSG("Writing to: %s.",buffer);
        //RUN(write_sequences_to_file(sb,buffer));

        free_sb(sb);
        return OK;
ERROR:
        return FAIL;
}

int create_n_random_sequences_norm_dist(struct seq_buffer* sb,struct parameters* param, int N)
{
        struct sequence* sequence = NULL;

        int i,j,c;
        double r;

        ASSERT(sb!= NULL, "No sequence buffer");
        ASSERT(param!= NULL, "No sequence buffer");

        ASSERT(sb->malloc_num >= N, "Not enough space in sequence buffer");
/* reset sb to nothing */
        for(i = 0; i < N;i++){
                sb->seqs[i]->seq_len = 0;
        }
        sb->num_seq = 0;
        for(i = 0; i < N;i++){
                sequence = sb->seqs[i];
                if(sequence->seq_len < 601){
                        sequence->malloc_len = 601;
                        MREALLOC(sequence->seq, sizeof(uint8_t) *sequence->malloc_len);
                }
                for(j = 0; j < 600;j++){
                        r = rk_double(&param->rndstate);
                        for(c = 0; c < 4;c++){
                                if(r <= background_norm[c]){
                                        sequence->seq[sequence->seq_len] = "ACGT"[c];
                                        sequence->seq_len++;

                                        break;
                                }
                        }
                }
                sequence->seq[sequence->seq_len] = 0;
                sequence->seq_len++;
                sb->num_seq++;
        }

        return OK;
ERROR:
        return FAIL;

}

int create_n_random_sequences_mix_dist(struct seq_buffer* sb,struct parameters* param, int N)
{
        struct sequence* sequence = NULL;

        int i,j,c;
        int state;
        double r;

        double t[2][2] = {
                {0.98,1.0},
                {0.02,1.0}
        };


        ASSERT(sb!= NULL, "No sequence buffer");
        ASSERT(param!= NULL, "No sequence buffer");

        ASSERT(sb->malloc_num >= N, "Not enough space in sequence buffer");
/* reset sb to nothing */
        for(i = 0; i < N;i++){
                sb->seqs[i]->seq_len = 0;
        }
        sb->num_seq = 0;


        for(i = 0; i < N;i++){
                sequence = sb->seqs[i];
                if(sequence->seq_len < 601){
                        sequence->malloc_len = 601;
                        MREALLOC(sequence->seq, sizeof(uint8_t) *sequence->malloc_len);
                }
                /* we begin in norm state */
                state =0;
                for(j = 0; j < 600;j++){
                        /* emission */
                        r = rk_double(&param->rndstate);
                        for(c = 0; c < 4;c++){
                                if(!state){
                                        if(r <= background_norm[c]){
                                                sequence->seq[sequence->seq_len] = "ACGT"[c];
                                                sequence->seq_len++;
                                                break;
                                        }
                                }else{
                                        if(r <= background_mix[c]){
                                                sequence->seq[sequence->seq_len] = "ACGT"[c];
                                                sequence->seq_len++;
                                                break;
                                        }
                                }
                        }
                        /* transition */
                        r = rk_double(&param->rndstate);
                        for(c = 0; c < 2;c++){
                                if(r <= t[state][c]){
                                        state = c;
                                        break;
                                }
                        }

                }
                sequence->seq[sequence->seq_len] = 0;
                sequence->seq_len++;
                sb->num_seq++;
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
        const char usage[] = " -in <fasta> -out <outfile>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-o","Target outdir" ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-r","Run name " ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-m","Mixing ratio." ,"[1.0]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-l","Motif location." ,"[100]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-n","Number of negative sequences." ,"[1000]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-s","Seed for random numbers." ,"[0]"  );

        return OK;
}
