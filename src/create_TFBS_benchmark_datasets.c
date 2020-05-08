
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <libgen.h>

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <float.h>


#include "tldevel.h"


#include "randomkit.h"
#include "outdir.h"

#include "matrix_io.h"

#include "benchmark_seq.h"

double background[4] = {0.25,0.5,0.75,1.0};

#define BUFFER_LEN 128


struct parameters{
        char* outdir;
        int error;
        int motif_location;
        int number_negative;
        unsigned long seed;
        rk_state rndstate;
};

/* Structs to hold sequences  */


int standard_MOTIF_challenge(struct parameters* param);
int create_n_random_sequences(struct seq_buffer* sb,struct parameters* param, int N);
int add_motif_to_all_sequences(struct seq_buffer* sb,struct parameters* param, char* motif);

static int print_help(char **argv);
static int free_parameters(struct parameters* param);




int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        //print_program_header(argv, "Generates standard challenge motif benchmarks.");

        MMALLOC(param, sizeof(struct parameters));

        param->outdir = NULL;

        param->number_negative = 1000;
        param->error = 0;
        param->motif_location = 100;
        param->seed = 0;

        while (1){
                static struct option long_options[] ={
                        {"error",required_argument,0,'e'},
                        {"location",required_argument,0,'l'},
                        {"seed",required_argument,0,'s'},
                        {"negative",required_argument,0,'n'},
                        {"out",required_argument,0,'o'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"he:o:l:s:n:",long_options, &option_index);

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

                case 'l':
                        param->motif_location = atoi(optarg);
                        break;
                case 'o':
                        param->outdir = optarg;
                        break;
                case 'e':
                        param->error = atoi(optarg);
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
        if(param->error > 15){
                RUN(print_help(argv));
                ERROR_MSG("Too many errors");
        }


        /* initialise random number generator  */
        if(param->seed){
                rk_seed(param->seed, &param->rndstate);
        }else{
                rk_randomseed(&param->rndstate);
        }
        standard_MOTIF_challenge(param);

        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
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
        const char usage[] = " -in <fasta> -out <outfile>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-o","Target outdir" ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-e","Number of errors." ,"[0]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-l","Motif location." ,"[100]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-n","Number of negative sequences." ,"[1000]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-s","Seed for random numbers." ,"[0]"  );

        return OK;
}
/*

From:

Redhead, Emma, and Timothy L. Bailey. "Discriminative motif discovery
in DNA and protein sequences using the DEME algorithm." BMC
bioinformatics 8.1 (2007): 385.


The idea is derived from the so-called "standard challenge problem"
introduced by Pevzner et al. [29] as a way of testing
non-discriminative motif finders. The standard challenge problem
specifies a synthetic DNA dataset consisting of 20 length-600
sequences, each containing an artificially generated motif
occurrence. The motif is represented a string of length 15, and each
occurrence contains exactly four mismatches. We augment the standard
challenge problem by including a set of negative sequences and
define four synthetic discriminative motif discovery problems.

Motif occurrences are generated by mutating exactly d positions in the
consensus, selected at independently and at random with a uniform
probability that a position is mutated. We examine values of d in the
range [0,5].


*/


int standard_MOTIF_challenge(struct parameters* param)
{
        struct seq_buffer* sb = NULL;

        char buffer[BUFFER_LEN];

        char motif[16];
        int i,c;

        double r;




        ASSERT(param != NULL,"No parameters");


        /* step one create output directories  */
        RUN(create_dir(param->outdir,1));


        RUNP(sb = alloc_seq_buffer( MACRO_MAX(20,param->number_negative)));

        while(1){
                /* make motif */
                for(i = 0; i < 15;i++){
                        r =  rk_double(&param->rndstate);
                        for(c = 0; c < 4;c++){
                                if(r <= background[c]){
                                        motif[i] = "ACGT"[c];
                                        break;
                                }
                        }
                }
                motif[15] = 0;
                /* create output filename */
                snprintf(buffer, BUFFER_LEN, "%s/train_%s_k%d.fa",param->outdir,  motif,param->error);

                /* if filename already exists skip and repeat */

                if(!my_file_exists(buffer)){

                        /* create training sequence */
                        RUN(create_n_random_sequences(sb,param, 20));
                        RUN(add_motif_to_all_sequences(sb, param, motif));


                        RUN(write_sequences_to_file(sb,buffer));
                        reset_sb(sb); // clear out stuff from before...

                        /* create testing sequence */
                        RUN(create_n_random_sequences(sb,param, param->number_negative));
                        RUN(add_motif_to_all_sequences(sb, param, motif));
                        snprintf(buffer, BUFFER_LEN, "%s/test_%s_k%d.fa", param->outdir,motif,param->error);
                        RUN(write_sequences_to_file(sb,buffer));



                        reset_sb(sb);// clear out stuff from before...


                        /* create negative sequences */
                        RUN(create_n_random_sequences(sb,param,param->number_negative));
                        snprintf(buffer, BUFFER_LEN, "%s/neg_%s_k%d.fa",param->outdir,motif,param->error);
                        RUN(write_sequences_to_file(sb,buffer));
                        break;
                }

        }


        //snprintf(buffer, BUFFER_LEN, "%s/%s_%s_mis_%d_%d.fa", outdir,"Standard_Challenge", motif,n_mismatches,lf*10);
        //LOG_MSG("Writing to: %s.",buffer);
        //RUN(write_sequences_to_file(sb,buffer));

        free_sb(sb);
        return OK;
ERROR:
        return FAIL;
}

int add_motif_to_all_sequences(struct seq_buffer* sb,struct parameters* param, char* motif)
{
        struct sequence* sequence = NULL;
        char motif_mutated[16];

        int i,j,c,f;
        double r;
        double sum;
        int sample_window;
        int upstream;

        ASSERT(sb!= NULL, "No sequence buffer");
        ASSERT(param!= NULL, "No sequence buffer");

        ASSERT(sb->num_seq >= 0, "Not sequences in buffer.");

        for(i = 0; i < sb->num_seq;i++){
                sequence = sb->seqs[i];
                /* copy original motif */
                for(j = 0; j < 15;j++){
                        motif_mutated[j] = motif[j];
                }
                /* mutate until exactly 4 mismatches */
                //f

                for(j = 0; j < param->error;j++){

                        r = rk_double(&param->rndstate);
                        sum = 1.0 / 15.0;
                        for(f = 0; f < 15;f++){

                                if(r <= sum){
                                        r = rk_double(&param->rndstate);
                                        for(c = 0; c < 4;c++){
                                                if(r <= background[c]){
                                                        fprintf(stdout,"Mutating: %d %c -> %c\n",f,motif_mutated[f], "ACGT"[c]);
                                                        motif_mutated[f] =  "ACGT"[c];
                                                        f = 100;
                                                        c = 100;

                                                }
                                        }

                                }
                                sum += 1.0 / 15.0;
                        }

                }

                fprintf(stdout, "Here we go:\n");
                for(c = 0; c < 15;c++){
                        fprintf(stdout, "%c", motif[c]);
                }
                fprintf(stdout, "\n");

                for(c = 0; c < 15;c++){
                        if(motif_mutated[c] == motif[c]){
                                fprintf(stdout, "|");
                        }else{
                                fprintf(stdout, " ");
                        }
                }

                fprintf(stdout, "\n");
                for(c = 0; c < 15;c++){
                        fprintf(stdout, "%c", motif_mutated[c]);
                }
                fprintf(stdout, "\n");

                sample_window = (int) ((600.0-15.0) * (float) ( param->motif_location) / 100.0f);
                upstream = (int) ((600-15) - sample_window)/2;
                c = (rk_random(&param->rndstate) % sample_window) + upstream;
                fprintf(stdout,"%d %d \n",c,sample_window);
                for(j = 0; j < 15;j++){
                        sequence->seq[c+j] = motif_mutated[j];
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int create_n_random_sequences(struct seq_buffer* sb,struct parameters* param, int N)
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

                for(j = 0; j < 600;j++){
                        r = rk_double(&param->rndstate);
                        for(c = 0; c < 4;c++){
                                if(r <= background[c]){
                                        sequence->seq[sequence->seq_len] = "ACGT"[c];
                                        sequence->seq_len++;
                                        if(sequence->seq_len == sequence->malloc_len){
                                                sequence->malloc_len = sequence->malloc_len << 1;
                                                MREALLOC(sequence->seq, sizeof(uint8_t) *sequence->malloc_len);
                                        }
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

