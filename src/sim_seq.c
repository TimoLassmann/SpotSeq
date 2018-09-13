
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

#include "tldevel.h"

#include "outdir.h"

#include "matrix_io.h"

#include "make_dot_file.h"

struct parameters{
        char* input;
        char* outdir;
        int len;
        int num_seq;
};

/* Structs to hold sequences  */

struct sequence{
        uint8_t* seq;
        char* name;
        int malloc_len;
        int seq_len;
};

struct seq_buffer{
        struct sequence** seqs;
        int malloc_num;
        int num_seq;
};

/* hmm strtuct */

#define STARTSTATE 0
#define ENDSTATE 1

/* sum of special states above.. */
#define NUM_ADDITIONAL_STATES 2

struct hmm{
        float** transitions;
        float** emissions;
        float* background;
        int num_states;
        int alphabet_len;
        int max_seq_len;
};


static int run_sim_seq(struct parameters* param);

int embedded(struct parameters* param, struct seq_buffer*sb);
int init_random(struct parameters* param, struct seq_buffer*sb );
int add_motif(struct seq_buffer*sb, char* motif);

int standard_challenge(char* outdir);
int standard_challenge_split(char* outdir, int gap);
int ACGT_concat_example(struct parameters* param, struct seq_buffer* sb,float mainres_emission,float self_transition);

int ACGT_embedded_example(struct parameters* param, struct seq_buffer* sb,float mainres_emission);

int two_state_example(struct parameters* param, struct seq_buffer* sb);

static int print_help(char **argv);
static int free_parameters(struct parameters* param);


/* seqbuffer stuff  */
static struct seq_buffer* alloc_seq_buffer(int num_seq);
static int write_sequences_to_file(struct seq_buffer* sb,char* filename);
static int reset_sb(struct seq_buffer* sb);
static void free_sb(struct seq_buffer* sb);

static struct hmm* malloc_hmm(int num_states, int alphabet_len, int max_seq_len);
static void free_hmm(struct hmm* hmm);

static int emit_sequence(struct hmm* hmm, struct sequence* sequence);

static int sanity_check_hmm(struct hmm*  hmm);
static int cumsum_hmm(struct hmm* hmm);

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        tlog.echo_build_config();

        MMALLOC(param, sizeof(struct parameters));
        param->num_seq = 1000;
        param->len = 1000;
        param->outdir = NULL;

        while (1){
                static struct option long_options[] ={
                        {"num",required_argument,0,'n'},
                        {"len",required_argument,0,'l'},
                        {"out",required_argument,0,'o'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"hn:l:o:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 'n':
                        param->num_seq = atoi(optarg);
                        break;

                case 'l':
                        param->len = atoi(optarg);
                        break;
                case 'o':
                        param->outdir = optarg;
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

        if(!param->num_seq){
                RUN(print_help(argv));
                ERROR_MSG("Numseq is 0! use --n 1 (or more!).");

        }
        if(!param->outdir){
                RUN(print_help(argv));
                ERROR_MSG("No output file! use -o <blah.fa>");

        }

        ASSERT(param->len > 10, "Simulated sequence length has to be > 10");

        standard_challenge(param->outdir);
        standard_challenge_split(param->outdir, 10);

        RUN(run_sim_seq(param));
        //RUN(run_build_ihmm(param));

        //RUN(seed_controller_thread(param));

        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}

int run_sim_seq(struct parameters* param)
{
        char buffer[BUFFER_LEN];
        struct seq_buffer* sb = NULL;
        FILE* f_ptr = NULL;
        int i;
        int j;
        ASSERT(param!= NULL, "No parameters found.");


        snprintf(buffer, BUFFER_LEN, "%s/" ,param->outdir);

        RUN(create_dir(buffer,1));

        //RUN(create_output_directories(param->outdir));

        snprintf(buffer, BUFFER_LEN, "%s/%s.log",param->outdir,"simlog");

        RUNP(f_ptr = fopen(buffer, "w"));
        fclose(f_ptr);

        tlog.set_logfile(buffer);

        // RUN(set_log_file(param->outdir,"sim_log"));

        //snprintf(buffer, BUFFER_LEN, "%s/%s/",param->outdir,OUTDIR_CHECKPOINTS);
        //DECLARE_CHK(MAIN_CHECK, buffer);

        /* Step one allocate seq struct.. */
        RUNP(sb = alloc_seq_buffer(param->num_seq));


        /* make sequence test sets; each should have a fasta file, the model
         * parameter matrix and a dot file */

        /* CpG example from black book */

        RUN(two_state_example(param,sb));
        RUN(reset_sb(sb));


        for(i = 25;i <= 100;i+=5){
                for(j = 0;j < 50;j+=5){
                        RUN(ACGT_concat_example(param,sb, (float) i / 100.0,(float) j / 100.0));
                        RUN(reset_sb(sb));
                }
        }

        for(i = 100;i <= 100;i+=5){

                RUN(ACGT_embedded_example(param,sb, (float) i / 100.0));
                RUN(reset_sb(sb));

        }
        RUN(reset_sb(sb));
        RUN(embedded(param,sb));
            free_sb(sb);
        //DESTROY_CHK(MAIN_CHECK);

        return OK;
ERROR:
        free_sb(sb);
        return FAIL;
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


int standard_challenge(char* outdir)
{
        struct seq_buffer*sb = NULL;
        struct sequence* sequence = NULL;
        float background[4] = {0.25f,0.5f,0.75f,1.0f};
        char buffer[BUFFER_LEN];
        char motif[16];
        char motif_mutated[16];
        int i,j,c;
        int n_mismatches;
        float r;
        int lf = 0;
        int sample_window;
        int upstream;

        ASSERT(outdir != NULL,"No parameters");
        RUNP(sb = alloc_seq_buffer(20));

        /* make motif */
        for(i = 0; i < 15;i++){
                r = random_float_zero_to_x(1.0);
                for(c = 0; c < 4;c++){
                        if(r <= background[c]){
                                motif[i] = "ACGT"[c];
                                break;
                        }
                }
        }
        motif[15] = 0;
        for(lf = 1;lf <= 10;lf++){

                /* Loop through mismatches */
                for(n_mismatches = 0; n_mismatches < 10; n_mismatches++ ){
                        /* reset sb to nothing */
                        for(i = 0; i < 20;i++){
                                sb->seqs[i]->seq_len = 0;
                        }
                        sb->num_seq = 0;
                        for(i = 0; i < 20;i++){
                                sequence = sb->seqs[i];

                                for(j = 0; j < 600;j++){
                                        r = random_float_zero_to_x(1.0);
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


                        for(i = 0; i < 20;i++){
                                sequence = sb->seqs[i];
                                /* copy motif */
                                for(j = 0; j < 15;j++){
                                        motif_mutated[j] = motif[j];
                                }
                                /* mutate until exactly 4 mismatches */
                                //f


                                for(c = 0; c < n_mismatches; c ++){
                                        j = random_int_zero_to_x(14);
                                        motif_mutated[j] = "ACGT"[random_int_zero_to_x(3)];
                                }
                                /*fprintf(stdout, "Here we go:\n");
                                for(j = 0; j < 15;j++){
                                        fprintf(stdout, "%c", motif[j]);
                                }
                                fprintf(stdout, "\n");

                                for(j = 0; j < 15;j++){
                                        if(motif_mutated[j] == motif[j]){
                                                fprintf(stdout, "|", motif[j]);
                                        }else{
                                                fprintf(stdout, " ", motif[j]);
                                        }
                                }

                                fprintf(stdout, "\n");
                                for(j = 0; j < 15;j++){
                                        fprintf(stdout, "%c", motif_mutated[j]);
                                }
                                fprintf(stdout, "\n");*/

                                sample_window = (int) ((600.0-15.0) * (float) (lf) / 10.0f);
                                upstream = (int) ((600-15) - sample_window)/2;
                                c = random_int_zero_to_x(sample_window) + upstream;

                                for(j = 0; j < 15;j++){
                                        sequence->seq[c+j] = motif_mutated[j];
                                }
                        }

                        snprintf(buffer, BUFFER_LEN, "%s/%s_%s_mis_%d_%d.fa", outdir,"Standard_Challenge", motif,n_mismatches,lf*10);
                        LOG_MSG("Writing to: %s.",buffer);
                        RUN(write_sequences_to_file(sb,buffer));
                }
        }
        free_sb(sb);
        return OK;
ERROR:
        return FAIL;
}

int standard_challenge_split(char* outdir, int gap)
{
        struct seq_buffer*sb = NULL;
        struct sequence* sequence = NULL;
        float background[4] = {0.25f,0.5f,0.75f,1.0f};
        char buffer[BUFFER_LEN];
        char motif[16];
        char motif_mutated[16];
        int i,j,c,g;
        int n_mismatches;
        float r;
        int sample_window;

        ASSERT(outdir != NULL,"No parameters");
        RUNP(sb = alloc_seq_buffer(20));

        /* make motif */
        for(i = 0; i < 15;i++){
                r = random_float_zero_to_x(1.0);
                for(c = 0; c < 4;c++){
                        if(r <= background[c]){
                                motif[i] = "ACGT"[c];
                                break;
                        }
                }
        }
        motif[15] = 0;
        /* Loop through mismatches */
        for(n_mismatches = 0; n_mismatches < 10; n_mismatches++ ){
                /* reset sb to nothing */
                for(i = 0; i < 20;i++){
                        sb->seqs[i]->seq_len = 0;
                }
                sb->num_seq = 0;
                for(i = 0; i < 20;i++){
                        sequence = sb->seqs[i];

                        for(j = 0; j < 600;j++){
                                r = random_float_zero_to_x(1.0);
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


                for(i = 0; i < 20;i++){
                        sequence = sb->seqs[i];
                        /* copy motif */
                        for(j = 0; j < 15;j++){
                                motif_mutated[j] = motif[j];
                        }
                        /* mutate until exactly 4 mismatches */
                        //

                        for(c = 0; c < n_mismatches; c ++){
                                j = random_int_zero_to_x(14);
                                motif_mutated[j] = "ACGT"[random_int_zero_to_x(3)];
                        }
                        /*fprintf(stdout, "Here we go:\n");
                        for(j = 0; j < 15;j++){
                                fprintf(stdout, "%c", motif[j]);
                        }
                        fprintf(stdout, "\n");

                        for(j = 0; j < 15;j++){
                                if(motif_mutated[j] == motif[j]){
                                        fprintf(stdout, "|", motif[j]);
                                }else{
                                        fprintf(stdout, " ", motif[j]);
                                }
                        }

                        fprintf(stdout, "\n");
                        for(j = 0; j < 15;j++){
                                fprintf(stdout, "%c", motif_mutated[j]);
                        }
                        fprintf(stdout, "\n");*/

                        sample_window = (int) ((600.0-(15.0+gap)));

                        c = random_int_zero_to_x(sample_window);

                        for(j = 0; j < 15;j++){
                                g = 0;
                                if(j >= 8){
                                        g = gap;
                                }
                                sequence->seq[c+j+g] = motif_mutated[j];
                        }
                }

                snprintf(buffer, BUFFER_LEN, "%s/%s_%s_mis_%d_%d.fa", outdir,"Standard_Challenge_gap", motif,n_mismatches,gap);
                LOG_MSG("Writing to: %s.",buffer);
                RUN(write_sequences_to_file(sb,buffer));
        }

        free_sb(sb);
        return OK;
ERROR:
        return FAIL;
}


int embedded(struct parameters* param, struct seq_buffer*sb)
{
        char buffer[BUFFER_LEN];
        char motif[BUFFER_LEN];

        ASSERT(sb != NULL,"No seq buffer");
        ASSERT(param != NULL,"No seq buffer");

         snprintf(buffer, BUFFER_LEN, "%s/%s.fa",param->outdir,"RANDOM");
        LOG_MSG("Writing to: %s.",buffer);

        RUN(reset_sb(sb));
        RUN(init_random(param,sb));
        RUN(write_sequences_to_file(sb,buffer));

        snprintf(motif, BUFFER_LEN, "%s","GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC");
        RUN(add_motif(sb,motif));
        snprintf(buffer, BUFFER_LEN, "%s/%s_%s.fa",param->outdir,"RANDOM",motif);
        LOG_MSG("Writing to: %s.",buffer);
        RUN(write_sequences_to_file(sb,buffer));


        snprintf(motif, BUFFER_LEN, "%s","TGCATGCATGCA");
        RUN(add_motif(sb,motif));
        snprintf(buffer, BUFFER_LEN, "%s/%s_%s.fa",param->outdir,"RANDOM",motif);
        LOG_MSG("Writing to: %s.",buffer);
        RUN(write_sequences_to_file(sb,buffer));



        return OK;
ERROR:
        return FAIL;

}

int init_random(struct parameters* param, struct seq_buffer*sb )
{
        struct sequence* sequence = NULL;
        int expected_length = param->len;
        int i,j,c;

        ASSERT(sb != NULL,"No seq buffer");

        ASSERT(sb->seqs[0]->seq_len == 0,"Need to reset seq buffer");


        for(i = 0; i < sb->malloc_num;i++){
                sequence = sb->seqs[i];
                sequence->seq_len = 0;
                for(j = 0; j < expected_length;j++){
                        c = random_int_zero_to_x(3);
                        sequence->seq[sequence->seq_len] = "ACGT"[c];
                        sequence->seq_len++;
                        if(sequence->seq_len == sequence->malloc_len){
                                sequence->malloc_len = sequence->malloc_len << 1;
                                MREALLOC(sequence->seq, sizeof(uint8_t) *sequence->malloc_len);
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

int add_motif(struct seq_buffer*sb, char* motif)
{
        struct sequence* sequence = NULL;
        int motif_len = 0;
        int i,j,c;
        ASSERT(sb != NULL,"No seq buffer");
        ASSERT(motif  != NULL,"No seq buffer");

        motif_len = strlen(motif);

        for(i = 0; i < sb->num_seq;i++){
                if(random_float_zero_to_x(1.0) < 1.0){
                        sequence = sb->seqs[i];
                        ASSERT(sequence->seq_len > motif_len,"Motif is longer than sequence");
                        c = random_int_zero_to_x(sequence->seq_len - (motif_len + 1));
                        for(j = 0; j < motif_len;j++){
                                sequence->seq[j+c] = motif[j];
                        }
                }
        }

        return OK;
ERROR:
        return FAIL;
}



int ACGT_concat_example(struct parameters* param, struct seq_buffer* sb,float mainres_emission,float self_transition)
{
        char buffer[BUFFER_LEN];
        struct hmm* hmm = NULL;
        float leave;
        float sum;
        int i,j;
        int expected_length = param->len;


        ASSERT(sb != NULL,"No seq buffer");

        ASSERT(sb->seqs[0]->seq_len == 0,"Need to reset seq buffer");

        ASSERT(mainres_emission <= 1.0,"Main emission has to be smaller than 1.0.");
        ASSERT(mainres_emission >= 0.0,"Main emission has to be greater than 0.0.");

        ASSERT(self_transition <= 1.0,"Main emission has to be smaller than 1.0.");
        ASSERT(self_transition >= 0.0,"Main emission has to be greater than 0.0.");


        RUNP(hmm = malloc_hmm(4,4,expected_length));

        /*
          NOTE This should be:
          leave = 1.0 - (float) ((expected_length) /(float)(1+   expected_length));
          BUT:
          I will substract one from the expected length because I co not allow a transition from start to stop (i.e. a zero length sequence
        */

        leave = 1.0 - (float) ((expected_length) /(float)(   expected_length+1));

        /* set transition to end to 1/ expected lengeth  */

        hmm->transitions[STARTSTATE][2] = 0.25f;
        hmm->transitions[STARTSTATE][3] = 0.25f;
        hmm->transitions[STARTSTATE][4] = 0.25f;
        hmm->transitions[STARTSTATE][5] = 0.25f;


        hmm->transitions[2][2] = (1.0-leave) * self_transition; /* self transition */
        hmm->transitions[2][3] = (1.0-leave) * (1.0- self_transition); /* to other state */
        hmm->transitions[2][4] = (1.0-leave) * 0.0; /* to other state */
        hmm->transitions[2][5] = (1.0-leave) * 0.0; /* to other state */
        hmm->transitions[2][ENDSTATE]  = leave;

        hmm->transitions[3][2] = (1.0-leave) * 0.0; /* self transition */
        hmm->transitions[3][3] = (1.0-leave) * self_transition; /* to other state */
        hmm->transitions[3][4] = (1.0-leave) * (1.0- self_transition); /* to other state */
        hmm->transitions[3][5] = (1.0-leave) * 0.0; /* to other state */
        hmm->transitions[3][ENDSTATE]  = leave;

        hmm->transitions[4][2] = (1.0-leave) * 0.0; /* self transition */
        hmm->transitions[4][3] = (1.0-leave) * 0.0; /* to other state */
        hmm->transitions[4][4] = (1.0-leave) * self_transition; /* to other state */
        hmm->transitions[4][5] = (1.0-leave) * (1.0 - self_transition); /* to other state */
        hmm->transitions[4][ENDSTATE]  = leave;

        hmm->transitions[5][2] = (1.0-leave) * (1.0- self_transition); /* self transition */
        hmm->transitions[5][3] = (1.0-leave) * 0.0; /* to other state */
        hmm->transitions[5][4] = (1.0-leave) * 0.0; /* to other state */
        hmm->transitions[5][5] = (1.0-leave) * self_transition; /* to other state */
        hmm->transitions[5][ENDSTATE]  = leave;


        for(i = 0; i < 4;i++){
                for(j = 0; j < 4;j++){
                        if(i ==j){
                                hmm->emissions[i+2][j] = mainres_emission;
                        }else{
                                hmm->emissions[i+2][j] = (1.0 - mainres_emission) / 3.0;
                        }
                        fprintf(stdout,"%f ",hmm->emissions[i+2][j]);
                }
                fprintf(stdout,"\n");
        }

        RUN(sanity_check_hmm(hmm));
        RUN(cumsum_hmm(hmm));


        sum =0.0;
        for(i =0; i < sb->malloc_num;i++){
                //while(sb->seqs[sb->num_seq]->seq_len > expected_length + 50 ||sb->seqs[sb->num_seq]->seq_len <  expected_length - 50 ){
                        RUN(emit_sequence(hmm,sb->seqs[sb->num_seq]));
                        //fprintf(stdout,"len: %d expected: %d\n", sb->seqs[sb->num_seq]->seq_len, expected_length );
                        //}
                sum += (float)(sb->seqs[sb->num_seq]->seq_len -1);
                sb->num_seq++;
        }
        LOG_MSG("Simulated %d seq of length %f\n",sb->malloc_num,sum / (float)sb->malloc_num);
        LOG_MSG("%f leave.",leave);

        /* write sequence  */
        snprintf(buffer, BUFFER_LEN, "%s/%s%0.2f_%s%0.2f.fa",param->outdir,"ACGT_states_RES",mainres_emission,"self",self_transition);
        LOG_MSG("Writing to: %s.",buffer);
        RUN(write_sequences_to_file(sb,buffer));

        free_hmm(hmm);
        return OK;
ERROR:
        free_hmm(hmm);
        return FAIL;

}

int ACGT_embedded_example(struct parameters* param, struct seq_buffer* sb,float mainres_emission)
{
        char buffer[BUFFER_LEN];
        struct hmm* hmm = NULL;
        float leave;
        float sum;
        int i,j;
        int expected_length = 500;


        ASSERT(sb != NULL,"No seq buffer");

        ASSERT(sb->seqs[0]->seq_len == 0,"Need to reset seq buffer");

        ASSERT(mainres_emission <= 1.0,"Main emission has to be smaller than 1.0.");
        ASSERT(mainres_emission >= 0.0,"Main emission has to be greater than 0.0.");


        RUNP(hmm = malloc_hmm(6,4,expected_length));

        /*
          NOTE This should be:
          leave = 1.0 - (float) ((expected_length) /(float)(1+   expected_length));
          BUT:
          I will substract one from the expected length because I co not allow a transition from start to stop (i.e. a zero length sequence
        */

        leave = 1.0 - (float) ((expected_length) /((float) expected_length+2.0));

        /* set transition to end to 1/ expected lengeth  */

        hmm->transitions[STARTSTATE][2] = 0.0f;
        hmm->transitions[STARTSTATE][3] = 0.0f;
        hmm->transitions[STARTSTATE][4] = 0.0f;
        hmm->transitions[STARTSTATE][5] = 0.0f;
        hmm->transitions[STARTSTATE][6] = 1.0f;
        hmm->transitions[STARTSTATE][7] = 0.0f;



        hmm->transitions[2][2] = 0.0; /* self transition */
        hmm->transitions[2][3] = 1.0; /* to other state */
        hmm->transitions[2][4] = 0.0; /* to other state */
        hmm->transitions[2][5] = 0.0; /* to other state */
        hmm->transitions[2][6] = 0.0; /* to other state */
        hmm->transitions[2][7] = 0.0;
        hmm->transitions[2][ENDSTATE]  = 0.0;

        hmm->transitions[3][2] = 0.0; /* self transition */
        hmm->transitions[3][3] = 0.0; /* to other state */
        hmm->transitions[3][4] = 1.0; /* to other state */
        hmm->transitions[3][5] = 0.0; /* to other state */
        hmm->transitions[3][6] = 0.0; /* to other state */
        hmm->transitions[3][7] = 0.0;
        hmm->transitions[3][ENDSTATE]  = 0.0;

        hmm->transitions[4][2] = 0.0; /* self transition */
        hmm->transitions[4][3] = 0.0; /* to other state */
        hmm->transitions[4][4] = 0.0; /* to other state */
        hmm->transitions[4][5] = 1.0; /* to other state */
        hmm->transitions[4][6] = 0.0; /* to other state */
        hmm->transitions[4][7] = 0.0;
        hmm->transitions[4][ENDSTATE]  = 0.0;

        hmm->transitions[5][2] = 0.0; /* self transition */
        hmm->transitions[5][3] = 0.0; /* to other state */
        hmm->transitions[5][4] = 0.0; /* to other state */
        hmm->transitions[5][5] = 0.0; /* to other state */
        hmm->transitions[5][6] = 0.0; /* to other state */
        hmm->transitions[5][7]  = 1.0;
        hmm->transitions[5][ENDSTATE]  = 0.0;

        hmm->transitions[6][2] = leave; /* self transition */
        hmm->transitions[6][3] = 0.0; /* to other state */
        hmm->transitions[6][4] = 0.0; /* to other state */
        hmm->transitions[6][5] = 0.0; /* to other state */
        hmm->transitions[6][6] = 1.0 - leave; /* to other state */
        hmm->transitions[6][7] = 0.0; /* to other state */
        hmm->transitions[6][ENDSTATE]  = 0.0;

        hmm->transitions[7][2] = 0.0; /* self transition */
        hmm->transitions[7][3] = 0.0; /* to other state */
        hmm->transitions[7][4] = 0.0; /* to other state */
        hmm->transitions[7][5] = 0.0; /* to other state */
        hmm->transitions[7][6] = (1.0 - leave)  * 0.0f; /* to other 5' loop state */
        hmm->transitions[7][7] = (1.0 - leave)  * 1.0f; /* to other state */
        hmm->transitions[7][ENDSTATE]  = leave;

        for(i = 0; i < 4;i++){
                hmm->emissions[6][i] = 0.25;
                hmm->emissions[7][i] = 0.25;
        }


        for(i = 0; i < 4;i++){
                for(j = 0; j < 4;j++){
                        if(i ==j){
                                hmm->emissions[i+2][j] = mainres_emission;
                        }else{
                                hmm->emissions[i+2][j] = (1.0 - mainres_emission) / 3.0;
                        }
                        fprintf(stdout,"%f ",hmm->emissions[i+2][j]);
                }
                fprintf(stdout,"\n");
        }

        RUN(sanity_check_hmm(hmm));
        RUN(cumsum_hmm(hmm));

        sum =0.0;
        for(i =0; i < sb->malloc_num;i++){
                RUN(emit_sequence(hmm,sb->seqs[sb->num_seq]));
                sum += (float)(sb->seqs[sb->num_seq]->seq_len -1);
                sb->num_seq++;
        }
        LOG_MSG("Simulated %d seq of length %f\n",sb->malloc_num,sum / (float)sb->malloc_num);
        LOG_MSG("%f leave.",leave);

        /* write sequence  */
        snprintf(buffer, BUFFER_LEN, "%s/%s%0.2f.fa",param->outdir,"ACGT_embedded_RES",mainres_emission);
        LOG_MSG("Writing to: %s.",buffer);
        RUN(write_sequences_to_file(sb,buffer));

        free_hmm(hmm);
        return OK;
ERROR:
        free_hmm(hmm);
        return FAIL;

}




int two_state_example(struct parameters* param,  struct seq_buffer* sb)
{
        char buffer[BUFFER_LEN];
        struct hmm* hmm = NULL;
        float leave;
        float sum;
        int i;
        int expected_length = 1000;

        ASSERT(sb != NULL,"No seq buffer");

        RUNP(hmm = malloc_hmm(2,4,expected_length));

        /*
          NOTE This should be:
          leave = 1.0 - (float) ((expected_length) /(float)(1+   expected_length));
          BUT:
          I will substract one from the expected length because I co not allow a transition from start to stop (i.e. a zero length sequence
        */

        leave = 1.0 - (float) ((expected_length-1) /(float)(   expected_length));

        /* set transition to end to 1/ expected lengeth  */

        hmm->transitions[STARTSTATE][2] = 0.5f;
        hmm->transitions[STARTSTATE][3] = 0.5f;

        hmm->transitions[2][2] = (1.0-leave) * 0.1; /* self transition */
        hmm->transitions[2][3] = (1.0-leave) * 0.9; /* to other state */
        hmm->transitions[2][ENDSTATE]  = leave;

        hmm->transitions[3][2] = (1.0-leave) * 0.9; /* to other state */
        hmm->transitions[3][3] = (1.0-leave) * 0.1; /* self transition */
        hmm->transitions[3][ENDSTATE]  = leave;


        hmm->emissions[2][0] = 0.7;
        hmm->emissions[2][1] = 0.1;
        hmm->emissions[2][2] = 0.1;
        hmm->emissions[2][3] = 0.1;


        hmm->emissions[3][0] = 0.1;
        hmm->emissions[3][1] =0.1;
        hmm->emissions[3][2] = 0.1;
        hmm->emissions[3][3] = 0.7;

        RUN(sanity_check_hmm(hmm));
        RUN(cumsum_hmm(hmm));

        sum =0.0;
        for(i =0; i < sb->malloc_num;i++){
                RUN(emit_sequence(hmm,sb->seqs[sb->num_seq]));
                sum += (float)(sb->seqs[sb->num_seq]->seq_len -1);
                sb->num_seq++;
        }
        LOG_MSG("Simulated %d seq of length %f\n",sb->malloc_num,sum / (float)sb->malloc_num);
        LOG_MSG("%f leave.",leave);

        /* write sequence  */
        snprintf(buffer, BUFFER_LEN, "%s/%s",param->outdir,"two_states.fa");
        LOG_MSG("Writing to: %s.",buffer);
        RUN(write_sequences_to_file(sb,buffer));

        free_hmm(hmm);
        return OK;
ERROR:
        free_hmm(hmm);
        return FAIL;

}


int emit_sequence(struct hmm* hmm, struct sequence* sequence)
{

        int state;
        float r;
        //int sum = 0;
        int i;
        ASSERT(sequence != NULL, "No sequence");
        sequence->seq_len = 0;
        state = STARTSTATE;
        while(state != ENDSTATE){
                /* transition */
                r = random_float_zero_to_x(1.0);
                for(i = 0 ; i < hmm->num_states;i++){
                        if(r <= hmm->transitions[state][i]){
                                state = i;
                                break;
                        }
                }
                /* emission */
                r = random_float_zero_to_x(1.0);
                for(i = 0 ; i < hmm->alphabet_len;i++){
                        if(r <= hmm->emissions[state][i]){
                                sequence->seq[sequence->seq_len] = "ACGT"[i];
                                sequence->seq_len++;
                                if(sequence->seq_len == sequence->malloc_len){
                                        sequence->malloc_len = sequence->malloc_len << 1;
                                        MREALLOC(sequence->seq, sizeof(uint8_t) *sequence->malloc_len);
                                }
                                break;
                        }
                }
                //fprintf(stdout,"%d %c\n",state,(char) sequence->seq[sequence->seq_len-1]);
        }
        sequence->seq[sequence->seq_len] = 0;
        sequence->seq_len++;

        if(sequence->seq_len == sequence->malloc_len){
                sequence->malloc_len = sequence->malloc_len << 1;
                MREALLOC(sequence->seq, sizeof(uint8_t) *sequence->malloc_len);
        }

        //fprintf(stdout,"%d\t%s\n",sequence->seq_len-1, sequence->seq);
        return OK;
ERROR:
        return FAIL;
}

int cumsum_hmm(struct hmm* hmm)
{
        int i,j;
         for(i = 0; i < hmm->num_states;i++){
                for(j = 1; j < hmm->num_states;j++){
                        hmm->transitions[i][j]+=hmm->transitions[i][j-1];

                }
        }

        for(i = 2; i < hmm->num_states;i++){

                for(j = 1; j < 4;j++){
                        hmm->emissions[i][j]+=hmm->emissions[i][j-1];
                }

        }
        return OK;
}

struct seq_buffer* alloc_seq_buffer(int num_seq)
{
        struct seq_buffer* sb  = NULL;
        struct sequence* sequence = NULL;
        int i;

        ASSERT(num_seq > 0, "No parameters.");

        MMALLOC(sb,sizeof(struct seq_buffer));
        sb->malloc_num = num_seq;
        sb->num_seq = 0;
        sb->seqs = NULL;
        MMALLOC(sb->seqs, sizeof(struct chromosome*) *sb->malloc_num );
        for(i = 0; i < sb->malloc_num;i++){
                sb->seqs[i] = NULL;
                sequence = NULL;
                MMALLOC(sequence,sizeof(struct sequence));
                sequence->seq = NULL;
                sequence->name = NULL;
                sequence->malloc_len = 128;
                sequence->seq_len = 0;
                MMALLOC(sequence->seq, sizeof(uint8_t) * sequence->malloc_len);
                MMALLOC(sequence->name, sizeof(char) * 256);
                sb->seqs[i] = sequence;
        }
        return sb;

        /*sequence->seq[sequence->seq_len] = line[i];
        sequence->seq_len++;
        if(sequence->seq_len == sequence->malloc_len){
                sequence->malloc_len = sequence->malloc_len << 1;
                MREALLOC(sequence->seq, sizeof(uint8_t) *sequence->malloc_len);
                }*/

ERROR:
        free_sb(sb);
        return NULL;
}

int write_sequences_to_file(struct seq_buffer* sb,char* filename)
{
        FILE* f_ptr = NULL;
        int i;
        ASSERT(sb != NULL," no sequences.");
        ASSERT(filename != NULL ," No filename.");
        RUNP(f_ptr = fopen(filename, "w"));

        for(i = 0; i < sb->num_seq;i++){
                fprintf(f_ptr,">Seq_%d\n%s\n",i+1,sb->seqs[i]->seq);
        }
        fclose(f_ptr);
        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        return FAIL;

}

int reset_sb(struct seq_buffer* sb)
{
        int i;
        sb->num_seq = 0;
        for(i = 0; i < sb->malloc_num;i++){
                sb->seqs[i]->seq_len = 0;
        }
        return OK;
}

void free_sb(struct seq_buffer* sb)
{
        int i;
        struct sequence* seq = NULL;
        if(sb){
                for(i =0; i < sb->malloc_num;i++){
                        seq = sb->seqs[i];
                        MFREE(seq->seq);
                        MFREE(seq->name);
                        MFREE(seq);
                }
                MFREE(sb->seqs);
                MFREE(sb);
        }
}


struct hmm* malloc_hmm(int num_states, int alphabet_len, int max_seq_len)
{
        struct hmm* hmm = NULL;


        MMALLOC(hmm, sizeof(struct hmm));
        hmm->alphabet_len =  alphabet_len;
        hmm->num_states = num_states + NUM_ADDITIONAL_STATES;
        hmm->max_seq_len = max_seq_len;
        hmm->emissions = NULL;
        hmm->transitions = NULL;

        hmm->background = NULL;

        MCALLOC(hmm->background, hmm->alphabet_len,sizeof(float));

	        hmm->emissions = malloc_2d_float(hmm->emissions, hmm->num_states, hmm->alphabet_len, 0.0f);
          hmm->transitions = malloc_2d_float(hmm->transitions, hmm->num_states, hmm->num_states,  0.0f);

          return hmm;
ERROR:
          free_hmm(hmm);
          return NULL;
}

void free_hmm(struct hmm* hmm)
{

        if(hmm){
                free_2d((void**)hmm->emissions );
                free_2d((void**)hmm->transitions );

                MFREE(hmm->background);
                MFREE(hmm);
        }
}

int sanity_check_hmm(struct hmm*  hmm)
{
        int i,j;

        init_logsum();
        float sum = 0.0;


        for(i = 0; i < hmm->num_states;i++){
                if(i != ENDSTATE){
                        sum = 0.0;
                        for(j = 0; j < hmm->num_states;j++){
                                sum +=  hmm->transitions[i][j];
                        }
                        ASSERT(fabs(1.0 - sum) < FLT_EPSILON,"State transitions from %d do not sum to 1:%20.20f",i,sum);
                }
        }

        for(i = 2; i < hmm->num_states;i++){
                sum = 0.0;
                for(j = 0; j < 4;j++){
                        sum += hmm->emissions[i][j];
                }
                ASSERT(fabs(1.0 - sum) < FLT_EPSILON,"State transitions from %d do not sum to 1:%20.20f",i,sum);
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
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--num","Number of sequences." ,"[1000]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--len","Average length of sequences." ,"[1000]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--out","Output directory path." ,"[1000]"  );
        return OK;
}

