
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

#include "outdir.h"

#include "matrix_io.h"

#include "make_dot_file.h"

struct parameters{       
        char* input;
        char* outdir;
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


int two_state_example(struct parameters* param, struct seq_buffer* sb);

static int print_help(char **argv);
static int free_parameters(struct parameters* param);






/* seqbuffer stuff  */
static struct seq_buffer* alloc_seq_buffer(struct parameters* param);
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
        param->outdir = NULL;
                
        while (1){	
                static struct option long_options[] ={
                        {"num",required_argument,0,'n'},
                        {"out",required_argument,0,'o'},			       
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"hn:o:",long_options, &option_index);
		
                if (c == -1){
                        break;
                }
                switch(c) {
                case 'n':
                        param->num_seq = atoi(optarg);
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
        
        int i;
       
        ASSERT(param!= NULL, "No parameters found.");
        


        RUN(create_output_directories(param->outdir));
        
        RUN(set_log_file(param->outdir,"scs_net"));

        snprintf(buffer, BUFFER_LEN, "%s/%s/",param->outdir,OUTDIR_CHECKPOINTS);
        DECLARE_CHK(MAIN_CHECK, buffer);
        
        /* Step one allocate seq struct.. */
        RUNP(sb = alloc_seq_buffer(param));


        /* make sequence test sets; each should have a fasta file, the model
         * parameter matrix and a dot file */

        /* CpG example from black book */

        RUN(two_state_example(param,sb));
        
        RUN(reset_sb(sb));
        
        free_sb(sb);
        DESTROY_CHK(MAIN_CHECK);
        
        return OK;
ERROR:
        free_sb(sb);
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
        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",param->outdir,OUTDIR_MODEL,"two_states.fa");
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
        int sum = 0;
        int i;

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

struct seq_buffer* alloc_seq_buffer(struct parameters* param)
{
        struct seq_buffer* sb  = NULL;
        struct sequence* sequence = NULL;
        int i;
        
        ASSERT(param != NULL, "No parameters.");
        
        MMALLOC(sb,sizeof(struct seq_buffer));
        sb->malloc_num = param->num_seq;
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
                sum = (int) (10000.0 * sum);
                
                ASSERT(sum == 10000.0,"State transitions from %d do not sum to 1:%20.20f",i,sum);
                }
        }
	
        for(i = 2; i < hmm->num_states;i++){
                sum = 0.0;
                for(j = 0; j < 4;j++){
                        sum += hmm->emissions[i][j];
                }
		
                sum = (int) (10000.0 * sum);
                ASSERT(sum == 10000.0 ,"State transitions from %d do not sum to 1:%20.20f",i,sum);
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
        return OK;
}

