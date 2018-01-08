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

struct parameters{       
        char* input;
        char* output;
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


static int run_build_ihmm(struct parameters* param);

static struct seq_buffer* load_sequences(struct parameters* param);
static void free_sb(struct seq_buffer* sb);

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
                
        while (1){	
                static struct option long_options[] ={
                        {"in",required_argument,0,'i'},
                        {"out",required_argument,0,'o'},			       
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"hi:o:",long_options, &option_index);
		
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
                ERROR_MSG("No output file! use -o <blah.fa>");
                
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
        ASSERT(param!= NULL, "No parameters found.");
        
        struct seq_buffer* sb = NULL;
        int i;

        /* Step one read in sequences */
        LOG_MSG("Loading sequences.");
        
        RUNP(sb = load_sequences(param));

        for(i = 0; i < MACRO_MIN(sb->num_seq, 10);i++){
                fprintf(stdout,"%i\t%s\n",i,sb->seqs[i]->name);
                fprintf(stdout,"%i\t%s\n",i,(char*) sb->seqs[i]->seq);

        }
        free_sb(sb);
                
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
//        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--seed-len","Length of seeds." ,"[12]"  );
        //    	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--seed-step","Distance between seeds." ,"[8]"  );
//      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthread","Number of threads." ,"[8]"  );
//      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--output","Output file name." ,"[?]"  );
        return OK;
}

struct seq_buffer* load_sequences(struct parameters* param)
{
       
        struct seq_buffer* sb  = NULL;
        struct sequence* sequence = NULL;
        FILE* f_ptr = NULL;
        char line[LINE_LEN];
        int i, seq_p;
        ASSERT(param != NULL, "No parameters.");

        ASSERT(param->input != NULL,"No input file specified - this should have been caught before!");

        MMALLOC(sb,sizeof(struct seq_buffer));
        sb->malloc_num = 1024;
        sb->num_seq = -1; 
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
        
        RUNP(f_ptr = fopen(param->input, "r" ));
        while(fgets(line, LINE_LEN, f_ptr)){
                //reading name .... 
                //fprintf(stderr,"%s",line);
                //if(line[0] != '#'){ 
                if(line[0] == '@' || line[0] == '>'){
                        line[strlen(line)-1] = 0;                       
                        for(i =0 ; i<strlen(line);i++){
                                if(isspace(line[i])){
                                        line[i] = 0;
                                                
                                }
                                        
                        }
                        sb->num_seq++; 

                        if(sb->num_seq == sb->malloc_num){
                                sb->malloc_num = sb->malloc_num << 1;
                                MREALLOC(sb->seqs,sizeof(struct sequence*) * sb->malloc_num);
                                for(i = sb->num_seq; i < sb->malloc_num;i++){
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
                                
                        }
                        sequence = sb->seqs[sb->num_seq];
                        snprintf(sequence->name,256,"%s",line+1);

                        seq_p = 1;

                }else if(line[0] == '+'){
                        seq_p = 0;

                }else{	
                
                        if(seq_p){

                                for(i = 0;i < LINE_LEN;i++){
                                        if(isalpha((int)line[i])){
                                                sequence->seq[sequence->seq_len] = line[i];
                                                sequence->seq_len++;
                                                if(sequence->seq_len == sequence->malloc_len){
                                                        sequence->malloc_len = sequence->malloc_len << 1;
                                                        MREALLOC(sequence->seq, sizeof(uint8_t) *sequence->malloc_len);
                                                }
                                        }
                                        if(iscntrl((int)line[i])){
                                                sequence->seq[sequence->seq_len] = 0;
                                                break;
                                        }
                                }
                        } /* here I would look for quality values.... */
                      
                }
        }

        return sb;
ERROR:
        return NULL;
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


/*
struct qs_struct* read_fasta_fastq2(FILE *file, struct qs_struct* qs,struct parameters* param) 
{
        int park_pos = -1;
        char line[LINE_LEN];
        unsigned char rc[LINE_LEN];
        int i,j;
        int seq_p = 0;
        int set = 0;
        int len = 0;
        int fastq = 0;
        qs->size = 0;
        for(i = 0; i < param->num_query ;i++){
                //qs->seq_info[i] = malloc(sizeof(struct seq_info));
                if(qs->seq_info[i]->sn){
                        free(qs->seq_info[i]->sn);//
                        free(qs->seq_info[i]->seq);
                        free(qs->seq_info[i]->reverse_seq);
                }
                if(qs->seq_info[i]->qual){
                        //fprintf(stderr,"%s\n",qs->seq_info[i]->qual);
                        free(qs->seq_info[i]->qual);
        s        }

                qs->seq_info[i]->sn = 0;
                qs->seq_info[i]->seq = 0;
                qs->seq_info[i]->reverse_seq = 0;
                qs->seq_info[i]->qual = 0;
		
		
                qs->seq_info[i]->last_f_test = -1000;
                qs->seq_info[i]->last_r_test = -1000;
                //qs->seq_info[i]->match_units= malloc(sizeof(unsigned long int) *  LIST_STORE_SIZE);//param->reported_seed_hits );
                for(j = 0; j < LIST_STORE_SIZE;j++){
                        qs->seq_info[i]->hits[j].chr = 0;
                        qs->seq_info[i]->hits[j].pos = 0;
                        qs->seq_info[i]->hits[j].strand = 0;
                        qs->seq_info[i]->hits[j].error = 0xF;
                        
                }
        }
                //}
        }
        return qs;
}

*/
