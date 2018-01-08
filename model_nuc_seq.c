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

#include "ihmm.h"

#include "ihmm_io.h"

#include "matrix_io.h"

struct parameters{       
        char* input;
        char* outdir;
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

int make_model(struct parameters* param, struct seq_buffer* sb);

int make_dot_files(struct parameters* param);


static int print_help(char **argv);
static int free_parameters(struct parameters* param);


int main (int argc, char *argv[]) 
{		
        struct parameters* param = NULL;
        int c;
        
        tlog.echo_build_config();
        
        MMALLOC(param, sizeof(struct parameters));
        param->input = NULL;
        param->outdir = NULL;
                
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

        if(!param->input){
                RUN(print_help(argv));
                ERROR_MSG("No input file! use --in <blah.fa>");
                
        }else{
                if(!my_file_exists(param->input)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->input);               
                }           
        }
        if(!param->outdir){
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
        char buffer[BUFFER_LEN];
        struct seq_buffer* sb = NULL;
        
        int i;
       
        ASSERT(param!= NULL, "No parameters found.");
        


        RUN(create_output_directories(param->outdir));
        
        RUN(set_log_file(param->outdir,"scs_net"));

        snprintf(buffer, BUFFER_LEN, "%s/%s/",param->outdir,OUTDIR_CHECKPOINTS);
        DECLARE_CHK(MAIN_CHECK, buffer);
        
        /* Step one read in sequences */
        LOG_MSG("Loading sequences.");

       


        RUNP(sb = load_sequences(param));

        for(i = 0; i < MACRO_MIN(sb->num_seq, 10);i++){
                fprintf(stdout,"%i\t%s\n",i,sb->seqs[i]->name);
                fprintf(stdout,"%i\t%s\n",i,(char*) sb->seqs[i]->seq);

        }
        

        LOG_MSG("Read %d sequences.",sb->num_seq);

      
       
	
        LOG_MSG("Make model.");
        snprintf(buffer, BUFFER_LEN, "make model out of %s.",param->input);
        RUN_CHECKPOINT(MAIN_CHECK,make_model(param,sb),buffer);

        LOG_MSG("Make dotfiles.");
        snprintf(buffer, BUFFER_LEN, "make model out of %s.",param->input);
        //RUN_CHECKPOINT(MAIN_CHECK,make_dot_files(param),buffer);
        RUN(make_dot_files(param));
        
        /*for(i = 0; i <= model->K;i++){
          model->sumM[i] = 0;
          }
          for(i = 0; i < iseq->num_seq;i++){
          for(j= 0; j <  iseq->len[i];j++){
          model->sumM[iseq->labels[i][j]]++;
          }
          }*/
	
        //print_transistion_matrix(model);
        //print_emisson_matrix(model);
	
        /*for(i =0;i < iseq->num_seq;i++){
          for(j = 0; j < iseq->len[i];j++){
          fprintf(stdout,"  %c ","ACGT"[(int)iseq->seq[i][j]]);
          }
          fprintf(stdout,"\n");
          for(j = 0; j < iseq->len[i];j++){
          fprintf(stdout,"%3d ",iseq->labels[i][j]);
          }
          fprintf(stdout,"\n");
          }*/
	
	
        //DPRINTF2("DONE PGAS.");
       
        
       
        free_sb(sb);
        DESTROY_CHK(MAIN_CHECK);
        
        return OK;
ERROR:
        free_sb(sb);
        return FAIL;
}





int make_model(struct parameters* param,struct seq_buffer* sb)
{
        char buffer[BUFFER_LEN];
        struct iHMM_model* model = NULL;
        char** tmp_seq_pointer = NULL;
        
        int i;
        ASSERT(param != NULL, "No param.");

        MMALLOC(tmp_seq_pointer,sizeof(char*) * sb->num_seq);
        for(i = 0; i < sb->num_seq;i++){
                tmp_seq_pointer[i] = (char*) sb->seqs[i]->seq;
        }
        
        RUNP(model = init_iHMM_model());
        /* Don't think I need three different variables here...  */
        model->numb = 100;
        model->nums = 1;
        model->numi = 1;

        //DPRINTF2("START PGAS.");
        RUN(particle_gibbs_with_ancestors_controller(model, tmp_seq_pointer,sb->num_seq));
        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",param->outdir,OUTDIR_MODEL,"iHMM_model_parameters.csv");
        LOG_MSG("Writing to: %s.",buffer);
        RUN(write_ihmm_parameters_to_file(model,buffer));
        free_iHMM_model(model);
        MFREE(tmp_seq_pointer);
        return OK;
ERROR:
        free_iHMM_model(model);
        MFREE(tmp_seq_pointer);
        return FAIL;
}


int make_dot_files(struct parameters* param)
{
        char buffer[BUFFER_LEN];
        struct double_matrix* matrix = NULL;
        FILE* f_ptr = NULL;

        double tmp_sum[4];

        double background[4];
        double IC;
        double max_IC;
        double tmp_prob;
        int ncol;
        int nrow;
        int i,j;

        int max_stack_height = 128;
        
        
        ASSERT(param != NULL, "No param.");

        
        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",param->outdir,OUTDIR_MODEL,"iHMM_model_parameters.csv");
        LOG_MSG("Reading from:  %s.",buffer);

        RUNP(matrix = read_double_matrix(buffer,1,1));
        print_double_matrix(matrix,stdout,1,1);
        ncol = matrix->ncol;
        nrow = matrix->nrow;

        for(i = 0; i < 4;i++){
                background[i] = scaledprob2prob( matrix->matrix[i][ncol-1]);
                fprintf(stdout,"%f\n",background[i]);
        }

       

        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",param->outdir,OUTDIR_VIZ,"test_dotfile.dot");
        RUNP(f_ptr = fopen(buffer, "w"));

        /* print dot header...  */

        fprintf(f_ptr,"digraph structs {\n");
        fprintf(f_ptr,"rankdir=LR;\n");
        fprintf(f_ptr,"overlap=false;\n");
        fprintf(f_ptr,"node [shape=circle];\n");//plaintext shape?
        max_IC = -1000.0;
        for(i = 0;i< 4;i++){
                IC = 0.0;
                
                for(j = 0; j < 4;j++){
                        if(j ==i){
                                tmp_prob = 1.0;
                        }else{
                                tmp_prob = 1e-7;//0.0;    
                        }
                        IC += tmp_prob * log2( tmp_prob / background[j]);
                        fprintf(stdout,"%f\t%f\t%f\n",tmp_prob  ,  background[j],log2(tmp_prob / background[j]));
                }

                fprintf(stdout,"%f\n",IC);
                if(IC> max_IC){
                        max_IC = IC;
                }
        }

        /* print start stop nodes  */
         
        fprintf(f_ptr,"%s [label=Start]\n", matrix->col_names[0]);
        fprintf(f_ptr,"%s [label=End]\n", matrix->col_names[1]);
         
         
        
        /* print nodes...  */
        for(i = 2;i < ncol-1;i++){
                IC = 0.0;

                for(j = 0; j < 4;j++){
                        IC += matrix->matrix[j][i] * log2( matrix->matrix[j][i] / background[j]);
                        fprintf(stdout,"%f\t%f\t%f\n", matrix->matrix[j][i],  background[j],log2( matrix->matrix[j][i] / background[j]));
                }
                fprintf(stdout,"%f\n",IC);
                tmp_prob = (double) max_stack_height / max_IC *IC;
                for(j = 0; j < 4;j++){
                        
                        tmp_sum[j] =  matrix->matrix[j][i] * tmp_prob ;
                        fprintf(stdout,"%f\n",tmp_sum[j]);
                }

                
                fprintf(f_ptr,"%s [label=<\n", matrix->col_names[i]);
                fprintf(f_ptr,"<TABLE CELLPADDING=\"0\" BORDER=\"0\" CELLSPACING=\"0\">\n");
                fprintf(f_ptr,"<TR>\n");
                fprintf(f_ptr,"<TD BGCOLOR=\"gray\"><FONT POINT-SIZE=\"%d\" COLOR=\"#8470ff\">A</FONT></TD>\n",(int)tmp_sum[0]);
                fprintf(f_ptr,"</TR>\n");
                fprintf(f_ptr,"<TR>\n");
                fprintf(f_ptr,"<TD BGCOLOR=\"gray\"><FONT POINT-SIZE=\"%d\"  COLOR=\"#f4a460\">C</FONT></TD>\n",(int)tmp_sum[1]);
                fprintf(f_ptr,"</TR>\n");
                fprintf(f_ptr,"<TR>\n");
                fprintf(f_ptr,"<TD BGCOLOR=\"gray\"><FONT POINT-SIZE=\"%d\" COLOR=\"#f08080\">G</FONT></TD>\n",(int)tmp_sum[2]);
                fprintf(f_ptr,"</TR>\n");
                fprintf(f_ptr,"<TR>\n");
                fprintf(f_ptr,"<TD BGCOLOR=\"gray\"><FONT POINT-SIZE=\"%d\" COLOR=\"#90ee90\">T</FONT></TD>\n",(int)tmp_sum[3]);
                fprintf(f_ptr,"</TR>\n");
                fprintf(f_ptr,"</TABLE>>];\n");
                
        }

        fprintf(f_ptr,"\n\n");
        /* print edges */
     
        for(i = 0; i < ncol-1;i++){

                for(j = 0;j < ncol-1;j++){
                        if(matrix->matrix[i+6][j] >= 1e-2){
                                fprintf(f_ptr,"%s -> %s[label=\"%0.2f\"];\n",matrix->col_names[j],matrix->col_names[i],matrix->matrix[i+6][j]);
                        }
                       
                }
        }
        
        /* print end of dot file  */

        fprintf(f_ptr,"}\n");


        fclose(f_ptr);

        LOG_MSG("To visualize: dot  -Tpdf  <.dot file>  -o  <blah.pdf>.");

        free_double_matrix(matrix);
        return OK;
ERROR:
        free_double_matrix(matrix);
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


        seq_p = 0;
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
        fclose(f_ptr);
        return sb;
ERROR:
        free_sb(sb);
        if(f_ptr){
                fclose(f_ptr);
        }
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
