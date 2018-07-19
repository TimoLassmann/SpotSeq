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


struct parameters{       
        char* input;
        char* output;
};


static int run_plot_ihmm(struct parameters* param);

static int make_dot_file(struct fast_hmm_param* ft, struct ihmm_model* model, struct parameters* param);

static int print_help(char **argv);
static int free_parameters(struct parameters* param);

int get_color(char* color, float x, float min_x, float max_x);

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
                ERROR_MSG("No output file! use --out <blah.fa>");
        }else{
                if(my_file_exists(param->output)){
                        WARNING_MSG("Will overwrite: %s.",param->output);          
                }   
        }

        
        RUN(run_plot_ihmm(param));

        //RUN(seed_controller_thread(param));
        
        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}

int run_plot_ihmm(struct parameters* param)
{
        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;
        
        int initial_states = 10;
        ASSERT(param!= NULL, "No parameters found.");
        

        
        RUNP(model = read_model(param->input));
        print_model_parameters(model);
        print_counts(model);
        RUNP(ft = alloc_fast_hmm_param(initial_states,model->L));
        RUN(fill_background_emission_from_model(ft,model));

        RUN(fill_fast_transitions(model,ft));
        //RUN(print_fast_hmm_params(ft));
        RUN(make_dot_file( ft, model, param));
                
        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        return OK;
ERROR:
        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        return FAIL;
}

int make_dot_file(struct fast_hmm_param* ft, struct ihmm_model* model, struct parameters* param)
{
        FILE* f_ptr = NULL;

        float* tmp_sum = NULL;
        int* total_counts = NULL;

        float* background;
        char color_buffer[BUFFER_LEN];
        double IC;
        double max_IC;
        double tmp_prob;
        double sum_usage; 
        int i,j;

        int max_stack_height = 128;

        ASSERT(ft != NULL, "No input matrix.");
     
        MMALLOC(tmp_sum, sizeof(float) * ft->L);
        MMALLOC(total_counts, sizeof(int) * model->num_states);
        background = ft->background_emission;
     
        
        RUNP(f_ptr = fopen(param->output, "w"));
        
        /* print dot header...  */
        
        fprintf(f_ptr,"digraph structs {\n");
        fprintf(f_ptr,"rankdir=LR;\n");
        fprintf(f_ptr,"overlap=false;\n");
        fprintf(f_ptr,"node [shape=circle];\n");//plaintext shape?
        max_IC = -1000.0;
        for(i = 0;i< 4;i++){
                IC = 0.0;
                
                for(j = 0; j < ft->L;j++){
                        if(j ==i){
                                tmp_prob = 1.0;
                        }else{
                                tmp_prob = 1e-7;//0.0;    
                        }
                        IC += tmp_prob * log2( tmp_prob / background[j]);
                        //       fprintf(stdout,"%f\t%f\t%f\n",tmp_prob  ,  background[j],log2(tmp_prob / background[j]));
                }

                //fprintf(stdout,"%f\n",IC);
                if(IC> max_IC){
                        max_IC = IC;
                }
        }

        /* print start stop nodes  */

        fprintf(f_ptr,"State%d [label=Start]\n", 0);
        fprintf(f_ptr,"State%d [label=End]\n", 1);

        /* I think here I try to figure out which state is used most -> has the
         * highest number of emissions. */
        sum_usage = 0;
        for(i = 0; i < model->num_states;i++){
                total_counts[i] = 0;
        }
        
        for(i = 0; i < model->L;i++){
                for(j = 0; j < model->num_states;j++){
                        total_counts[j] += model->emission_counts[i][j];
                }
        }
        for(j = 2; j < model->num_states;j++){
                if(total_counts[j] >  sum_usage){
                        sum_usage = total_counts[j];
                }
                fprintf(stdout,"%d: %d\n",j,  total_counts[j]);
        }
        fprintf(stdout,"%f\n",sum_usage);

        
        
        /* print nodes...  */
        for(i = 2; i < model->num_states;i++){
                IC = 0.0;

                for(j = 0; j < ft->L;j++){
                        IC += ft->emission[j][i] * log2(ft->emission[j][i] / background[j]);
                        
                        //                    fprintf(stdout,"%f\t%f\t%f\n", matrix->matrix[j][i],  background[j],log2( matrix->matrix[j][i] / background[j]));
                }
//                fprintf(stdout,"%f\n",IC);
                /* scale total bar height by information content */
                tmp_prob =  (double) max_stack_height;// * ((double) total_counts[i]) / (double)sum_usage);// / model->emission_counts[0][i];// *((double) model->emission_counts[0][i]  /sum_usage);// max_IC *IC      ;

                /* Further scale bar height by usage of state.  */
                //tmp_prob = tmp_prob * (double) matrix->matrix[4][i] / sum_usage;
                fprintf(stdout,"%d %f %f\n",i,tmp_prob ,(double) model->emission_counts[0][i] / sum_usage);
                LOG_MSG("Len:%d %d",ft->L,model->num_states);

                //ncol =0;  
                for(j = 0; j < ft->L;j++){
                        
                        tmp_sum[j] =  ft->emission[j][i] * tmp_prob  ;
                        fprintf(stdout,"%d %c %f\n", i,"ACDEFGHIKLMNPQRSTVWY"[j],ft->emission[j][i]  );
                        //LOG_MSG("i:%d j:%d  %d  (wetfg%d)",i,j,ft->L, ncol);
                        //ncol++;
                        
                        //       fprintf(stdout,"\t%f afa\n",tmp_sum[j]);
                }
                fprintf(f_ptr,"State%d [label=<\n",i);
                
                fprintf(f_ptr,"<TABLE CELLPADDING=\"0\" BORDER=\"0\" CELLSPACING=\"0\">\n");

                fprintf(f_ptr,"<TR>\n");
                                fprintf(f_ptr,"<TD BGCOLOR=\"gray\"><FONT POINT-SIZE=\"%d\">%d</FONT></TD>\n",12,total_counts[i]);
                                fprintf(f_ptr,"</TR>\n");

                
                if(ft->L == ALPHABET_DNA){
                        for(j = 0; j < ft->L;j++){
                                fprintf(f_ptr,"<TR>\n");
                                fprintf(f_ptr,"<TD BGCOLOR=\"gray\"><FONT POINT-SIZE=\"%d\">%c</FONT></TD>\n",MACRO_MAX( (int)tmp_sum[j],1),"ACGTN"[j]);
                                fprintf(f_ptr,"</TR>\n");
                        }
                }else if(ft->L == ALPHABET_PROTEIN ){
                        for(j = 0; j < ft->L;j++){
                                fprintf(f_ptr,"<TR>\n");
                                fprintf(f_ptr,"<TD BGCOLOR=\"gray\"><FONT POINT-SIZE=\"%d\">%c</FONT></TD>\n",MACRO_MAX ((int)tmp_sum[j],1),"ACDEFGHIKLMNPQRSTVWY"[j]);
                                fprintf(f_ptr,"</TR>\n");
                        }
                }
             
                
                
                /*
                
                
                fprintf(f_ptr,"<TR>\n");
                fprintf(f_ptr,"<TD BGCOLOR=\"gray\"><FONT POINT-SIZE=\"%d\"  COLOR=\"#f4a460\">C</FONT></TD>\n",(int)tmp_sum[1]);
                fprintf(f_ptr,"</TR>\n");
                
                fprintf(f_ptr,"<TR>\n");
                fprintf(f_ptr,"<TD BGCOLOR=\"gray\"><FONT POINT-SIZE=\"%d\" COLOR=\"#f08080\">G</FONT></TD>\n",(int)tmp_sum[2]);
                fprintf(f_ptr,"</TR>\n");
                
                fprintf(f_ptr,"<TR>\n");
                fprintf(f_ptr,"<TD BGCOLOR=\"gray\"><FONT POINT-SIZE=\"%d\" COLOR=\"#90ee90\">T</FONT></TD>\n",(int)tmp_sum[3]);
                fprintf(f_ptr,"</TR>\n");
                */
                
                fprintf(f_ptr,"</TABLE>>];\n");
                LOG_MSG("i:%d",i);
                
                
        }
        fprintf(f_ptr,"\n\n");
        /* print edges */
     
       
        for(i = 0;i < ft->last_state;i++){
                for(j = 0;j < ft->last_state;j++){
                        if(ft->transition[i][j] >= 1e-5){
                                RUN(get_color(color_buffer,ft->transition[i][j], 0.0f,1.0f ));
                                fprintf(f_ptr,"State%d -> State%d[label=\"%0.2f\",color=\"%s\", penwidth=%d];\n",i,j, ft->transition[i][j] , color_buffer, (int) (ft->transition[i][j] *10)+1 );
                        }
                       
                }
        }
        
        /* print end of dot file  */
        fprintf(f_ptr,"}\n");



        fclose(f_ptr);

        LOG_MSG("To visualize: dot  -Tpdf  <.dot file>  -o  <blah.pdf>.");
        MFREE(tmp_sum);
        MFREE(total_counts);
        return OK;
ERROR:
        if(tmp_sum){
                MFREE(tmp_sum);
        }
        if(total_counts){
                MFREE(total_counts);
        }
        return FAIL;
}


int get_color(char* color, float x, float min_x, float max_x)
{
        float range = 0;
        float r,g,b;
        
        ASSERT(color != NULL, "No colour buffer.");

        ASSERT(max_x >= min_x, "Max smaller than min.");

        /* Sanity */
        if(x < min_x){
                x = min_x;
        }
        if(x > max_x){
                x = max_x;
        }
        range = max_x - min_x;
        r = 1.0;
        g = 1.0;
        b = 1.0;
        
        if (x < (min_x + 0.25 * range)) {
                r = 0;
                g = 4.0f * (x - min_x) / range;
        } else if (x < (min_x + 0.5 * range)) {
                r = 0.0f;
                b = 1.0f + 4.0f * (min_x + 0.25f * range - x) / range;
        } else if (x < (min_x + 0.75 * range)) {
                r = 4.0f * (x - min_x - 0.5f * range) / range;
                b = 0;
        } else {
                g = 1.0f + 4.0f * (min_x + 0.75f * range - x) / range;
                b = 0;
        }

        /* Second go (from kalignvu 2002 code.. ) */
        //	( *  255 (* -1 (  log (* 2 0.6))))
        r = 0.0f;
        g = 0.0f;
        b = 0.0f;
        
        if(((x - min_x) / range) < 0.5){
                b = -1 * log(2.0*((x - min_x) / range)  );
                
        }else{
                r += log(2.0*((x - min_x) / range)  );
        }
        if(b > 1.0){
                b = 1.0;
        }
        if(r > 1.0){
                r = 1.0;
        }

        /*if(x < 0.5){
                b = 255;
        }else if(x < 1.0){
                b = (int)(-log(x)  * 255.0);
        }else if(x < 2.0){
                r = (int)(log(x) * 255.0);
        }else{
                r = 255;
                }*/
        snprintf(color, BUFFER_LEN,"#%02x%02x%02x",(unsigned int) (r * 255.0f), (unsigned int) (g * 255.0f),(unsigned int) (b * 255.0f));
                
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

        //    	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--seed-step","Distance between seeds." ,"[8]"  );
        return OK;
}


