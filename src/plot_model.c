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
        float edge_threshold;
        int node_count_cutoff;
};


static int run_plot_ihmm(struct parameters* param);
static int run_plot_positional_state_distribution(struct parameters* param);
static int make_dot_file(struct fast_hmm_param* ft, struct ihmm_model* model, struct parameters* param);

static int plot_model_entropy(struct parameters* param);

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
        param->edge_threshold = 0.5f;
        param->node_count_cutoff = 0;
        while (1){
                static struct option long_options[] ={
                        {"in",required_argument,0,'i'},
                        {"out",required_argument,0,'o'},
                        {"ethres",required_argument,0,'e'},
                        {"nthres",required_argument,0,'n'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"hi:o:e:n:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 'n':
                        param->node_count_cutoff = atoi(optarg);
                        break;
                case 'e':
                        param->edge_threshold = atof(optarg);
                        break;
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

        RUN(plot_model_entropy(param));

        RUN(run_plot_positional_state_distribution(param));

        //RUN(seed_controller_thread(param));

        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}

int hdf5_testing(struct parameters* param)
{

        ASSERT(param != NULL, "No parameters");
        RUN(write_model_hdf5(NULL,"TEST.h5"));
        return OK;
ERROR:
        return FAIL;
}

int run_plot_ihmm(struct parameters* param)
{
        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;
        struct seq_buffer* sb = NULL;

        int initial_states = 10;
        ASSERT(param!= NULL, "No parameters found.");

        RUNP(model = read_model_hdf5(param->input));

        RUNP(ft = alloc_fast_hmm_param(initial_states,model->L));
        RUN(fill_background_emission_from_model(ft,model));

        RUN(fill_fast_transitions_only_matrices(model,ft));
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

int run_plot_positional_state_distribution(struct parameters* param)
{
        struct seq_buffer* sb = NULL;
        struct ihmm_model* model = NULL;
        float** matrix = NULL;
        float* state_sums = NULL;

        int i,j,c,index;
        float l;

        FILE* fptr = NULL;
        ASSERT(param != NULL, "no parameters");
        RUNP(sb = get_sequences_from_hdf5_model(param->input));
        RUNP(model = read_model_hdf5(param->input));


        RUNP(matrix = malloc_2d_float(matrix, model->num_states , 201, 0.0f));
        MMALLOC(state_sums, sizeof(float) *  model->num_states);
        for(i = 0; i < model->num_states;i++){
                state_sums[i] = 0.0f;
        }

        for(i = 0; i < sb->num_seq;i++){
                l = (float) sb->sequences[i]->seq_len;
                for(j = 0; j < sb->sequences[i]->seq_len;j++){
                        c = (float) sb->sequences[i]->label[j];
                        index = roundf(200.0f * ((float) j / l));
                        matrix[c][index] += 1.0f;
                        state_sums[c] += 1.0f;
                }
        }
        /* First two states are START/ STOP */

        RUNP(fptr = fopen("test.csv", "w"));
        fprintf(fptr, "Pos");
        for(i = 2; i < model->num_states;i++){
                if(state_sums[i]){
                        fprintf(fptr,",State%d",i);
                }
        }
        fprintf(fptr,"\n");
        for(j = 0; j < 200;j++){
                fprintf(fptr, "%d",j);
                for(i = 2; i < model->num_states;i++){
                        if(state_sums[i]){
                                fprintf(fptr,",%0.6f",matrix[i][j]);
                        }
                }
                fprintf(fptr,"\n");
        }
        fclose(fptr);
        MFREE(state_sums);
        free_2d((void**) matrix);
        free_ihmm_model(model);
        free_ihmm_sequences(sb);
        return OK;
ERROR:
        if(state_sums){
                MFREE(state_sums);
        }
        if(sb){
                free_ihmm_sequences(sb);
        }
        free_2d((void**) matrix);
        free_ihmm_model(model);
        return FAIL;
}

int plot_model_entropy(struct parameters* param)
{
        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;

        int initial_states = 10;
        int iter;
        int i,j,c;
        float** s1 = NULL;
        float** s2 = NULL;

        int iterations = 1000;

        ASSERT(param!= NULL, "No parameters found.");


        RUNP(model = read_model_hdf5(param->input));
        RUNP(ft = alloc_fast_hmm_param(initial_states,model->L));

        /* first index is state * letter ; second is sample (max = 100) */

        s1 = malloc_2d_float(s1, model->num_states, model->L, 0.0);
        s2 = malloc_2d_float(s2, model->num_states, model->L, 0.0);


        for(iter=  0;iter < iterations;iter++){
                RUN(fill_fast_transitions_only_matrices(model,ft));
                for(i = 0;i < model->num_states;i++){
                        for(c = 0; c < model->L;c++){
                                s1[i][c] += ft->emission[c][i];
                                s2[i][c] += (ft->emission[c][i] * ft->emission[c][i]);
                        }
                }
        }

        for(i = 0; i < model->num_states;i++){
                for(j = 0; j < model->L;j++){
                        s2[i][j] = sqrt(  ((double) iterations * s2[i][j] - s1[i][j] * s1[i][j])/ ((double) iterations * ((double) iterations -1.0)));
                        s1[i][j] = s1[i][j] / (double) iterations;
                }
        }
        free_2d((void**) s1);
        free_2d((void**) s2);


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
                if(total_counts[i] >= param->node_count_cutoff){
                        IC = 0.0;

                        for(j = 0; j < ft->L;j++){
                                IC += ft->emission[j][i] * log2(ft->emission[j][i] / background[j]);

                                //                    fprintf(stdout,"%f\t%f\t%f\n", matrix->matrix[j][i],  background[j],log2( matrix->matrix[j][i] / background[j]));
                        }
//                fprintf(stdout,"%f\n",IC);
                        /* scale total bar height by information content */
                        tmp_prob =  (double) max_stack_height;// * (((double) total_counts[i]) / (double)sum_usage);// / model->emission_counts[0][i];// *((double) model->emission_counts[0][i]  /sum_usage);// max_IC *IC      ;

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
                        fprintf(f_ptr,"<TD BGCOLOR=\"gray\"><FONT POINT-SIZE=\"%d\">%d: %d</FONT></TD>\n",12,i,total_counts[i]);
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
        }
        fprintf(f_ptr,"\n\n");
        /* print edges */


        for(i = 0;i < ft->last_state;i++){
                if(total_counts[i] >= param->node_count_cutoff){
                for(j = 0;j < ft->last_state;j++){
                          if(total_counts[j] >= param->node_count_cutoff){
                        if(ft->transition[i][j] >= param->edge_threshold ){
                                RUN(get_color(color_buffer,ft->transition[i][j], 0.0f,1.0f ));
                                fprintf(f_ptr,"State%d -> State%d[label=\"%0.2f\",color=\"%s\", penwidth=%d];\n",i,j,  ft->transition[i][j] , color_buffer, (int) (ft->transition[i][j] *10)+1 );
                        }
                          }

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
        const char usage[] = " -i <input model> -o <output dot file>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthres","Distance between seeds." ,"[0]"  );

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--ethres","Edge threshold." ,"[0.5]"  );
        return OK;
}


