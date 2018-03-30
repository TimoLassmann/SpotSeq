#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>


#include "tldevel.h"


#include "matrix_io.h"

#include "make_dot_file.h"

int make_dot_from_matrix(char* in_matrix, char* out_dot)
{
        struct double_matrix* matrix = NULL;
        FILE* f_ptr = NULL;

        double tmp_sum[4];

        double background[4];
        double IC;
        double max_IC;
        double tmp_prob;
        double sum_usage; 
        int ncol;
        int i,j;

        int max_stack_height = 128;

        ASSERT(in_matrix != NULL, "No input matrix.");
        ASSERT(out_dot != NULL, "No output dot file.");
        
        RUNP(matrix = read_double_matrix(in_matrix,1,1));
        print_double_matrix(matrix,stdout,1,1);
        ncol = matrix->ncol;

        for(i = 0; i < 4;i++){
                background[i] = scaledprob2prob( matrix->matrix[i][ncol-1]);
                fprintf(stdout,"%f\n",background[i]);
        }
        
        RUNP(f_ptr = fopen( out_dot, "w"));
        
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
                        //       fprintf(stdout,"%f\t%f\t%f\n",tmp_prob  ,  background[j],log2(tmp_prob / background[j]));
                }

                //fprintf(stdout,"%f\n",IC);
                if(IC> max_IC){
                        max_IC = IC;
                }
        }

        /* print start stop nodes  */
         
        fprintf(f_ptr,"%s [label=Start]\n", matrix->col_names[0]);
        fprintf(f_ptr,"%s [label=End]\n", matrix->col_names[1]);
         
        sum_usage = 0;
        
        for(i = 2;i < ncol-1;i++){
                if(matrix->matrix[4][i] >  sum_usage){
                        sum_usage = matrix->matrix[4][i];
                }
        }
        fprintf(stdout,"%f\n",sum_usage);

        
        
        /* print nodes...  */
        for(i = 2;i < ncol-1;i++){
                IC = 0.0;

                for(j = 0; j < 4;j++){
                        IC += matrix->matrix[j][i] * log2( matrix->matrix[j][i] / background[j]);
                        //                    fprintf(stdout,"%f\t%f\t%f\n", matrix->matrix[j][i],  background[j],log2( matrix->matrix[j][i] / background[j]));
                }
//                fprintf(stdout,"%f\n",IC);
                /* scale total bar height by information content */
                tmp_prob = (double) max_stack_height *((double) matrix->matrix[4][i] /sum_usage);// max_IC *IC      ;

                /* Further scale bar height by usage of state.  */
                //tmp_prob = tmp_prob * (double) matrix->matrix[4][i] / sum_usage;
                fprintf(stdout,"%d %f %f\n",i,tmp_prob ,(double) matrix->matrix[4][i] / sum_usage);


                        
                for(j = 0; j < 4;j++){
                        
                        tmp_sum[j] =  matrix->matrix[j][i] * tmp_prob  ;
                                               fprintf(stdout,"\t%f\n",tmp_sum[j]);
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
                        if(matrix->matrix[i+6][j] >= 1e-5){
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

