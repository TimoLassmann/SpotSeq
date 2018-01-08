#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"

#include "matrix_io.h"
#include "ihmm.h"



int write_ihmm_parameters_to_file(struct iHMM_model* model,char* filename)
{
        FILE* f_ptr = NULL;
        struct double_matrix* matrix = NULL;
        int i,j;
        
        ASSERT(model != NULL,"No model.");

        /* allocate a matrix to hold emission parameter on top followed by
         * transistion parameters */

        /* oh oh annoying + 1's again = necessary here as there are <=
         * infinityghost states */
        RUNP(matrix = alloc_double_matrix(model->infinityghost +1, model->infinityghost + 6 + 1, 32));

        /* Column names */
        for(i = 0; i < model->infinityghost; i++){
                snprintf(matrix->col_names[i],32,"State%d",i+1);
        }
        /* First rows for emission and other  */
        for(i = 0; i < 4;i++){
                snprintf(matrix->row_names[i],32,"Emission%c","ACGT"[i]);
        }
        snprintf(matrix->row_names[4],32,"Used");
        snprintf(matrix->row_names[5],32,"Empty");

        /* Fill out remainign row names with transitions */
        for(i = 0; i<= model->infinityghost ;i++){
                snprintf(matrix->row_names[i+6],32,"Trans%d", i+1);
        }

        /* fill in emission and usage...  */
        for(i = 0; i <= model->infinityghost;i++){
                matrix->matrix[4][i] = model->sumM[i];
                for(j = 0;j < model->L;j++){
                        matrix->matrix[j][i] = model->emission[i][j];
                }
        }
        
        /* Fill out remaining matrix with transistion probabilities...  */
        for(i = 0; i <= model->infinityghost;i++){

                for(j = 0;j <= model->infinityghost;j++){
                        matrix->matrix[j+6][i] = model->transition[i][j];
                       
                }
        }
        
        RUNP(f_ptr = fopen(filename,"w"));
        RUN(print_double_matrix(matrix, f_ptr,1,1));

        fclose(f_ptr);

        free_double_matrix(matrix);
        
        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        free_double_matrix(matrix);
        return FAIL; 
}


struct iHMM_model* read_ihmm_parameters_from_file(char* filename)
{
        struct iHMM_model* model = NULL;
        ASSERT(filename != NULL,"No filename.");
        return model;
ERROR:
        return NULL;
}
