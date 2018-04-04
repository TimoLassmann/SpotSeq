#include "model.h"

/* Function to write and read  a model (hyperparameters and counts) */
int write_model(struct ihmm_model* model, char* filename)
{
        FILE* f_ptr = NULL;
        int i,j;
        
        
        ASSERT(model != NULL, "No model.");

        RUNP(f_ptr = fopen(filename, "w"));
        fprintf(f_ptr,"Number of states: %d\n", model->num_states);
        fprintf(f_ptr,"Number of letters: %d\n", model->L);
        
        fprintf(f_ptr,"Gamma: %f\n", model->gamma);
        fprintf(f_ptr,"gamma_a: %f\n", model->gamma_a);
        fprintf(f_ptr,"gamma_b: %f\n", model->gamma_b);
        
        fprintf(f_ptr,"alpha: %f\n", model->alpha);
        fprintf(f_ptr,"alpha0_a: %f\n", model->alpha0_a);
        fprintf(f_ptr,"alpha0_b: %f\n", model->alpha0_b);

        
        fprintf(f_ptr,"Beta:\n");
        for(i = 0; i < model->num_states;i++){
                fprintf(f_ptr,"%f\n", model->beta[i]);
        }
        fprintf(f_ptr,"Transition counts:\n");
        for(i = 0; i < model->num_states;i++){
                for(j = 0; j < model->num_states;j++){
                        fprintf(f_ptr,"%f\n", model->transition_counts[i][j]);
                }
        }
        fprintf(f_ptr,"Emission counts:\n");
        for(i = 0; i < model->L;i++){
                for(j = 0; j < model->num_states;j++){
                        fprintf(f_ptr,"%f\n", model->emission_counts[i][j]);
                }
        }
        
        fclose(f_ptr);
        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        return FAIL;
}

struct ihmm_model* read_model(char* filename)
{
        struct ihmm_model* model = NULL;
        FILE* f_ptr = NULL;
        int a,b;
        int i,j;
        ASSERT(filename != NULL, "No filename");
        ASSERT(my_file_exists(filename) != 0,"File %s does not exist.",filename);
        RUNP(f_ptr = fopen(filename, "r"));
        fscanf(f_ptr,"Number of states: %d\n",&a);
        fscanf(f_ptr,"Number of letters: %d\n", &b);
        fprintf(stdout,"Number of states after reading! : %d\n",a);
        

     
        RUNP(model = alloc_ihmm_model(a, b));
        model->num_states = a;
        fscanf(f_ptr,"Gamma: %f\n", &model->gamma);
        fscanf(f_ptr,"gamma_a: %f\n", &model->gamma_a);
        fscanf(f_ptr,"gamma_b: %f\n", &model->gamma_b);
        
        fscanf(f_ptr,"alpha: %f\n", &model->alpha);
        fscanf(f_ptr,"alpha0_a: %f\n", &model->alpha0_a);
        fscanf(f_ptr,"alpha0_b: %f\n", &model->alpha0_b);

        
        fscanf(f_ptr, "%*[^\n]\n");
        for(i = 0; i < model->num_states;i++){
                fscanf(f_ptr,"%f\n", &model->beta[i]);
        }
        fscanf(f_ptr, "%*[^\n]\n");
        for(i = 0; i < model->num_states;i++){
                for(j = 0; j < model->num_states;j++){
                        fscanf(f_ptr,"%f\n", &model->transition_counts[i][j]);
                }
        }
        fscanf(f_ptr, "%*[^\n]\n");
         for(i = 0; i < model->L;i++){
                for(j = 0; j < model->num_states;j++){
                        fscanf(f_ptr,"%f\n", &model->emission_counts[i][j]);
                }
        }
        rk_randomseed(&model->rndstate);
        fclose(f_ptr);
        return model;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        return NULL;
}
