#include "model.h"
#include "hdf5_glue.h"

struct ihmm_model* read_model_hdf5(char* filename)
{
        struct ihmm_model* model = NULL;
        struct hdf5_data* hdf5_data = NULL;
        int a,b;
        int i;
        ASSERT(filename != NULL, "No filename");
        ASSERT(my_file_exists(filename) != 0,"File %s does not exist.",filename);


        hdf5_data = hdf5_create();

        hdf5_open_file(filename,hdf5_data);
        hdf5_read_attributes(hdf5_data,hdf5_data->file);
        print_attributes(hdf5_data);
        get_group_names(hdf5_data);
        fprintf(stdout,"Groups:\n");
        for(i = 0; i < hdf5_data->grp_names->num_names;i++){
                fprintf(stdout,"%d %s\n",i,hdf5_data->grp_names->names[i]);
        }

        hdf5_open_group("imodel",hdf5_data);
        hdf5_read_attributes(hdf5_data, hdf5_data->group);
        ASSERT(hdf5_data->num_attr != 0 , "Could not find attributes");
        print_attributes(hdf5_data);
        a = 0;
        b = 0;
        for(i = 0; i < hdf5_data->num_attr;i++){
                if(!strncmp("Number of states", hdf5_data->attr[i]->attr_name, 16)){
                        a = hdf5_data->attr[i]->int_val;
                }
                if(!strncmp("Number of letters", hdf5_data->attr[i]->attr_name, 17)){
                        b = hdf5_data->attr[i]->int_val;
                }

        }
        ASSERT(a!=0, "No states???");
        ASSERT(b!=0, "No letters???");
        RUNP(model = alloc_ihmm_model(a, b));
        free_2d((void**) model->emission_counts);
        free_2d((void**) model->transition_counts);
        MFREE(model->beta);
        for(i = 0; i < hdf5_data->num_attr;i++){
                if(!strncmp("Number of states", hdf5_data->attr[i]->attr_name, 16)){
                        model->num_states = hdf5_data->attr[i]->int_val;}
                if(!strncmp("Number of letters", hdf5_data->attr[i]->attr_name, 17)){
                        model->L = hdf5_data->attr[i]->int_val;}

                if(!strncmp("Gamma", hdf5_data->attr[i]->attr_name, 5)){
                        model->gamma = hdf5_data->attr[i]->float_val;
                }
                if(!strncmp("gamma_a", hdf5_data->attr[i]->attr_name, 7)){
                        model->gamma_a = hdf5_data->attr[i]->float_val;
                }
                if(!strncmp("gamma_b", hdf5_data->attr[i]->attr_name, 7)){
                        model->gamma_b = hdf5_data->attr[i]->float_val;
                }
                if(!strncmp("Alpha", hdf5_data->attr[i]->attr_name, 5)){
                        model->alpha = hdf5_data->attr[i]->float_val;
                }
                if(!strncmp("alpha0_a", hdf5_data->attr[i]->attr_name, 8)){
                        model->alpha0_a = hdf5_data->attr[i]->float_val;
                }
                if(!strncmp("alpha0_b", hdf5_data->attr[i]->attr_name, 8)){
                        model->alpha0_b = hdf5_data->attr[i]->float_val;
                }
        }

        hdf5_read_dataset("Beta",hdf5_data);

        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 1, "Could not read beta");
        model->beta = (float*)hdf5_data->data;


        hdf5_read_dataset("transition_counts",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read transition_counts");
        model->transition_counts = (float**)hdf5_data->data;

        hdf5_read_dataset("emission_counts",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read emission_counts");
        model->emission_counts = (float**)hdf5_data->data;

        hdf5_close_group(hdf5_data);

        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);

        rk_randomseed(&model->rndstate);
        return model;
ERROR:
        return NULL;
}

int add_fhmm(char* filename,float** transition,float** emission, int N, int L)
{
        struct hdf5_data* hdf5_data = NULL;
        RUNP(hdf5_data = hdf5_create());

        hdf5_open_file(filename,hdf5_data);

        hdf5_create_group("fmodel",hdf5_data);

        hdf5_data->rank = 2;
        hdf5_data->dim[0] = N;
        hdf5_data->dim[1] = L;
        hdf5_data->chunk_dim[0] = N;
        hdf5_data->chunk_dim[1] = L;
        hdf5_write_2D_float("emission",emission, hdf5_data);

        hdf5_data->rank = 2;
        hdf5_data->dim[0] = N;
        hdf5_data->dim[1] = N;
        hdf5_data->chunk_dim[0] = N;
        hdf5_data->chunk_dim[1] = N;
        hdf5_write_2D_float("transition", transition, hdf5_data);

        hdf5_close_group(hdf5_data);
        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);
        return OK;

ERROR:
        if(hdf5_data){
                hdf5_close_file(hdf5_data);
                hdf5_free(hdf5_data);
        }
        return FAIL;
}

int add_background_emission(char* filename,float* background,int L)
{
        struct hdf5_data* hdf5_data = NULL;
        RUNP(hdf5_data = hdf5_create());

        hdf5_open_file(filename,hdf5_data);

        hdf5_create_group("SequenceInformation",hdf5_data);

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = L;
        hdf5_data->dim[1] = -1;
        hdf5_data->chunk_dim[0] = L;
        hdf5_data->chunk_dim[1] = -1;
        hdf5_write_1D_float("background",background, hdf5_data);


        hdf5_close_group(hdf5_data);
        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);
        return OK;
ERROR:
        if(hdf5_data){
                hdf5_close_file(hdf5_data);
                hdf5_free(hdf5_data);
        }
        return FAIL;
}

int add_annotation( char* filename, char* name, char* value)
{
        struct hdf5_data* hdf5_data = NULL;
        RUNP(hdf5_data = hdf5_create());
        hdf5_add_attribute(hdf5_data, name, value, 0, 0.0f, HDF5GLUE_CHAR);
        hdf5_open_file(filename,hdf5_data);
        hdf5_write_attributes(hdf5_data, hdf5_data->file);
        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);
        return OK;
ERROR:
        if(hdf5_data){
                hdf5_close_file(hdf5_data);
                hdf5_free(hdf5_data);
        }
        return FAIL;
}

int write_model_hdf5(struct ihmm_model* model, char* filename)
{
        char buffer[BUFFER_LEN];
        struct hdf5_data* hdf5_data = NULL;

        RUNP(hdf5_data = hdf5_create());
        snprintf(buffer, BUFFER_LEN, "%s",PACKAGE_NAME );
        hdf5_add_attribute(hdf5_data, "Program", buffer, 0, 0.0f, HDF5GLUE_CHAR);
        snprintf(buffer, BUFFER_LEN, "%s", PACKAGE_VERSION);
        hdf5_add_attribute(hdf5_data, "Version", buffer, 0, 0.0f, HDF5GLUE_CHAR);


        hdf5_create_file(filename,hdf5_data);

        hdf5_write_attributes(hdf5_data, hdf5_data->file);

        hdf5_data->num_attr = 0;

        hdf5_add_attribute(hdf5_data, "Number of states", "", model->num_states, 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "Number of letters", "", model->L, 0.0f, HDF5GLUE_INT);

        hdf5_add_attribute(hdf5_data, "Gamma", "", 0, model->gamma, HDF5GLUE_FLOAT);
        hdf5_add_attribute(hdf5_data, "gamma_a","",0, model->gamma_a, HDF5GLUE_FLOAT);
        hdf5_add_attribute(hdf5_data, "gamma_b","",0, model->gamma_b, HDF5GLUE_FLOAT);

        hdf5_add_attribute(hdf5_data, "Alpha",    "",0, model->alpha, HDF5GLUE_FLOAT);
        hdf5_add_attribute(hdf5_data, "alpha0_a", "",0, model->alpha0_a, HDF5GLUE_FLOAT);
        hdf5_add_attribute(hdf5_data, "alpha0_b", "",0, model->alpha0_b, HDF5GLUE_FLOAT);

        hdf5_create_group("imodel",hdf5_data);
        hdf5_write_attributes(hdf5_data, hdf5_data->group);

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = model->alloc_num_states;
        hdf5_data->dim[1] = -1;
        hdf5_data->chunk_dim[0] = model->alloc_num_states;
        hdf5_data->chunk_dim[1] = -1;
        hdf5_write_1D_float("Beta",model->beta, hdf5_data);

        hdf5_data->rank = 2;
        hdf5_data->dim[0] = model->alloc_num_states;
        hdf5_data->dim[1] = model->alloc_num_states;
        hdf5_data->chunk_dim[0] = model->alloc_num_states;
        hdf5_data->chunk_dim[1] = model->alloc_num_states;

        hdf5_write_2D_float("transition_counts",model->transition_counts, hdf5_data);

        hdf5_data->rank = 2;
        hdf5_data->dim[0] = model->L;
        hdf5_data->dim[1] = model->alloc_num_states;
        hdf5_data->chunk_dim[0] = model->L;
        hdf5_data->chunk_dim[1] = model->alloc_num_states;

        hdf5_write_2D_float("emission_counts",model->emission_counts, hdf5_data);



        /*
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
        */

        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);
        return OK;
ERROR:
        if(hdf5_data){
                hdf5_close_file(hdf5_data);
                hdf5_free(hdf5_data);
        }
        return FAIL;
}



/* Function to write and read  a model (hyperparameters and counts) */
int write_model(struct ihmm_model* model, char* filename)
{
        FILE* f_ptr = NULL;
        int i,j;


        ASSERT(model != NULL, "No model.");

        RUNP(f_ptr = fopen(filename, "w"));
        fprintf(f_ptr,"Number of states: %d\n", model->num_states);
        fprintf(f_ptr,"Number of letters: %d\n", model->L);
        fprintf(f_ptr,"Target length: %d\n", model->target_len);

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
        fscanf(f_ptr,"Target length: %d\n", &model->target_len);
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
