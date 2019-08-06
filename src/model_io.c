#include "model.h"
#include "hdf5_glue.h"
#include "thread_data.h"

static int add_RNG_state(struct hdf5_data* hdf5_data, char* group,rk_state* a);

static int read_RNG_state(struct hdf5_data* hdf5_data, char* group,rk_state* a);


static struct ihmm_model* read_model_hdf5(struct hdf5_data* hdf5_data,char* group);
static int write_model_hdf5(struct hdf5_data* hdf5_data,struct ihmm_model* model, char* group);

struct ihmm_model* read_best_imodel(char* filename, int* best_model)
{
        struct ihmm_model* m = NULL;
        struct model_bag* model_bag = NULL;

        RUNP( model_bag = read_model_bag_hdf5(filename));
        m = model_bag->models[model_bag->best_model];
        model_bag->models[model_bag->best_model] = NULL;
        *best_model = model_bag->best_model;
        free_model_bag(model_bag);
        return m;
ERROR:
        return NULL;
}

struct fhmm* read_best_fmodel(char* filename, int* best_model)
{

        struct fhmm* m = NULL;
        struct model_bag* model_bag = NULL;

        RUNP( model_bag = read_model_bag_hdf5(filename));
        m = model_bag->finite_models[model_bag->best_model];
        model_bag->finite_models[model_bag->best_model] = NULL;
        *best_model = model_bag->best_model;
        free_model_bag(model_bag);
        return m;
ERROR:
        return NULL;
}

struct model_bag* read_model_bag_hdf5(char* filename)
{
        char buffer[BUFFER_LEN];
        struct model_bag* bag = NULL;
        struct hdf5_data* hdf5_data = NULL;

        int i;
        ASSERT(filename != NULL, "No filename");
        ASSERT(my_file_exists(filename) != 0,"File %s does not exist.",filename);



        MMALLOC(bag, sizeof(struct model_bag));

        bag->models = NULL;
        bag->finite_models = NULL;
        bag->min_u = NULL;
        bag->num_models = 0;
        bag->max_num_states = 0;
        bag->best_model = -1;

        hdf5_data = hdf5_create();

        hdf5_open_file(filename,hdf5_data);
        hdf5_read_attributes(hdf5_data,hdf5_data->file);
        //print_attributes(hdf5_data);

        for(i = 0; i < hdf5_data->num_attr;i++){
                if(!strcmp("MaxStates", hdf5_data->attr[i]->attr_name)){
                        bag->max_num_states = hdf5_data->attr[i]->int_val;
                }
                if(!strcmp("Number of models", hdf5_data->attr[i]->attr_name)){
                        bag->num_models = hdf5_data->attr[i]->int_val;
                }
                if(!strcmp("Seed", hdf5_data->attr[i]->attr_name)){
                        bag->seed = hdf5_data->attr[i]->int_val;
                }
        }
        ASSERT(bag->num_models > 0, "No models!");
        MMALLOC(bag->models, sizeof(struct ihmm_model*)* bag->num_models);


        MMALLOC(bag->finite_models, sizeof(struct fhmm*)* bag->num_models);
        MMALLOC(bag->min_u , sizeof(double) * bag->num_models);

        //get_group_names(hdf5_data);
        //fprintf(stdout,"Groups:\n");
        //for(i = 0; i < hdf5_data->grp_names->num_names;i++){
        //        fprintf(stdout,"%d %s\n",i,hdf5_data->grp_names->names[i]);
        //}
        RUN(hdf5_open_group("/",hdf5_data ));
        if((hdf5_data->dataset = H5Dopen(hdf5_data->group, "BestModel",H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        //printf ("H5Dopen returns: %d\n", hdf5_data->dataset);
        hdf5_read_attributes(hdf5_data,hdf5_data->dataset);
        hdf5_data->datatype  = H5Dget_type(hdf5_data->dataset );     /* datatype handle */
        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        hdf5_data->rank      = H5Sget_simple_extent_ndims(hdf5_data->dataspace);
        hdf5_data->status  = H5Sget_simple_extent_dims(hdf5_data->dataspace,hdf5_data->dim , NULL);
        hdf5_data->status = H5Dread(hdf5_data->dataset, hdf5_data->datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bag->best_model);
        if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");
        RUN(hdf5_close_group(hdf5_data));



        RUN(read_RNG_state(hdf5_data, "RNG",&bag->rndstate));

        for(i = 0; i < bag->num_models;i++){
                bag->models[i] = NULL;
                snprintf(buffer, BUFFER_LEN, "models/m%d",i+1);
                RUNP(bag->models[i] = read_model_hdf5(hdf5_data, buffer));
                if(bag->max_num_states < bag->models[i]->num_states){
                        bag->max_num_states = bag->models[i]->num_states;
                }


                bag->finite_models[i] = NULL;

                snprintf(buffer, BUFFER_LEN, "models/m%d/fhmm",i+1);
                /* get HMM parameters  */
                RUNP(bag->finite_models[i] = read_fhmm_parameters(hdf5_data, buffer));

        }

        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);
        return bag;
ERROR:
        return NULL;
}

int write_best_model(char* filename, int best_model)
{
        struct hdf5_data* hdf5_data = NULL;

        void *ptr;
        ASSERT(filename != NULL, "No filename");
        ASSERT(my_file_exists(filename) != 0,"File %s does not exist.",filename);


        hdf5_data = hdf5_create();

        hdf5_open_file(filename,hdf5_data);

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = 1;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->native_type = H5T_NATIVE_INT;
        ptr = (void*) &best_model;
        RUN(hdf5_open_group("/",hdf5_data));
        if((hdf5_data->dataset = H5Dopen(hdf5_data->group,"BestModel",H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        if((hdf5_data->status  = H5Dwrite(hdf5_data->dataset,hdf5_data->native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr)) < 0) ERROR_MSG("H5Dwrite failed");
        RUN(hdf5_close_group(hdf5_data));

        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);
        return OK;
ERROR:
        return FAIL;

}

int write_model_bag_hdf5(struct model_bag* bag, char* filename)
{
        char buffer[BUFFER_LEN];
        struct hdf5_data* hdf5_data = NULL;
        int model_i;
        int best_model = -1;
        RUNP(hdf5_data = hdf5_create());
        snprintf(buffer, BUFFER_LEN, "%s",PACKAGE_NAME );
        hdf5_add_attribute(hdf5_data, "Program", buffer, 0, 0.0f, HDF5GLUE_CHAR);
        snprintf(buffer, BUFFER_LEN, "%s", PACKAGE_VERSION);
        hdf5_add_attribute(hdf5_data, "Version", buffer, 0, 0.0f, HDF5GLUE_CHAR);
        hdf5_add_attribute(hdf5_data, "Number of models", "", bag->num_models, 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "MaxStates", "", bag->max_num_states, 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "Seed", "", bag->seed , 0.0f, HDF5GLUE_INT);


        RUN(hdf5_create_file(filename,hdf5_data));

        hdf5_write_attributes(hdf5_data, hdf5_data->file);

        hdf5_data->num_attr = 0;

        //LOG_MSG("header done");


        RUN(hdf5_create_group("/",hdf5_data));
        hdf5_data->rank = 1;
        hdf5_data->dim[0] = 1;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->native_type = H5T_NATIVE_INT;
        RUN(hdf5_write("BestModel",&best_model, hdf5_data));
        RUN(hdf5_close_group(hdf5_data));

        RUN(add_RNG_state(hdf5_data, "RNG",&bag->rndstate));
        //LOG_MSG("RNG top  done");
        RUN(hdf5_create_group("models",hdf5_data));

        RUN(hdf5_close_group(hdf5_data));
        //LOG_MSG("Create model group done.");

        for(model_i = 0; model_i < bag->num_models;model_i++){
                snprintf(buffer, BUFFER_LEN, "models/m%d", model_i+1);
                //LOG_MSG("Writing model %d.",model_i);
                RUN(write_model_hdf5(hdf5_data, bag->models[model_i], buffer));

                snprintf(buffer, BUFFER_LEN, "models/m%d/fhmm", model_i+1);
                RUN(add_fhmm(hdf5_data, bag->finite_models[model_i] , buffer));
                //LOG_MSG("Done");
        }
        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);


        /*struct model_bag* bag2 = NULL;
        bag2 = read_model_bag_hdf5(filename);
        int i;
        for(i = 0;i < 10;i++){
                fprintf(stdout,"MAIN:%f %f\n",rk_double(&bag->rndstate),rk_double(&bag2->rndstate));
        }
        for(i = 0; i < bag2->num_models;i++){
                fprintf(stdout,"%d:%f %f\n",i,rk_double(&bag->models[i]->rndstate ),rk_double(&bag2->models[i]->rndstate));
        }

        free_model_bag(bag2);*/
        return OK;
ERROR:
        if(hdf5_data){
                hdf5_close_file(hdf5_data);
                hdf5_free(hdf5_data);
        }
        return FAIL;
}


struct ihmm_model* read_model_hdf5(struct hdf5_data* hdf5_data,char* group)
{
        char buffer[BUFFER_LEN+5];
        struct ihmm_model* model = NULL;
        int a,b;
        int i;
        //ASSERT(filename != NULL, "No filename");
        //ASSERT(my_file_exists(filename) != 0,"File %s does not exist.",filename);


        //hdf5_data = hdf5_create();

        //hdf5_open_file(filename,hdf5_data);
        //hdf5_read_attributes(hdf5_data,hdf5_data->file);
        //print_attributes(hdf5_data);
        //get_group_names(hdf5_data);
        //fprintf(stdout,"Groups:\n");
        //for(i = 0; i < hdf5_data->grp_names->num_names;i++){
        //        fprintf(stdout,"%d %s\n",i,hdf5_data->grp_names->names[i]);
        //}

        hdf5_open_group(group,hdf5_data);
        hdf5_read_attributes(hdf5_data, hdf5_data->group);
        ASSERT(hdf5_data->num_attr != 0 , "Could not find attributes");
        //print_attributes(hdf5_data);
        a = 0;
        b = 0;
        for(i = 0; i < hdf5_data->num_attr;i++){
                if(!strcmp("Number of states", hdf5_data->attr[i]->attr_name)){
                        a = hdf5_data->attr[i]->int_val;
                }
                if(!strcmp("Number of letters", hdf5_data->attr[i]->attr_name)){
                        b = hdf5_data->attr[i]->int_val;
                }
        }
        ASSERT(a!=0, "No states???");
        ASSERT(b!=0, "No letters???");
        RUNP(model = alloc_ihmm_model(a, b,0));
        gfree(model->emission_counts);
        gfree(model->transition_counts);
        MFREE(model->beta);
        for(i = 0; i < hdf5_data->num_attr;i++){
                if(!strcmp("Number of states", hdf5_data->attr[i]->attr_name)){
                        model->num_states = hdf5_data->attr[i]->int_val;
                }
                if(!strcmp("Number of letters", hdf5_data->attr[i]->attr_name)){
                        model->L = hdf5_data->attr[i]->int_val;
                }

                if(!strcmp("Gamma", hdf5_data->attr[i]->attr_name)){
                        model->gamma = hdf5_data->attr[i]->double_val;
                }
                if(!strcmp("gamma_a", hdf5_data->attr[i]->attr_name)){
                        model->gamma_a = hdf5_data->attr[i]->double_val;
                }
                if(!strcmp("gamma_b", hdf5_data->attr[i]->attr_name)){
                        model->gamma_b = hdf5_data->attr[i]->double_val;
                }
                if(!strcmp("Alpha", hdf5_data->attr[i]->attr_name)){
                        model->alpha = hdf5_data->attr[i]->double_val;
                }
                if(!strcmp("alpha0_a", hdf5_data->attr[i]->attr_name)){
                        model->alpha_a = hdf5_data->attr[i]->double_val;
                }
                if(!strcmp("alpha0_b", hdf5_data->attr[i]->attr_name)){
                        model->alpha_b = hdf5_data->attr[i]->double_val;
                }
                if(!strcmp("Iteration", hdf5_data->attr[i]->attr_name)){
                        model->training_iterations = hdf5_data->attr[i]->int_val;
                }
                if(!strcmp("Seed", hdf5_data->attr[i]->attr_name)){
                        model->seed = hdf5_data->attr[i]->int_val;
                }
                if(!strcmp("Alpha_limit", hdf5_data->attr[i]->attr_name)){
                        model->alpha_limit =  hdf5_data->attr[i]->double_val;
                }
                if(!strcmp("Gamma_limit", hdf5_data->attr[i]->attr_name)){
                        model->gamma_limit =  hdf5_data->attr[i]->double_val;
                }
        }

        hdf5_read_dataset("Beta",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 1, "Could not read beta");
        model->beta = (double*)hdf5_data->data;

        hdf5_read_dataset("transition_counts",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read transition_counts");
        model->transition_counts = (double**)hdf5_data->data;

        hdf5_read_dataset("emission_counts",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read emission_counts");
        model->emission_counts = (double**)hdf5_data->data;

        hdf5_close_group(hdf5_data);

        //hdf5_close_file(hdf5_data);
        //hdf5_free(hdf5_data);
        /* stretch matrices... */
        model->alloc_num_states = model->num_states;
        RUN(resize_ihmm_model(model, model->num_states+1));

        snprintf(buffer, BUFFER_LEN+5 , "%s/RNG", group);
        //LOG_MSG("Trying to create group: %s", buffer);
        RUN(read_RNG_state(hdf5_data, buffer,&model->rndstate));
        //WARNING_MSG("Each time a run is continued a new RNG seed is selected...");

        return model;
ERROR:
        return NULL;
}



int write_model_hdf5(struct hdf5_data* hdf5_data,struct ihmm_model* model, char* group)
{
        //struct hdf5_data* hdf5_data = NULL;
        char buffer[BUFFER_LEN+5];
        double** tmp = NULL;
        int i,j;

        //RUNP(hdf5_data = hdf5_create());
        //snprintf(buffer, BUFFER_LEN, "%s",PACKAGE_NAME );
        //hdf5_add_attribute(hdf5_data, "Program", buffer, 0, 0.0f, HDF5GLUE_CHAR);
        //snprintf(buffer, BUFFER_LEN, "%s", PACKAGE_VERSION);
        //hdf5_add_attribute(hdf5_data, "Version", buffer, 0, 0.0f, HDF5GLUE_CHAR);


        //RUN(hdf5_create_file(filename,hdf5_data));

        //hdf5_write_attributes(hdf5_data, hdf5_data->file);

        hdf5_data->num_attr = 0;

        hdf5_add_attribute(hdf5_data, "Number of states", "", model->num_states, 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "Number of letters", "", model->L, 0.0f, HDF5GLUE_INT);

        hdf5_add_attribute(hdf5_data, "Gamma", "", 0, model->gamma, HDF5GLUE_DOUBLE);
        hdf5_add_attribute(hdf5_data, "Gamma_limit", "", 0, model->gamma_limit, HDF5GLUE_DOUBLE);
        hdf5_add_attribute(hdf5_data, "gamma_a","",0, model->gamma_a, HDF5GLUE_DOUBLE);
        hdf5_add_attribute(hdf5_data, "gamma_b","",0, model->gamma_b, HDF5GLUE_DOUBLE);

        hdf5_add_attribute(hdf5_data, "Alpha",    "",0, model->alpha, HDF5GLUE_DOUBLE);
        hdf5_add_attribute(hdf5_data, "Alpha_limit",    "",0, model->alpha_limit , HDF5GLUE_DOUBLE);
        hdf5_add_attribute(hdf5_data, "alpha0_a", "",0, model->alpha_a, HDF5GLUE_DOUBLE);
        hdf5_add_attribute(hdf5_data, "alpha0_b", "",0, model->alpha_b, HDF5GLUE_DOUBLE);

        hdf5_add_attribute(hdf5_data, "Iteration", "",model->training_iterations, 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "Seed", "",model->seed, 0.0f, HDF5GLUE_INT);

        RUN(hdf5_create_group(group,hdf5_data));
        hdf5_write_attributes(hdf5_data, hdf5_data->group);
        hdf5_data->num_attr = 0;

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = model->num_states;
        hdf5_data->dim[1] = -1;
        hdf5_data->chunk_dim[0] = model->num_states;
        hdf5_data->chunk_dim[1] = -1;
        hdf5_data->native_type = H5T_NATIVE_DOUBLE;
        RUN(hdf5_write("Beta",&model->beta[0], hdf5_data));

        RUNP(tmp = galloc(tmp,  model->num_states,  model->num_states, 0.0));
        for(i = 0; i < model->num_states;i++){
                for(j = 0; j < model->num_states;j++){
                        tmp[i][j] = model->transition_counts[i][j];
                }
        }
        hdf5_data->rank = 2;
        hdf5_data->dim[0] = model->num_states;
        hdf5_data->dim[1] = model->num_states;
        hdf5_data->chunk_dim[0] = model->num_states;
        hdf5_data->chunk_dim[1] = model->num_states;
        hdf5_data->native_type = H5T_NATIVE_DOUBLE;
        RUN(hdf5_write("transition_counts",&tmp[0][0], hdf5_data));

        gfree(tmp);
        tmp= NULL;

        RUNP(tmp = galloc(tmp,  model->L ,  model->num_states, 0.0));
        for(i = 0; i < model->L;i++){
                for(j = 0; j < model->num_states;j++){
                        tmp[i][j] = model->emission_counts[i][j];
                }
        }
        hdf5_data->rank = 2;
        hdf5_data->dim[0] = model->L;
        hdf5_data->dim[1] = model->num_states;
        hdf5_data->chunk_dim[0] = model->L;
        hdf5_data->chunk_dim[1] = model->num_states;
        hdf5_data->native_type = H5T_NATIVE_DOUBLE;
        RUN(hdf5_write("emission_counts",&tmp[0][0], hdf5_data));

        gfree(tmp);
        tmp= NULL;

        RUN(hdf5_close_group(hdf5_data));

        snprintf(buffer, BUFFER_LEN+5 , "%s/RNG", group);
        //LOG_MSG("Trying to create group: %s", buffer);
        RUN(add_RNG_state(hdf5_data, buffer,&model->rndstate));

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

        //hdf5_close_file(hdf5_data);
        //hdf5_free(hdf5_data);
        return OK;
ERROR:
        //i/f(hdf5_data){
        //        hdf5_close_file(hdf5_data);
        //        hdf5_free(hdf5_data);
        //}
        return FAIL;
}




struct wims_thread_data** read_thread_data_to_hdf5(char* filename)
{
        struct wims_thread_data** td = NULL;
        struct hdf5_data* hdf5_data = NULL;
        char buffer[BUFFER_LEN];
        int num_threads = 0;
        int max_len = 0;
        int max_K = 0;
        int i;
        RUNP(hdf5_data = hdf5_create());

        hdf5_open_file(filename,hdf5_data);

        hdf5_open_group("thread_data",hdf5_data);
        hdf5_read_attributes(hdf5_data, hdf5_data->group);
        ASSERT(hdf5_data->num_attr != 0 , "Could not find attributes");
        //print_attributes(hdf5_data);
        for(i = 0; i < hdf5_data->num_attr;i++){
                if(!strcmp("Nthreads", hdf5_data->attr[i]->attr_name)){
                        num_threads = hdf5_data->attr[i]->int_val;
                }
                if(!strcmp("Max_L", hdf5_data->attr[i]->attr_name)){
                        max_len = hdf5_data->attr[i]->int_val;
                }
                if(!strcmp("Max_K", hdf5_data->attr[i]->attr_name)){
                        max_K = hdf5_data->attr[i]->int_val;
                }
        }
        ASSERT(num_threads!=0, "No threads???");
        ASSERT(max_len!=0, "No len???");
        ASSERT(max_K!=0, "No states???");

        hdf5_close_group(hdf5_data);
        RUNP(td = create_wims_thread_data(&num_threads, max_len, max_K, NULL));


        for(i = 0 ;i < num_threads;i++){
                snprintf(buffer, BUFFER_LEN , "thread_data/RNG%d",i);
                //LOG_MSG("Trying to create group: %s", buffer);
                RUN(read_RNG_state(hdf5_data, buffer, &td[i]->rndstate));
        }
        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);
        return td;
ERROR:
        return NULL;

}


int write_thread_data_to_hdf5(char* filename,struct wims_thread_data** td,int num_threads,int max_len,int max_K)
{
        struct hdf5_data* hdf5_data = NULL;
        char buffer[BUFFER_LEN];
        //unsigned int* seeds = NULL;
        int i;
        RUNP(hdf5_data = hdf5_create());

        hdf5_open_file(filename,hdf5_data);


        hdf5_data->num_attr = 0;

        hdf5_add_attribute(hdf5_data, "Nthreads", "",num_threads, 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "Max_L", "",max_len, 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "Max_K", "",max_K, 0.0f, HDF5GLUE_INT);
        hdf5_create_group("thread_data",hdf5_data);
        hdf5_write_attributes(hdf5_data, hdf5_data->group);
        hdf5_data->num_attr = 0;

        hdf5_close_group(hdf5_data);


        /* MMALLOC(seeds , sizeof(unsigned int)* num_threads); */
        /* for(i = 0; i < num_threads;i++){ */
        /*         seeds[i] = td[i]->seed; */
        /* } */




        for(i = 0; i < num_threads;i++){
                snprintf(buffer, BUFFER_LEN , "thread_data/RNG%d",i);
                //LOG_MSG("Trying to create group: %s", buffer);
                RUN(add_RNG_state(hdf5_data, buffer, &td[i]->rndstate));
        }

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

int add_fhmm(struct hdf5_data* hdf5_data, struct fhmm* fhmm, char* group)
{
        hdf5_create_group(group,hdf5_data);

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = 1;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->native_type = H5T_NATIVE_INT;
        RUN(hdf5_write("K",&fhmm->K, hdf5_data));

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = 1;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->native_type = H5T_NATIVE_INT;
        RUN(hdf5_write("L",&fhmm->L, hdf5_data));

        hdf5_data->rank = 2;
        hdf5_data->dim[0] = fhmm->K;
        hdf5_data->dim[1] = fhmm->L;
        hdf5_data->chunk_dim[0] = fhmm->K;
        hdf5_data->chunk_dim[1] = fhmm->L;
        hdf5_data->native_type = H5T_NATIVE_DOUBLE;
        hdf5_write("emission",  &fhmm->e[0][0], hdf5_data);

        hdf5_data->rank = 2;
        hdf5_data->dim[0] = fhmm->K;
        hdf5_data->dim[1] = fhmm->K;
        hdf5_data->chunk_dim[0] = fhmm->K;
        hdf5_data->chunk_dim[1] = fhmm->K;
        hdf5_data->native_type = H5T_NATIVE_DOUBLE;
        hdf5_write("transition", &fhmm->t[0][0] , hdf5_data);

        hdf5_data->rank = 2;
        hdf5_data->dim[0] = fhmm->K;
        hdf5_data->dim[1] = fhmm->K+1;
        hdf5_data->chunk_dim[0] = fhmm->K;
        hdf5_data->chunk_dim[1] = fhmm->K+1;
        hdf5_data->native_type = H5T_NATIVE_INT;
        hdf5_write("transition_index", &fhmm->tindex[0][0] , hdf5_data);

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = fhmm->L;
        hdf5_data->dim[1] = -1;
        hdf5_data->chunk_dim[0] = fhmm->L;
        hdf5_data->chunk_dim[1] = -1;
        hdf5_data->native_type = H5T_NATIVE_DOUBLE;
        hdf5_write("background", &fhmm->background[0], hdf5_data);

        hdf5_close_group(hdf5_data);

        return OK;
ERROR:
        if(hdf5_data){
                hdf5_close_file(hdf5_data);
                hdf5_free(hdf5_data);
        }
        return FAIL;
}

int add_background_emission(char* filename,double* background,int L)
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
        hdf5_data->native_type = H5T_NATIVE_DOUBLE;
        hdf5_write("background",&background[0], hdf5_data);


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




int read_RNG_state(struct hdf5_data* hdf5_data, char* group,rk_state* a)
{
        int i;

        ulong* tmp = NULL;

        //hdf5_read_attributes(hdf5_data,hdf5_data->file);
        //print_attributes(hdf5_data);
        //get_group_names(hdf5_data);
        //fprintf(stdout,"Groups:\n");
        //for(i = 0; i < hdf5_data->grp_names->num_names;i++){
        //        fprintf(stdout,"%d %s\n",i,hdf5_data->grp_names->names[i]);
        //}

        RUN(hdf5_open_group(group,hdf5_data));


        RUN(hdf5_read_dataset("key",hdf5_data));
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 1, "Could not read key");
        tmp= (ulong*)hdf5_data->data;
        for(i = 0; i < RK_STATE_LEN;i++){
                a->key[i] = tmp[i];
        }
        MFREE(tmp);


        if((hdf5_data->dataset = H5Dopen(hdf5_data->group, "gauss",H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        //printf ("H5Dopen returns: %d\n", hdf5_data->dataset);
        hdf5_read_attributes(hdf5_data,hdf5_data->dataset);
        hdf5_data->datatype  = H5Dget_type(hdf5_data->dataset );     /* datatype handle */
        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        hdf5_data->rank      = H5Sget_simple_extent_ndims(hdf5_data->dataspace);
        hdf5_data->status  = H5Sget_simple_extent_dims(hdf5_data->dataspace,hdf5_data->dim , NULL);
        hdf5_data->status = H5Dread(hdf5_data->dataset, hdf5_data->datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a->gauss);
        if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");


        if((hdf5_data->dataset = H5Dopen(hdf5_data->group, "psave",H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        //printf ("H5Dopen returns: %d\n", hdf5_data->dataset);
        hdf5_read_attributes(hdf5_data,hdf5_data->dataset);
        hdf5_data->datatype  = H5Dget_type(hdf5_data->dataset );     /* datatype handle */
        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        hdf5_data->rank      = H5Sget_simple_extent_ndims(hdf5_data->dataspace);
        hdf5_data->status  = H5Sget_simple_extent_dims(hdf5_data->dataspace,hdf5_data->dim , NULL);
        hdf5_data->status = H5Dread(hdf5_data->dataset, hdf5_data->datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a->psave);
        if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");


        if((hdf5_data->dataset = H5Dopen(hdf5_data->group, "has_binomial",H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        //printf ("H5Dopen returns: %d\n", hdf5_data->dataset);
        hdf5_read_attributes(hdf5_data,hdf5_data->dataset);
        hdf5_data->datatype  = H5Dget_type(hdf5_data->dataset );     /* datatype handle */
        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        hdf5_data->rank      = H5Sget_simple_extent_ndims(hdf5_data->dataspace);
        hdf5_data->status  = H5Sget_simple_extent_dims(hdf5_data->dataspace,hdf5_data->dim , NULL);
        hdf5_data->status = H5Dread(hdf5_data->dataset, hdf5_data->datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a->has_binomial);
        if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");


        if((hdf5_data->dataset = H5Dopen(hdf5_data->group, "has_gauss",H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        //printf ("H5Dopen returns: %d\n", hdf5_data->dataset);
        hdf5_read_attributes(hdf5_data,hdf5_data->dataset);
        hdf5_data->datatype  = H5Dget_type(hdf5_data->dataset );     /* datatype handle */
        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        hdf5_data->rank      = H5Sget_simple_extent_ndims(hdf5_data->dataspace);
        hdf5_data->status  = H5Sget_simple_extent_dims(hdf5_data->dataspace,hdf5_data->dim , NULL);
        hdf5_data->status = H5Dread(hdf5_data->dataset, hdf5_data->datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a->has_gauss);
        if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");




        if((hdf5_data->dataset = H5Dopen(hdf5_data->group, "nsave",H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        //printf ("H5Dopen returns: %d\n", hdf5_data->dataset);
        hdf5_read_attributes(hdf5_data,hdf5_data->dataset);
        hdf5_data->datatype  = H5Dget_type(hdf5_data->dataset );     /* datatype handle */
        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        hdf5_data->rank      = H5Sget_simple_extent_ndims(hdf5_data->dataspace);
        hdf5_data->status  = H5Sget_simple_extent_dims(hdf5_data->dataspace,hdf5_data->dim , NULL);
        hdf5_data->status = H5Dread(hdf5_data->dataset, hdf5_data->datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a->nsave);
        if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");

        if((hdf5_data->dataset = H5Dopen(hdf5_data->group, "pos",H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        //printf ("H5Dopen returns: %d\n", hdf5_data->dataset);
        hdf5_read_attributes(hdf5_data,hdf5_data->dataset);
        hdf5_data->datatype  = H5Dget_type(hdf5_data->dataset );     /* datatype handle */
        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        hdf5_data->rank      = H5Sget_simple_extent_ndims(hdf5_data->dataspace);
        hdf5_data->status  = H5Sget_simple_extent_dims(hdf5_data->dataspace,hdf5_data->dim , NULL);
        hdf5_data->status = H5Dread(hdf5_data->dataset, hdf5_data->datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a->pos);
        if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");

        hdf5_close_group(hdf5_data);
        return OK;
ERROR:
        return FAIL;
}

int add_RNG_state(struct hdf5_data* hdf5_data, char* group,rk_state* a)
{
        //struct hdf5_data* hdf5_data = NULL;
        RUN(hdf5_create_group(group,hdf5_data));

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = RK_STATE_LEN;
        hdf5_data->chunk_dim[0] = RK_STATE_LEN;
        hdf5_data->native_type = H5T_NATIVE_ULONG;
        RUN(hdf5_write("key",&a->key, hdf5_data));

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = 1;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->native_type = H5T_NATIVE_DOUBLE;
        RUN(hdf5_write("gauss",&a->gauss, hdf5_data));

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = 1;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->native_type = H5T_NATIVE_DOUBLE;
        RUN(hdf5_write("psave",&a->psave, hdf5_data));

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = 1;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->native_type = H5T_NATIVE_INT;
        RUN(hdf5_write("has_binomial",&a->has_binomial, hdf5_data));


        hdf5_data->rank = 1;
        hdf5_data->dim[0] = 1;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->native_type = H5T_NATIVE_INT;
        RUN(hdf5_write("has_gauss",&a->has_gauss, hdf5_data));


        hdf5_data->rank = 1;
        hdf5_data->dim[0] = 1;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->native_type = H5T_NATIVE_INT;
        RUN(hdf5_write("pos",&a->pos, hdf5_data));

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = 1;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->native_type = H5T_NATIVE_LONG;
        RUN(hdf5_write("nsave",&a->nsave, hdf5_data));

        RUN(hdf5_close_group(hdf5_data));
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

        fprintf(f_ptr,"Gamma: %f\n", model->gamma);
        fprintf(f_ptr,"gamma_a: %f\n", model->gamma_a);
        fprintf(f_ptr,"gamma_b: %f\n", model->gamma_b);

        fprintf(f_ptr,"alpha: %f\n", model->alpha);
        fprintf(f_ptr,"alpha_a: %f\n", model->alpha_a);
        fprintf(f_ptr,"alpha_b: %f\n", model->alpha_b);

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



        RUNP(model = alloc_ihmm_model(a, b,0));
        model->num_states = a;

        fscanf(f_ptr,"Gamma: %lf\n", &model->gamma);
        fscanf(f_ptr,"gamma_a: %lf\n", &model->gamma_a);
        fscanf(f_ptr,"gamma_b: %lf\n", &model->gamma_b);

        fscanf(f_ptr,"alpha: %lf\n", &model->alpha);
        fscanf(f_ptr,"alpha0_a: %lf\n", &model->alpha_a);
        fscanf(f_ptr,"alpha0_b: %lf\n", &model->alpha_b);


        fscanf(f_ptr, "%*[^\n]\n");
        for(i = 0; i < model->num_states;i++){
                fscanf(f_ptr,"%lf\n", &model->beta[i]);
        }
        fscanf(f_ptr, "%*[^\n]\n");
        for(i = 0; i < model->num_states;i++){
                for(j = 0; j < model->num_states;j++){
                        fscanf(f_ptr,"%lf\n", &model->transition_counts[i][j]);
                }
        }
        fscanf(f_ptr, "%*[^\n]\n");
         for(i = 0; i < model->L;i++){
                for(j = 0; j < model->num_states;j++){
                        fscanf(f_ptr,"%lf\n", &model->emission_counts[i][j]);
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
