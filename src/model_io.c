
#include "tldevel.h"
//#include "model.h"
#include "model_struct.h"
#include "model_io.h"
#include "tlmisc.h"
#include <string.h>

#include "tlhdf5wrap.h"
#include "model_alloc.h"

#include "randomkit_io.h"

#include "finite_hmm_io.h"

#define BUFFER_LEN 256


struct ihmm_model* read_model_hdf5(struct hdf5_data* hdf5_data,char* group);
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

        //hdf5_data = hdf5_create();

        RUN(open_hdf5_file(&hdf5_data,filename));

        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/","MaxStates",&bag->max_num_states));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/","Number of models",&bag->num_models));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/","Seed",&bag->seed));
        ASSERT(bag->num_models > 0, "No models!");
        MMALLOC(bag->models, sizeof(struct ihmm_model*)* bag->num_models);


        MMALLOC(bag->finite_models, sizeof(struct fhmm*)* bag->num_models);
        MMALLOC(bag->min_u , sizeof(double) * bag->num_models);

        //get_group_names(hdf5_data);
        //fprintf(stdout,"Groups:\n");
        //for(i = 0; i < hdf5_data->grp_names->num_names;i++){
        //        fprintf(stdout,"%d %s\n",i,hdf5_data->grp_names->names[i]);
        //}
        RUN(HDFWRAP_READ_DATA(hdf5_data,"/","BestModel",&bag->best_model));

        /*
        RUN(hdf5_open_group("/",hdf5_data ));
        if((hdf5_data->dataset = H5Dopen(hdf5_data->group, "BestModel",H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        //printf ("H5Dopen returns: %d\n", hdf5_data->dataset);
        hdf5_read_attributes(hdf5_data,hdf5_data->dataset);
        hdf5_data->datatype  = H5Dget_type(hdf5_data->dataset );
        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        hdf5_data->rank      = H5Sget_simple_extent_ndims(hdf5_data->dataspace);
        hdf5_data->status  = H5Sget_simple_extent_dims(hdf5_data->dataspace,hdf5_data->dim , NULL);
        hdf5_data->status = H5Dread(hdf5_data->dataset, hdf5_data->datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bag->best_model);
        if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");
        RUN(hdf5_close_group(hdf5_data));*/

        RUN(read_RNG_state(hdf5_data, "/RNG",&bag->rndstate));

        for(i = 0; i < bag->num_models;i++){
                bag->models[i] = NULL;
                snprintf(buffer, BUFFER_LEN, "/models/m%d",i+1);
                RUNP(bag->models[i] = read_model_hdf5(hdf5_data, buffer));
                //if(bag->max_num_states < bag->models[i]->num_states){
                //      bag->max_num_states = bag->models[i]->num_states;
                //}
                //LOG_MSG("MAX: %d %d\n", bag->max_num_states,bag->models[i]->num_states);
                bag->finite_models[i] = NULL;
                //snprintf(buffer, BUFFER_LEN, "/models/m%d/fhmm",i+1);
                /* get HMM parameters  */
                //RUNP(bag->finite_models[i] = read_fhmm_parameters(hdf5_data, buffer));
        }
        RUN(close_hdf5_file(&hdf5_data));
        return bag;
ERROR:
        return NULL;
}


int write_best_model(char* filename, int best_model)
{
        struct hdf5_data* hdf5_data = NULL;

        //void *ptr;
        ASSERT(filename != NULL, "No filename");
        ASSERT(my_file_exists(filename) != 0,"File %s does not exist.",filename);

        RUN(open_hdf5_file(&hdf5_data,filename));

        RUN(HDFWRAP_WRITE_DATA(hdf5_data,"/","BestModel",best_model));
        RUN(close_hdf5_file(&hdf5_data));
        /*hdf5_data->rank = 1;
        hdf5_data->dim[0] = 1;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->native_type = H5T_NATIVE_INT;
        ptr = (void*) &best_model;
        RUN(hdf5_open_group("/",hdf5_data));
        if((hdf5_data->dataset = H5Dopen(hdf5_data->group,"BestModel",H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        if((hdf5_data->status  = H5Dwrite(hdf5_data->dataset,hdf5_data->native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr)) < 0) ERROR_MSG("H5Dwrite failed");
        RUN(hdf5_close_group(hdf5_data));

        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);*/
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
        /* Create File  */
        RUN(open_hdf5_file(&hdf5_data,filename));

        /* write top level meta-data */
        snprintf(buffer, BUFFER_LEN, "%s",PACKAGE_NAME );
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/","Program",buffer));
        snprintf(buffer, BUFFER_LEN, "%s", PACKAGE_VERSION);
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/","Version",buffer));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/","Number of models",bag->num_models));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/","MaxStates",bag->max_num_states ));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/","Seed",bag->seed));


        /*
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
        */
        //LOG_MSG("header done");
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/","BestModel",best_model ));


        //RUN(hdf5_create_group("/",hdf5_data));
        //hdf5_data->rank = 1;
        //hdf5_data->dim[0] = 1;
        //hdf5_data->chunk_dim[0] = 1;
        //hdf5_data->native_type = H5T_NATIVE_INT;
        //RUN(hdf5_write("BestModel",&best_model, hdf5_data));
        //RUN(hdf5_close_group(hdf5_data));

        RUN(add_RNG_state(hdf5_data, "/RNG",&bag->rndstate));
        //LOG_MSG("RNG top  done");
        //RUN(hdf5_create_group("models",hdf5_data));

        //RUN(hdf5_close_group(hdf5_data));
        //LOG_MSG("Create model group done.");

        for(model_i = 0; model_i < bag->num_models;model_i++){
                snprintf(buffer, BUFFER_LEN, "/models/m%d", model_i+1);
                //LOG_MSG("Writing model %d.",model_i);
                RUN(write_model_hdf5(hdf5_data, bag->models[model_i], buffer));

                //snprintf(buffer, BUFFER_LEN, "/models/m%d/fhmm", model_i+1);
                //RUN(add_fhmm(hdf5_data, bag->finite_models[model_i] , buffer));
                //LOG_MSG("Done");
        }
        close_hdf5_file(&hdf5_data);
        //hdf5_close_file(hdf5_data);
        //hdf5_free(hdf5_data);


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
                close_hdf5_file(&hdf5_data);

        }
        return FAIL;
}


struct ihmm_model* read_model_hdf5(struct hdf5_data* hdf5_data,char* group)
{
        char buffer[BUFFER_LEN+5];
        struct ihmm_model* model = NULL;
        int a,b,c;
        //int i;
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

        //hdf5_open_group(group,hdf5_data);
        //hdf5_read_attributes(hdf5_data, hdf5_data->group);
        //ASSERT(hdf5_data->num_attr != 0 , "Could not find attributes");
        //print_attributes(hdf5_data);
        a = 0;
        b = 0;
        c = 0;
        //LOG_MSG("Reading from: %s", group);
        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "Numberofstates", &a);
        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "Numberofletters", &b);
        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "MaxNumberofstates", &c);
        //LOG_MSG("K:%d L:%d", a,b);
        ASSERT(a!=0, "No states???");
        ASSERT(b!=0, "No letters???");
        ASSERT(c!=0, "No max num states???");
        int maxK = c;
        RUNP(model = alloc_ihmm_model(maxK, b,0));
        gfree(model->emission_counts);
        gfree(model->transition_counts);
        gfree(model->beta);

        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "Numberofstates", &model->num_states);
        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "Numberofletters", &model->L);
        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "Gamma", &model->gamma);
        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "gamma_a", &model->gamma_a);
        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "gamma_b", &model->gamma_b);

        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "Alpha", &model->alpha);
        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "alpha_a", &model->alpha_a);
        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "alpha_b", &model->alpha_b);

        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "Iteration", &model->training_iterations);
        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "Seed", &model->seed);
        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "Alpha_limit", &model->alpha_limit);
        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "Gamma_limit", &model->gamma_limit);



        //hdf5_close_file(hdf5_data);
        //hdf5_free(hdf5_data);
        /* stretch matrices... */

        //RUN(resize_ihmm_model(model, model->num_states+1));

        RUN(HDFWRAP_READ_DATA(hdf5_data, group, "Beta", &model->beta));
        RUN(HDFWRAP_READ_DATA(hdf5_data ,group,"Background",&model->background));
        RUN(HDFWRAP_READ_DATA(hdf5_data, group, "transition_counts", &model->transition_counts));
        RUN(HDFWRAP_READ_DATA(hdf5_data, group, "emission_counts", &model->emission_counts));
        snprintf(buffer, BUFFER_LEN+5 , "%s/RNG", group);
        //LOG_MSG("Trying to create group: %s", buffer);
        RUN(read_RNG_state(hdf5_data, buffer,&model->rndstate));
        //WARNING_MSG("Each time a run is continued a new RNG seed is selected...");
        RUN(get_dim1(model->beta, &model->alloc_num_states));
        return model;
ERROR:
        return NULL;
}



int write_model_hdf5(struct hdf5_data* hdf5_data,struct ihmm_model* model, char* group)
{
        //struct hdf5_data* hdf5_data = NULL;
        char buffer[BUFFER_LEN+5];
        //double** tmp = NULL;
        //int i,j;
        int best_model = -1;
        //RUNP(hdf5_data = hdf5_create());
        //snprintf(buffer, BUFFER_LEN, "%s",PACKAGE_NAME );
        //hdf5_add_attribute(hdf5_data, "Program", buffer, 0, 0.0f, HDF5GLUE_CHAR);
        //snprintf(buffer, BUFFER_LEN, "%s", PACKAGE_VERSION);
        //hdf5_add_attribute(hdf5_data, "Version", buffer, 0, 0.0f, HDF5GLUE_CHAR);


        //RUN(hdf5_create_file(filename,hdf5_data));

        //hdf5_write_attributes(hdf5_data, hdf5_data->file);
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"BestModel",best_model ));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"Numberofstates", model->num_states));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"MaxNumberofstates", model->alloc_num_states));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"Numberofletters",model->L));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"Gamma", model->gamma));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"Gamma_limit",model->gamma_limit));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"gamma_a",model->gamma_a));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"gamma_b",model->gamma_b));

        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"Alpha",model->alpha));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"Alpha_limit",model->alpha_limit));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"alpha_a",model->alpha_a));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"alpha_b",model->alpha_b));

        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"Iteration",model->training_iterations));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"Seed",model->seed));



        /*
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
        */
        /*hdf5_data->rank = 1;
        hdf5_data->dim[0] = model->num_states;
        hdf5_data->dim[1] = -1;
        hdf5_data->chunk_dim[0] = model->num_states;
        hdf5_data->chunk_dim[1] = -1;
        hdf5_data->native_type = H5T_NATIVE_DOUBLE;
        RUN(hdf5_write("Beta",&model->beta[0], hdf5_data));*/
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"Beta",model->beta));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"Background",model->background));
        //RUNP(tmp = galloc(tmp,  model->num_states,  model->num_states, 0.0));
        /* TODO: this may not be necessary as I use galloc to allocate model->transition_counts  */

        /*RUN(galloc(&tmp,  model->num_states,  model->num_states));

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
        tmp= NULL;*/
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"transition_counts",model->transition_counts));
        //RUNP(tmp = galloc(tmp,  model->L ,  model->num_states, 0.0));
        /*RUN(galloc(&tmp,  model->L ,  model->num_states));

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
        tmp= NULL;*/
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"emission_counts",model->emission_counts));

        /* Don't think I need to close group with new hdf api */
        //RUN(hdf5_close_group(hdf5_data));

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







