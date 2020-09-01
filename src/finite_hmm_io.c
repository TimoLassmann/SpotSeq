#include "tldevel.h"

#include "model_struct.h"
#include "model_io.h"
#include "model_alloc.h"
#include "finite_hmm_struct.h"
#include "finite_hmm_alloc.h"

#define FINITE_HMM_IO_IMPORT
#include "finite_hmm_io.h"


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



int add_fhmm(struct hdf5_data* hdf5_data, struct fhmm* fhmm, char* group)
{


        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"H",fhmm->H));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"lambda",fhmm->lambda));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"tau",fhmm->tau));


        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"K",fhmm->K));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"alloc_K",fhmm->alloc_K));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"L",fhmm->L));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"emission",fhmm->e));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"transition",fhmm->t));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"transition_index",fhmm->tindex));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"background",fhmm->background));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"ModelCompoBack",fhmm->m_comp_back ));
        return OK;
ERROR:
        return FAIL;
}


struct fhmm*  read_fhmm_parameters(struct hdf5_data* hdf5_data, char* group)
{
        struct fhmm* fhmm = NULL;

        ASSERT(hdf5_data != NULL, "No filename");

        RUNP(fhmm = alloc_fhmm());


        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data ,group,"H",&fhmm->H));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data ,group,"lambda",&fhmm->lambda));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data ,group,"tau",&fhmm->tau));

        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "K", &fhmm->K));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data ,group,"alloc_K",&fhmm->alloc_K));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "L", &fhmm->L));

        RUN(HDFWRAP_READ_DATA(hdf5_data, group, "emission", &fhmm->e));

        RUN(HDFWRAP_READ_DATA(hdf5_data, group, "transition", &fhmm->t));

        RUN(HDFWRAP_READ_DATA(hdf5_data, group, "transition_index", &fhmm->tindex));

        RUN(HDFWRAP_READ_DATA(hdf5_data, group, "background", &fhmm->background));

        RUN(HDFWRAP_READ_DATA(hdf5_data ,group,"ModelCompoBack",&fhmm->m_comp_back ));


        return fhmm;
ERROR:
        return NULL;

}

int write_searchfhmm(char* filename, struct fhmm* fhmm)
{

        struct hdf5_data* hdf5_data = NULL;
        float** tmp = NULL;
        int**tmp_int = NULL;
        int i,j;
        /* Create File  */
        RUN(open_hdf5_file(&hdf5_data,filename));

        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,"/bestfhmm","H",fhmm->H));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,"/bestfhmm","lambda",fhmm->lambda));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,"/bestfhmm","tau",fhmm->tau));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,"/bestfhmm","K",fhmm->K));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,"/bestfhmm","alloc_K",fhmm->K));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,"/bestfhmm","L",fhmm->L));

        /* transition */
        RUN(galloc(&tmp, fhmm->K,fhmm->K));
        for(i = 0; i < fhmm->K;i++){
                for(j = 0; j < fhmm->K;j++){
                        tmp[i][j] = fhmm->t[i][j];
                }
        }
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/bestfhmm","transition",tmp));
        gfree(tmp);
        tmp= NULL;

        /* emission */
        RUN(galloc(&tmp, fhmm->K,fhmm->L));
        for(i = 0; i < fhmm->K;i++){
                for(j = 0; j < fhmm->L;j++){
                        tmp[i][j] = fhmm->e[i][j];
                }
        }

        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/bestfhmm","emission",tmp));
        gfree(tmp);
        tmp= NULL;

        /* t-index */
        RUN(galloc(&tmp_int, fhmm->K,fhmm->K+1));
        for(i = 0; i < fhmm->K;i++){
                for(j = 0; j < fhmm->K+1;j++){
                        tmp_int[i][j] = fhmm->tindex[i][j];
                }
        }

        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/bestfhmm","transition_index",tmp_int));
        gfree(tmp);
        tmp= NULL;


        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/bestfhmm","background",fhmm->background));
        close_hdf5_file(&hdf5_data);
        return OK;
ERROR:
        if(hdf5_data){
                close_hdf5_file(&hdf5_data);

        }
        return FAIL;
}

int read_searchfhmm(char* filename,struct fhmm** ret)
{

        struct hdf5_data* hdf5_data = NULL;
        struct fhmm* fhmm = NULL;
        /* Create File  */
        RUN(open_hdf5_file(&hdf5_data,filename));



        ASSERT(hdf5_data != NULL, "No filename");

        RUNP(fhmm = alloc_fhmm());


        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data ,"/bestfhmm","H",&fhmm->H));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data ,"/bestfhmm","lambda",&fhmm->lambda));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data ,"/bestfhmm","tau",&fhmm->tau));

        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, "/bestfhmm", "K", &fhmm->K));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data ,"/bestfhmm","alloc_K",&fhmm->alloc_K));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, "/bestfhmm", "L", &fhmm->L));

        RUN(HDFWRAP_READ_DATA(hdf5_data, "/bestfhmm", "emission", &fhmm->e));

        RUN(HDFWRAP_READ_DATA(hdf5_data, "/bestfhmm", "transition", &fhmm->t));

        RUN(HDFWRAP_READ_DATA(hdf5_data, "/bestfhmm", "transition_index", &fhmm->tindex));

        RUN(HDFWRAP_READ_DATA(hdf5_data, "/bestfhmm", "background", &fhmm->background));


        close_hdf5_file(&hdf5_data);

        *ret = fhmm;
        return OK;
ERROR:
        if(hdf5_data){
                close_hdf5_file(&hdf5_data);

        }
        return FAIL;
}
