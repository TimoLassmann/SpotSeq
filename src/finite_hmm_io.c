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
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"K",fhmm->K));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,group,"L",fhmm->L));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"emission",fhmm->e));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"transition",fhmm->t));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"transition_index",fhmm->tindex));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"background",fhmm->background));
        return OK;
ERROR:
        return FAIL;
}

struct fhmm*  read_fhmm_parameters(struct hdf5_data* hdf5_data, char* group)
{
        struct fhmm* fhmm = NULL;

        ASSERT(hdf5_data != NULL, "No filename");

        RUNP(fhmm = alloc_fhmm());

        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "K", &fhmm->K);

        HDFWRAP_READ_ATTRIBUTE(hdf5_data, group, "L", &fhmm->L);

        HDFWRAP_READ_DATA(hdf5_data, group, "emission", &fhmm->e);

        HDFWRAP_READ_DATA(hdf5_data, group, "transition", &fhmm->t);

        HDFWRAP_READ_DATA(hdf5_data, group, "transition_index", &fhmm->tindex);

        HDFWRAP_READ_DATA(hdf5_data, group, "background", &fhmm->background);


        return fhmm;
ERROR:
        return NULL;

}
