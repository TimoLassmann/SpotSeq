#include "tldevel.h"

#include "tllogsum.h"

#include "tlhdf5wrap.h"

#include "finite_hmm_alloc.h"

#define BIASMODEL_IMPORT
#include "bias_model.h"


/* This is the bias model adopted from the HMMER 3.0 code. */

int build_bias_model(struct fhmm* model, struct fhmm** bias_model)
{
        struct fhmm* b = NULL;

        float M8;
        float sanity;
        int i;
        int j;

        RUNP(b = alloc_fhmm());

        b->K = 2;
        b->L = model->L;

        RUN(galloc(&b->e, b->K, b->L));

        RUN(galloc(&b->t, b->K, b->K));
        /* In Hmmer code  */
        M8 = (float)model->L / 1.0F;

        b->t[0][0] = prob2scaledprob(400.0F / 401.0F);
        b->t[0][1] = prob2scaledprob(  1.0F / 401.0F);
        b->t[1][1] = prob2scaledprob(    M8 / (M8 +1.0F));
        b->t[1][0] = prob2scaledprob(  1.0F / (M8 +1.0F));

        /* not used but allocates / initialised to prevent issues with hdf5 write */
        RUN( galloc(&b->background, b->L));
        RUN( galloc(&b->m_comp_back, b->L));
        RUN( galloc(&b->tindex, b->K+1,b->K+1));
        for(i = 0; i < b->L;i++){
                b->background[i] = model->background[i];
                b->m_comp_back[i] = model->m_comp_back[i];
        }

        for(i = 0; i <= b->K;i++){
                for(j = 0; j <= b->K;j++){
                        b->tindex[i][j] = 0;
                }
        }

        for(i = 0; i < b->L;i++){
                b->e[0][i] = prob2scaledprob(model->background[i] /model->background[i]);
                b->e[1][i] = prob2scaledprob(model->m_comp_back[i]/ model->background[i]);
        }

        sanity = prob2scaledprob(0.0F);
        for(i = 0; i < b->L;i++){
                sanity = logsum(sanity, prob2scaledprob(model->background[i]));
        }
        LOG_MSG("State 1 emission sum: %f", scaledprob2prob (sanity));
        sanity = prob2scaledprob(0.0F);
        for(i = 0; i < b->L;i++){
                sanity = logsum(sanity, prob2scaledprob(model->m_comp_back[i]));
        }
        LOG_MSG("State 1 emission sum: %f", scaledprob2prob(sanity));

        *bias_model = b;
        return OK;
ERROR:
        free_fhmm(b);
        return FAIL;
}




int write_biashmm(char* filename, struct fhmm* fhmm)
{

        struct hdf5_data* hdf5_data = NULL;
        float** tmp = NULL;
        int**tmp_int = NULL;
        int i,j;
        /* Create File  */
        RUN(open_hdf5_file(&hdf5_data,filename));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,"/biasfhmm","H",fhmm->H));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,"/biasfhmm","lambda",fhmm->lambda));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,"/biasfhmm","tau",fhmm->tau));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,"/biasfhmm","K",fhmm->K));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,"/biasfhmm","alloc_K",fhmm->K));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data ,"/biasfhmm","L",fhmm->L));

        /* transition */
        RUN(galloc(&tmp, fhmm->K,fhmm->K));
        for(i = 0; i < fhmm->K;i++){
                for(j = 0; j < fhmm->K;j++){
                        tmp[i][j] = fhmm->t[i][j];
                }
        }
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/biasfhmm","transition",tmp));
        gfree(tmp);
        tmp= NULL;

        /* emission */
        RUN(galloc(&tmp, fhmm->K,fhmm->L));
        for(i = 0; i < fhmm->K;i++){
                for(j = 0; j < fhmm->L;j++){
                        tmp[i][j] = fhmm->e[i][j];
                }
        }

        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/biasfhmm","emission",tmp));
        gfree(tmp);
        tmp= NULL;

        /* t-index */
        RUN(galloc(&tmp_int, fhmm->K,fhmm->K+1));
        for(i = 0; i < fhmm->K;i++){
                for(j = 0; j < fhmm->K+1;j++){
                        tmp_int[i][j] = fhmm->tindex[i][j];
                }
        }

        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/biasfhmm","transition_index",tmp_int));
        gfree(tmp);
        tmp= NULL;

        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/biasfhmm","background",fhmm->background));
        close_hdf5_file(&hdf5_data);

        return OK;
ERROR:
        if(hdf5_data){
                close_hdf5_file(&hdf5_data);

        }
        return FAIL;
}

int read_biasfhmm(char* filename,struct fhmm** ret)
{

        struct hdf5_data* hdf5_data = NULL;
        struct fhmm* fhmm = NULL;
        /* Create File  */
        RUN(open_hdf5_file(&hdf5_data,filename));



        ASSERT(hdf5_data != NULL, "No filename");

        RUNP(fhmm = alloc_fhmm());


        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data ,"/biasfhmm","H",&fhmm->H));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data ,"/biasfhmm","lambda",&fhmm->lambda));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data ,"/biasfhmm","tau",&fhmm->tau));

        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, "/biasfhmm", "K", &fhmm->K));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data ,"/biasfhmm","alloc_K",&fhmm->alloc_K));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, "/biasfhmm", "L", &fhmm->L));

        RUN(HDFWRAP_READ_DATA(hdf5_data, "/biasfhmm", "emission", &fhmm->e));

        RUN(HDFWRAP_READ_DATA(hdf5_data, "/biasfhmm", "transition", &fhmm->t));

        RUN(HDFWRAP_READ_DATA(hdf5_data, "/biasfhmm", "transition_index", &fhmm->tindex));

        RUN(HDFWRAP_READ_DATA(hdf5_data, "/biasfhmm", "background", &fhmm->background));


        close_hdf5_file(&hdf5_data);

        *ret = fhmm;
        return OK;
ERROR:
        if(hdf5_data){
                close_hdf5_file(&hdf5_data);

        }
        return FAIL;
}
