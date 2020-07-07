#include "tldevel.h"
#include "tlhdf5wrap.h"
#include "pst_structs.h"

#define PST_IO_IMPORT
#include "pst_io.h"

int write_pst_hdf5(struct pst* p, char* filename)
{
        struct hdf5_data* hdf5_data = NULL;

        RUN(open_hdf5_file(&hdf5_data,filename));

        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/PstModel","p_min",p->p_min));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/PstModel","gamma",p->gamma_min));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/PstModel","a",p->a));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/PstModel","b",p->b));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/PstModel","var",p->var));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/PstModel","L",p->L ));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/PstModel","len",p->len ));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/PstModel","Size",p->fpst_root->l ));
        //LOG_MSG("p->max_ob: %d", p->max_observed_len);
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/PstModel","MaxObsLen",p->max_observed_len));
        //LOG_MSG("p->fit: %p", p->fit);
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/PstModel","Fit", p->fit));
        //LOG_MSG("p->fit: %p", p->fit_index);
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/PstModel","FitIndex", p->fit_index));


        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/PstModel","Background", p->background));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/PstModel","LogBackground", p->lbg ));

        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/PstModel","Links", p->fpst_root->links));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/PstModel","Probabilities", p->fpst_root->prob));

        close_hdf5_file(&hdf5_data);

        return OK;
ERROR:
        if(hdf5_data){
                close_hdf5_file(&hdf5_data);

        }
        return FAIL;
}

int read_pst_hdf5(struct pst** pst, char* filename)
{
        struct hdf5_data* hdf5_data = NULL;
        struct pst* p = NULL;
        struct fpst*f = NULL;

        MMALLOC(p, sizeof(struct pst));

        p->L = 0;
        p->background = NULL;
        p->fpst_root = NULL;
        p->gamma_min = 0.0f;
        p->p_min = 0.0f;
        p->len = 0;
        p->lbg = NULL;

        MMALLOC(f, sizeof(struct fpst));
        f->prob = NULL;
        f->links = NULL;
        f->m = 0;
        f->l= 0;
        p->fpst_root = f;

        RUN(open_hdf5_file(&hdf5_data,filename));

        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/PstModel","p_min",&p->p_min));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/PstModel","gamma",&p->gamma_min));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/PstModel","a",&p->a));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/PstModel","b",&p->b));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/PstModel","var",&p->var));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/PstModel","MaxObsLen",&p->max_observed_len));

        RUN(HDFWRAP_READ_DATA(hdf5_data ,"/PstModel","Fit", &p->fit));
        RUN(HDFWRAP_READ_DATA(hdf5_data ,"/PstModel","FitIndex", &p->fit_index));

        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/PstModel","L",&p->L ));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/PstModel","len",&p->len ));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/PstModel","Size",&p->fpst_root->l ));

        RUN(HDFWRAP_READ_DATA(hdf5_data ,"/PstModel","Background", &p->background));
        RUN(HDFWRAP_READ_DATA(hdf5_data ,"/PstModel","LogBackground", &p->lbg ));

        RUN(HDFWRAP_READ_DATA(hdf5_data ,"/PstModel","Links", &p->fpst_root->links));
        RUN(HDFWRAP_READ_DATA(hdf5_data ,"/PstModel","Probabilities", &p->fpst_root->prob));


        close_hdf5_file(&hdf5_data);

        *pst = p;
        return OK;
ERROR:
        if(hdf5_data){
                close_hdf5_file(&hdf5_data);

        }

        if(p){
                if(p->fpst_root){
                        MFREE(p->fpst_root);
                }
                MFREE(p);
        }
        return FAIL;
}
