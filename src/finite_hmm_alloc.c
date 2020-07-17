#include "tldevel.h"

#include "finite_hmm_struct.h"

#define FINITE_HMM_ALLOC_IMPORT
#include "finite_hmm_alloc.h"




struct fhmm* alloc_fhmm(void)
{

        struct fhmm* fhmm = NULL;

        MMALLOC(fhmm, sizeof(struct fhmm));
        fhmm->F_matrix = NULL;
        fhmm->B_matrix = NULL;
        fhmm->F_NBECJ = NULL;
        fhmm->B_NBECJ = NULL;
        fhmm->e = NULL;
        fhmm->t = NULL;
        fhmm->tindex = NULL;
        fhmm->background = NULL;
        fhmm->K = 0;
        fhmm->L = 0;
        fhmm->f_score = 0.0;
        fhmm->b_score = 0.0;
        fhmm->r_score = 0.0;

        fhmm->alloc_matrix_len = 0;

        fhmm->alloc_K = 0;
        return fhmm;
ERROR:
        free_fhmm(fhmm);
        return NULL;
}

void free_fhmm(struct fhmm* fhmm)
{
        if(fhmm){
                if(fhmm->F_matrix){
                        gfree(fhmm->F_matrix);
                        //free_2d((void**) fhmm->F_matrix);
                }
                if(fhmm->B_matrix){
                        gfree(fhmm->B_matrix);
                        //free_2d((void**) fhmm->B_matrix);
                }

                if(fhmm->e){
                        gfree(fhmm->e);
                        //free_2d((void**) fhmm->e);
                }
                if(fhmm->t){
                        gfree(fhmm->t);
                        //free_2d((void**) fhmm->t);
                }
                if(fhmm->background){
                        gfree(fhmm->background);
                }
                if(fhmm->tindex){
                        gfree(fhmm->tindex);
                        //free_2d((void**)fhmm->tindex);
                }

                MFREE(fhmm);
        }
}
