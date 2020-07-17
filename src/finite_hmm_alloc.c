#include "tldevel.h"

#include "finite_hmm_struct.h"

#define FINITE_HMM_ALLOC_IMPORT
#include "finite_hmm_alloc.h"



int alloc_fhmm_dyn_mat(struct fhmm_dyn_mat** mat,int L,int K)
{
        struct fhmm_dyn_mat* dm = NULL;
        int i,j;
        MMALLOC(dm,sizeof(struct fhmm_dyn_mat));

        dm->F_matrix = NULL;
        dm->B_matrix = NULL;
        dm->F_NBECJ = NULL;
        dm->B_NBECJ = NULL;

        dm->alloc_K = K;
        dm->alloc_matrix_len = L;

        RUN(galloc(&dm->F_matrix, dm->alloc_matrix_len, dm->alloc_K));
        RUN(galloc(&dm->B_matrix, dm->alloc_matrix_len, dm->alloc_K));
        for(i = 0; i < dm->alloc_matrix_len;i++){
                for(j = 0;j < dm->alloc_K;j++){
                        dm->F_matrix[i][j] = 0.0;
                        dm->B_matrix[i][j] = 0.0;
                }
        }
        RUN(galloc(&dm->F_NBECJ, dm->alloc_matrix_len, 5));
        RUN(galloc(&dm->B_NBECJ, dm->alloc_matrix_len, 5));
        for(i = 0; i < dm->alloc_matrix_len;i++){
                for(j = 0;j < 5;j++){
                        dm->F_NBECJ[i][j] = 0.0;
                        dm->B_NBECJ[i][j] = 0.0;
                }
        }

        *mat = dm;
        return OK;
ERROR:
        return FAIL;
}

int resize_fhmm_dyn_mat(struct fhmm_dyn_mat* dm,int new_len)
{
        int i,j;
        ASSERT(dm != NULL, "No model");
        ASSERT(new_len > 0, "newlen has to be > 0");
        ASSERT(dm->alloc_matrix_len > 0, "No matrix allocated yet...");

        if(dm->alloc_matrix_len < new_len){
                while(dm->alloc_matrix_len < new_len){
                        dm->alloc_matrix_len = dm->alloc_matrix_len +  dm->alloc_matrix_len / 2;
                }
                RUN(galloc(&dm->F_matrix, dm->alloc_matrix_len, dm->alloc_K));
                RUN(galloc(&dm->B_matrix, dm->alloc_matrix_len, dm->alloc_K));
                for(i = 0; i < dm->alloc_matrix_len;i++){
                        for(j = 0;j < dm->alloc_K;j++){
                                dm->F_matrix[i][j] = 0.0;
                                dm->B_matrix[i][j] = 0.0;
                        }
                }
                RUN(galloc(&dm->F_NBECJ, dm->alloc_matrix_len, 5));
                RUN(galloc(&dm->B_NBECJ, dm->alloc_matrix_len, 5));
                for(i = 0; i < dm->alloc_matrix_len;i++){
                        for(j = 0;j < 5;j++){
                                dm->F_NBECJ[i][j] = 0.0;
                                dm->B_NBECJ[i][j] = 0.0;
                        }
                }

        }
        return OK;
ERROR:
        return FAIL;
}

int free_fhmm_dyn_mat(struct fhmm_dyn_mat* dm)
{
        if(dm){
                if(dm->F_matrix){
                        gfree(dm->F_matrix);
                }
                if(dm->B_matrix){
                        gfree(dm->B_matrix);
                }
                if(dm->F_NBECJ){
                        gfree(dm->F_NBECJ);
                }
                if(dm->B_NBECJ){
                        gfree(dm->B_NBECJ);
                }

                MFREE(dm);

        }

        return OK;
}


struct fhmm* alloc_fhmm(void)
{

        struct fhmm* fhmm = NULL;

        MMALLOC(fhmm, sizeof(struct fhmm));
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
