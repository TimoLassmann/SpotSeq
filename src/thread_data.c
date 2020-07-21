#include "tldevel.h"
#include "randomkit.h"
#include "randomkit_tl_add.h"

#include "finite_hmm_alloc.h"

#define THREAD_DATA_IMPORT
#include "thread_data.h"

//#include "global.h"



struct seqer_thread_data** create_seqer_thread_data(int* num_threads, int max_len, int K,rk_state* random,int mode)
{
        struct seqer_thread_data** td = NULL;
        int i,j,c;
        int local_num_treads;

        ASSERT(*num_threads> 0, "no threads");
        local_num_treads = *num_threads;

        MMALLOC(td, sizeof(struct seqer_thread_data*) * local_num_treads);
        for(i = 0; i < local_num_treads;i++){
                td[i] = NULL;
                MMALLOC(td[i], sizeof(struct seqer_thread_data));
                td[i]->dyn = NULL;
                td[i]->fhmm = NULL;
                RUN(alloc_fhmm_dyn_mat(&td[i]->fmat, max_len, K));
                //td[i]->F_matrix = NULL;
                //td[i]->B_matrix = NULL;
                //td[i]->t = NULL;
                //td[i]->e = NULL;

                //  RUNP(matrix = malloc_2d_float(matrix,sb->max_len+1, ft->last_state, 0.0f));

                RUN(galloc(&td[i]->dyn, max_len, K));

                //RUN(galloc(&td[i]->F_matrix, max_len, K));
                //RUN(galloc(&td[i]->B_matrix, max_len, K));
                for(j = 0; j < max_len;j++){
                        for(c = 0; c < K;c++){
                                td[i]->dyn[j][c] = 0.0;
                                //td[i]->F_matrix[j][c] = 0.0;
                                //td[i]->B_matrix[j][c] = 0.0;
                        }
                }

                /* was initailized to -INFINITY  */
                /*if(mode == THREAD_DATA_FULL){
                        RUN(galloc(&td[i]->t, K,K));
                        RUN(galloc(&td[i]->e,ALPHABET_PROTEIN,K));
                        for(j = 0; j < K;j++){
                                for(c = 0; c < K;c++){
                                        td[i]->t[j][c] = -INFINITY;
                                }
                        }
                        for(j = 0; j < ALPHABET_PROTEIN;j++){
                                for(c = 0; c < K;c++){
                                        td[i]->e[j][c] = -INFINITY;
                                }
                        }
                        RUN(galloc(&td[i]->F_matrix, max_len, K));
                        RUN(galloc(&td[i]->B_matrix, max_len, K));
                        for(j = 0; j < max_len;j++){
                                for(c = 0; c < K;c++){
                                        td[i]->F_matrix[j][c] = 0.0;
                                        td[i]->B_matrix[j][c] = 0.0;
                                }
                        }
                }*/


                td[i]->ft = NULL;
                td[i]->sb = NULL;
                td[i]->thread_ID = i;
                td[i]->num_threads = local_num_treads;
                if(random){
                        td[i]->seed =  rk_ulong(random);
                        rk_seed(td[i]->seed, &td[i]->rndstate);
                }else{
                        td[i]->seed = 42; /* placeholder! if no RNG
                                           * given as argument a RNG
                                           * statewill be readfrom
                                           * hdf5 */
                }
                //fprintf(stdout,"thread:%d seed: %d\n",i, td[i]->seed);
        }

        *num_threads = local_num_treads;
        return td;
ERROR:
        free_seqer_thread_data(td);
        return NULL;
}


int resize_seqer_thread_data(struct seqer_thread_data** td,int* num_threads, int max_len, int K)
{
        int i,j,c;
        int local_num_treads;
        int cur_threads;

        ASSERT(*num_threads> 0, "no threads");
        local_num_treads = *num_threads;
        cur_threads =  *num_threads;



        for(i = local_num_treads; i < cur_threads;i++){
                gfree(td[i]->dyn);
                //free_2d((void**) td[i]->dyn);
                MFREE(td[i]);
        }

        //LOG_MSG("mallocing auxiliary datastructures to %d %d", max_len,K);
        for(i = 0; i < local_num_treads;i++){
                RUN(galloc(&td[i]->dyn, max_len, K));
                LOG_MSG("Alloc: %d %d", max_len,K);
                RUN(resize_fhmm_dyn_mat(td[i]->fmat , max_len, K));
                //RUN(galloc(&td[i]->F_matrix, max_len, K));
                //RUN(galloc(&td[i]->B_matrix, max_len, K));
                for(j = 0; j < max_len;j++){
                        for(c = 0; c < K;c++){
                                td[i]->dyn[j][c] = 0.0;
                                //td[i]->F_matrix[j][c] = 0.0;
                                //td[i]->B_matrix[j][c] = 0.0;
                        }
                }
                /*RUN(galloc(&td[i]->t, K,K));
                RUN(galloc(&td[i]->e,ALPHABET_PROTEIN,K));
                for(j = 0; j < K;j++){
                        for(c = 0; c < K;c++){
                                td[i]->t[j][c] = -INFINITY;
                        }
                }
                for(j = 0; j < ALPHABET_PROTEIN;j++){
                        for(c = 0; c < K;c++){
                                td[i]->e[j][c] = -INFINITY;
                        }
                }*/
                td[i]->num_threads = local_num_treads;
        }
        *num_threads = local_num_treads;
        return OK;
ERROR:
        return FAIL;
}

int compare_wims_data(struct seqer_thread_data** a , struct seqer_thread_data** b, int num)
{
        int i;
        for(i = 0; i < num;i++){
                ASSERT(a[i]->seed == b[i]->seed,"Seeds differ: %d %d",a[i]->seed,b[i]->seed);
                //ASSERT(a[i]->num_seq == b[i]->num_seq,"Num_Seqs differ: %d %d",a[i]->num_seq,b[i]->num_seq);
                ASSERT(a[i]->thread_ID == b[i]->thread_ID,"Thread_IDs differ: %d %d",a[i]->thread_ID,b[i]->thread_ID);
                ASSERT(a[i]->num_threads == b[i]->num_threads,"Num_Threadss differ: %d %d",a[i]->num_threads,b[i]->num_threads);

                compare_rk_state(&a[i]->rndstate,&b[i]->rndstate);
        }
        return OK;
ERROR:
        return FAIL;
}

void free_seqer_thread_data(struct seqer_thread_data** td)
{
        int i;
        if(td){
                int num_threads = td[0]->num_threads;
                for(i = 0; i < num_threads;i++){
                        gfree(td[i]->dyn);
                        free_fhmm_dyn_mat(td[i]->fmat);
                        //gfree(td[i]->F_matrix);
                        //gfree(td[i]->B_matrix);
                        //gfree(td[i]->t);
                        //gfree(td[i]->e);
                        MFREE(td[i]);
                }
                MFREE(td);
        }
}
