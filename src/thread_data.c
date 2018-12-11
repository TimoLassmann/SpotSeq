

#include "thread_data.h"

#include "global.h"




struct spotseq_thread_data** create_spotseq_thread_data(int* num_threads, int max_len, int K,rk_state* random)
{
        struct spotseq_thread_data** td = NULL;
        int i;
        int local_num_treads;
        size_t mem_needed;
        ASSERT(*num_threads> 0, "no threads");
        local_num_treads = *num_threads;


        mem_needed = sizeof(double) * local_num_treads * max_len * K;

        while(mem_needed  > GB){
                local_num_treads--;
                ASSERT(local_num_treads != 0, "No space! %d asked for but the limit is %d",mem_needed,GB);
                mem_needed = sizeof(double) * local_num_treads * max_len * K;

        }

        MMALLOC(td, sizeof(struct spotseq_thread_data*) * local_num_treads);
        for(i = 0; i < local_num_treads;i++){
                td[i] = NULL;
                MMALLOC(td[i], sizeof(struct spotseq_thread_data));
                td[i]->dyn = NULL;
                td[i]->F_matrix = NULL;
                td[i]->B_matrix = NULL;
                td[i]->t = NULL;
                td[i]->e = NULL;
                td[i]->fhmm = NULL;
                //  RUNP(matrix = malloc_2d_float(matrix,sb->max_len+1, ft->last_state, 0.0f));
                RUNP(td[i]->dyn = galloc(td[i]->dyn, max_len, K, 0.0));
                RUNP(td[i]->F_matrix = galloc(td[i]->F_matrix, max_len, K, 0.0));
                RUNP(td[i]->B_matrix = galloc(td[i]->B_matrix, max_len, K, 0.0));

                RUNP(td[i]->t = galloc(td[i]->t, K,K, -INFINITY));
                RUNP(td[i]->e = galloc(td[i]->e,ALPHABET_PROTEIN,K,-INFINITY));


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
        free_spotseq_thread_data(td, *num_threads);
        return NULL;
}


int resize_spotseq_thread_data(struct spotseq_thread_data** td,int* num_threads, int max_len, int K)
{
        int i;
        int local_num_treads;
        int cur_threads;
        size_t mem_needed;
        ASSERT(*num_threads> 0, "no threads");
        local_num_treads = *num_threads;
        cur_threads =  *num_threads;

        mem_needed = sizeof(double) * local_num_treads * max_len * K;

        while(mem_needed  > GB){
                local_num_treads--;
                ASSERT(local_num_treads != 0, "No space! %d asked for but the limit is %d",mem_needed,GB);
                mem_needed = sizeof(double) * local_num_treads * max_len * K;
        }


        for(i = local_num_treads; i < cur_threads;i++){

                free_2d((void**) td[i]->dyn);
                MFREE(td[i]);
        }

        //LOG_MSG("mallocing auxiliary datastructures to %d %d", max_len,K);
        for(i = 0; i < local_num_treads;i++){
                RUNP(td[i]->dyn = galloc(td[i]->dyn, max_len, K, 0.0));

                RUNP(td[i]->F_matrix = galloc(td[i]->F_matrix, max_len, K, 0.0));
                RUNP(td[i]->B_matrix = galloc(td[i]->B_matrix, max_len, K, 0.0));

                RUNP(td[i]->t = galloc(td[i]->t, K,K, -INFINITY));
                RUNP(td[i]->e = galloc(td[i]->e,ALPHABET_PROTEIN,K,-INFINITY));

                td[i]->num_threads = local_num_treads;
        }
        *num_threads = local_num_treads;
        return OK;
ERROR:
        return FAIL;
}


void free_spotseq_thread_data(struct spotseq_thread_data** td, int num_threads)
{
        int i;
        if(td){
                for(i = 0; i < num_threads;i++){
                        gfree(td[i]->dyn);
                        gfree(td[i]->F_matrix);
                        gfree(td[i]->B_matrix);
                        gfree(td[i]->t);
                        gfree(td[i]->e);
                        MFREE(td[i]);
                }
                MFREE(td);
        }
}
