

#include "thread_data.h"

#include "global.h"



struct spotseq_thread_data** create_spotseq_thread_data(int* num_threads, int max_len, int K)
{
        struct spotseq_thread_data** td = NULL;
        int i;
        int local_num_treads;
        size_t mem_needed;
        ASSERT(*num_threads> 0, "no threads");
        local_num_treads = *num_threads;


        mem_needed = sizeof(float) * local_num_treads * max_len * K;

        while(mem_needed  > GB){
                local_num_treads--;
                ASSERT(local_num_treads != 0, "No space! %d asked for but the limit is %d",mem_needed,GB);
                mem_needed = sizeof(float) * local_num_treads * max_len * K;

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
                RUNP(td[i]->dyn = malloc_2d_float(td[i]->dyn, max_len, K, 0.0f));
                RUNP(td[i]->F_matrix = malloc_2d_float(td[i]->F_matrix, max_len, K, 0.0f));
                RUNP(td[i]->B_matrix = malloc_2d_float(td[i]->B_matrix, max_len, K, 0.0f));

                RUNP(td[i]->t = malloc_2d_float(td[i]->t, K,K, -INFINITY));
                RUNP(td[i]->e = malloc_2d_float(td[i]->e,ALPHABET_PROTEIN,K,-INFINITY));


                td[i]->ft = NULL;
                td[i]->sb = NULL;
                td[i]->thread_ID = i;
                td[i]->num_threads = local_num_treads;
                td[i]->seed =  time(NULL) * (i+1);


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

        mem_needed = sizeof(float) * local_num_treads * max_len * K;

        while(mem_needed  > GB){
                local_num_treads--;
                ASSERT(local_num_treads != 0, "No space! %d asked for but the limit is %d",mem_needed,GB);
                mem_needed = sizeof(float) * local_num_treads * max_len * K;
        }


        for(i = local_num_treads; i < cur_threads;i++){

                free_2d((void**) td[i]->dyn);
                MFREE(td[i]);
        }

        //LOG_MSG("mallocing auxiliary datastructures to %d %d", max_len,K);
        for(i = 0; i < local_num_treads;i++){
                RUNP(td[i]->dyn = malloc_2d_float(td[i]->dyn, max_len, K, 0.0f));

                RUNP(td[i]->F_matrix = malloc_2d_float(td[i]->F_matrix, max_len, K, 0.0f));
                RUNP(td[i]->B_matrix = malloc_2d_float(td[i]->B_matrix, max_len, K, 0.0f));

                RUNP(td[i]->t = malloc_2d_float(td[i]->t, K,K, -INFINITY));
                RUNP(td[i]->e = malloc_2d_float(td[i]->e,ALPHABET_PROTEIN,K,-INFINITY));

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
                        free_2d((void**) td[i]->dyn);
                        free_2d((void**) td[i]->F_matrix);
                        free_2d((void**) td[i]->B_matrix);
                        free_2d((void**) td[i]->t);
                        free_2d((void**) td[i]->e);
                        MFREE(td[i]);
                }
                MFREE(td);
        }
}
