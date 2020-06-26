
#include "tldevel.h"
#include "tlalphabet.h"

#include "sequence_struct.h"

#define SEQUENCE_ALLOC_IMPORT
#include "sequence_alloc.h"

static void free_ihmm_sequence(struct ihmm_sequence* sequence);
static int alloc_multi_model_label_and_u(struct ihmm_sequence* sequence,int max_len, int num_models);

int alloc_sequence_buffer(struct seq_buffer** seq_buf,int num_seq)
{
        struct seq_buffer* sb = NULL;
        int i;
        MMALLOC(sb,sizeof(struct seq_buffer));
        sb->malloc_num = num_seq;
        sb->num_seq = -1;
        sb->sequences = NULL;
        sb->max_len = 0;
        sb->L = -1;
        sb->org_num_seq = -1;
        sb->alphabet = NULL;
        sb->num_state_arr = NULL;
        MMALLOC(sb->sequences, sizeof(struct ihmm_sequence*) *sb->malloc_num );
        for(i = 0; i < sb->malloc_num;i++){
                sb->sequences[i] = NULL;
                RUN(alloc_ihmm_seq(&sb->sequences[i]));
        }

        *seq_buf = sb;
        return OK;
ERROR:
        return FAIL;
}

int alloc_ihmm_seq(struct ihmm_sequence** is)
{
        struct ihmm_sequence* sequence = NULL;
        MMALLOC(sequence,sizeof(struct ihmm_sequence));
        sequence->has_path = NULL;
        sequence->seq = NULL;
        sequence->u = NULL;
        sequence->label = NULL;
        sequence->name = NULL;
        sequence->malloc_len = 128;
        sequence->seq_len = 0;
        sequence->score = 1.0f;
        sequence->r_score = 1.0f;
        MMALLOC(sequence->seq, sizeof(uint8_t) * sequence->malloc_len);
        MMALLOC(sequence->u, sizeof(double) * (sequence->malloc_len+1));
        MMALLOC(sequence->label , sizeof(int) * sequence->malloc_len);
        MMALLOC(sequence->name, sizeof(char) * MAX_SEQUENCE_NAME_LEN);
        sequence->label_arr = NULL;
        sequence->tmp_label_arr = NULL;
        sequence->u_arr = NULL;
        sequence->score_arr = NULL;

        *is = sequence;
        return OK;
ERROR:
        free_ihmm_sequence(sequence);
        return FAIL;
}

int realloc_ihmm_seq(struct ihmm_sequence* sequence, int new_len)
{
        ASSERT(sequence != NULL, "No Sequence.");

        while(sequence->malloc_len <= new_len){
                sequence->malloc_len = sequence->malloc_len + sequence->malloc_len / 2;
        }
        MREALLOC(sequence->seq, sizeof(uint8_t) *sequence->malloc_len);
        MREALLOC(sequence->u, sizeof(double) * (sequence->malloc_len+1));
        MREALLOC(sequence->label , sizeof(int) * sequence->malloc_len);

        return OK;
ERROR:
        return FAIL;
}


void free_ihmm_sequence(struct ihmm_sequence* sequence)
{
        if(sequence){
                if(sequence->has_path){
                        MFREE(sequence->has_path);
                }
                if(sequence->seq){
                        MFREE(sequence->seq);
                }
                if(sequence->name){
                        MFREE(sequence->name);
                }
                if(sequence->u){
                        MFREE(sequence->u);
                }
                if(sequence->label){
                        MFREE(sequence->label);
                }
                if(sequence->u_arr){
                        gfree(sequence->u_arr);
                }
                if(sequence->label_arr ){
                        gfree(sequence->label_arr);
                }

                if(sequence->tmp_label_arr ){
                        gfree(sequence->tmp_label_arr);
                }



                if(sequence->score_arr){
                        gfree(sequence->score_arr);
                }
                MFREE(sequence);
        }
}


void free_ihmm_sequences(struct seq_buffer* sb)
{
        int i;
        if(sb){
                for(i =0; i < sb->malloc_num;i++){
                        free_ihmm_sequence(sb->sequences[i]);
                }
                MFREE(sb->sequences);
                if(sb->alphabet){
                        free_alphabet(sb->alphabet);
                }
                if(sb->num_state_arr){
                        MFREE(sb->num_state_arr);
                }
                //if(sb->background){
                //gfree(sb->background);
                //}
                MFREE(sb);
        }
}


int add_multi_model_label_and_u(struct seq_buffer* sb,int num_models)
{
        int i;

        ASSERT(sb != NULL,"No sequence buffer");

        for(i = 0; i < sb->num_seq;i++){
                RUN(alloc_multi_model_label_and_u(sb->sequences[i],sb->max_len,   num_models));

        }
        return OK;
ERROR:
        return FAIL;
}

int alloc_multi_model_label_and_u(struct ihmm_sequence* sequence,int max_len, int num_models)
{
        int i,j;

        ASSERT(sequence != NULL, "No sequence");

        RUN(galloc(&sequence->u_arr, num_models, max_len+1));

        RUN(galloc(&sequence->label_arr, num_models, max_len+1));

        RUN(galloc(&sequence->tmp_label_arr, num_models, max_len+1));

        for(i =0; i < num_models;i++){
                for(j = 0; j < max_len+1;j++){
                        sequence->u_arr[i][j] = 0.0;
                        sequence->label_arr[i][j] = 0;
                        sequence->tmp_label_arr[i][j] = 0;
                }
        }

        RUN(galloc(&sequence->score_arr, num_models));

        for(i = 0; i < num_models;i++){
                sequence->score_arr[i] = 1.0; /* default weight is one.  */
        }

        MMALLOC(sequence->has_path,sizeof(uint8_t) * num_models);
        return OK;
ERROR:
        free_ihmm_sequence(sequence);
        return FAIL;
}
