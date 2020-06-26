#include "tldevel.h"
#include "tlalphabet.h"
#include "tlseqbuffer.h"
#include "tlrng.h"

#include "global.h"             /* needed for MAX_NUM_STATES */

#include "sequence_struct.h"
#include "sequence_alloc.h"

#define SEQUENCE_PREP_IMPORT
#include "sequence_prep.h"


static int init_labelling(struct seq_buffer* sb, struct rng_state* rng, int num_models,int num_states, double sigma);
static int random_label_sequences(struct seq_buffer* sb, int K, int model_index, struct rng_state* rng);

static int set_alphabet_and_convert(struct seq_buffer* sb, struct rng_state* rng);

/* This function should prepare sequences to be used in beam sampleing
   This involves:
   - adding an alphabet (and sanity checks !!! )
   - convert sequences from char to numbered
   - add multi_u etc arrays
   - initial random labelling

 */
int prep_sequences(struct seq_buffer* sb, struct rng_state* rng, int num_models,int num_states, double sigma)
{

        ASSERT(sb != NULL, "No sequence buffer");
        ASSERT(num_models > 0, "No num_models");

        /* alphabet */
        RUN(set_alphabet_and_convert(sb,rng));

        /* add u and label */
        RUN(add_multi_model_label_and_u(sb, num_models));

        /* randomly label sequences  */
        /* I know the number of models ; unknown are number of states  */
        RUN(init_labelling(sb, rng, num_models, num_states, sigma));

        return OK;
ERROR:
        return FAIL;
}

int init_labelling(struct seq_buffer* sb, struct rng_state* rng, int num_models,int num_states, double sigma)
{
        int* num_state_arr = NULL;

        double average_sequence_len = 0.0;

        int i;

        if(sigma == 0){
                sigma = 1.0;
        }
        MMALLOC(num_state_arr, sizeof(int) * num_models);

        if(!num_states){
                average_sequence_len = 0.0;
                for(i = 0; i < sb->num_seq;i++){
                        average_sequence_len += sb->sequences[i]->seq_len;
                }

                average_sequence_len /= (double) sb->num_seq;
                average_sequence_len *= 0.5;
                num_states = round(average_sequence_len);
        }
        for(i = 0; i < num_models;i++){
                num_state_arr[i] = tl_random_gaussian(rng, num_states,sigma);
                num_state_arr[i] = MACRO_MIN(num_state_arr[i], MAX_NUM_STATES-2);
                num_state_arr[i] = MACRO_MAX(num_state_arr[i], 10);

                num_state_arr[i] = 16;
        }

        for(i = 0; i < num_models;i++){
                random_label_sequences(sb, num_state_arr[i], i, rng);
        }


        MFREE(num_state_arr);

        return OK;
ERROR:
        return FAIL;
}

int random_label_sequences(struct seq_buffer* sb, int K, int model_index, struct rng_state* rng)
{
        int i,j;
        uint16_t* label;
        int len;
        int min,max;
        ASSERT(sb != NULL, "No sequences");

        min = 999999;
        max = -1;
        for(i = 0;i< sb->num_seq;i++){
                label = sb->sequences[i]->label_arr[model_index];
                len = sb->sequences[i]->seq_len;
                for(j = 0;j < len;j++){
                        label[j] = tl_random_int(rng, K-2) +2;
                        if(label[j] < min){
                                min = label[j];
                        }
                        if(label[j] > max){
                                max = label[j];
                        }
                }

        }
        ASSERT(min >= 2, "Min is too small");
        ASSERT(max < K, "Max is too big");
        return OK;
ERROR:
        return FAIL;
}


int set_alphabet_and_convert(struct seq_buffer* sb, struct rng_state* rng)
{
        struct alphabet* a = NULL;
        int max = -1;
        int i;


        /* create alphabet */
        switch (sb->L) {
        case TL_SEQ_BUFFER_DNA: {
                RUN(create_alphabet(&a, rng, TLALPHABET_NOAMBIGUOUS_DNA));
                break;
        }
        case TL_SEQ_BUFFER_PROTEIN: {
                RUN(create_alphabet(&a, rng, TLALPHABET_NOAMBIGIOUS_PROTEIN));

                break;
        }
        default:
                ERROR_MSG("Alphabet type not recognized: %d", sb->L);
                break;
        }


        sb->alphabet = a;

        /* convert sequences and keep track of "biggest" letter */
        for(i = 0; i < sb->num_seq;i++){
                int len;
                uint8_t* seq;
                len = sb->sequences[i]->seq_len;
                seq = sb->sequences[i]->seq;
                RUN(convert_to_internal(a, seq, len));

        }
        //fprintf(stdout,"MAX:%d", max);
        /* Make sure the conversion worked - i.e. all letters should be in the appropriate ranges */
        switch (sb->L) {
        case TL_SEQ_BUFFER_DNA: {
                ASSERT(max < 4,"Problem with alphabet. Expecting 0-3 but got %d", max);
                break;
        }
        case TL_SEQ_BUFFER_PROTEIN: {
                ASSERT(max < 4,"Problem with alphabet. Expecting 0-19 but got %d", max);
                break;
        }
        default:
                ERROR_MSG("Alphabet type not recognized: %d (should never happen because checked above)", sb->L);
                break;

        }
        return OK;
ERROR:
        return FAIL;
}
