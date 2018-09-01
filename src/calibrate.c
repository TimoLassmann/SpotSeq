#include "calibrate.h"

#include "ihmm_seq.h"
#include "emit_random.h"
#include "run_score.h"
#include "finite_hmm.h"


int calibrate(char* model_file,int num_threads, float* lambda, float* beta)
{
        struct seq_buffer* sb = NULL;
        struct seq_buffer* sb_in = NULL;
        struct fhmm* fhmm = NULL;

        int num_seq = 10;
        float* scores = NULL;
        int i;



        /* Step 1: emit sequences with same background and length
         * distribution as training sequences */
        LOG_MSG("Generating random sequences.");
        RUNP(sb_in = get_sequences_from_hdf5_model(model_file));

        RUNP(sb = emit_sequences_from_random_model(sb_in, num_seq));
        free_ihmm_sequences(sb_in);
        sb_in = NULL;

        //sb = emit_sequences_from_random_model(struct seq_buffer *sb_in, int num)

        /* Step 2: score random sequence, get scores into an array */

        LOG_MSG("Loading model.");
        RUNP(fhmm = init_fhmm(model_file));


        RUN(run_score_sequences(fhmm, sb,num_threads));

        /* copy scores over */
        MMALLOC(scores, sizeof(float)* num_seq);
        for(i = 0; i < sb->num_seq;i++){
                scores[i] = sb->sequences[i]->score;
        }

        /* Step 3: fit gumbel distribution  */


        free_ihmm_sequences(sb_in);
        free_ihmm_sequences(sb);
        MFREE(scores);

        return OK;
ERROR:
        free_ihmm_sequences(sb);
        free_ihmm_sequences(sb_in);
        MFREE(scores);
        return FAIL;
}

//#ifdef ITEST



int main (int argc,char * argv[])
{
        float lambda, beta;
        fprintf(stdout,"Running libhdf5glue sanity tests\n");
        RUN(calibrate("Standard_Challenge_GCCTAGGAGTCGGTT_mis_0_10.fa.h5", 8, &lambda, &beta));
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


//#endif
