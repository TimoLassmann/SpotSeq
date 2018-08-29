#include "emit_random.h"
#include "ihmm_seq.h"
#include "distributions.h"

int emit_sequences_from_random_model(char* filename, int num)
{
        struct seq_buffer* sb_in  = NULL;
        struct seq_buffer* sb_out = NULL;

        int i,j;
        double s1,s2;
        double s1_t,s2_t;

        double r;

        rk_state rndstate;

        RUN(rk_randomseed(&rndstate));

        RUNP(sb_in = get_sequences_from_hdf5_model(filename));

        ASSERT(sb_in != NULL, "No sequence Buffer");
        s1 = 0.0;
        s2 = 0.0;
        for(i = 0; i < sb_in->num_seq;i++){
                //sb_in->sequences[i]->seq_len = 10 + (int)(rk_double(&rndstate)*10.0) - 5.0;
                s1 += sb_in->sequences[i]->seq_len;
                s2 += (sb_in->sequences[i]->seq_len * sb_in->sequences[i]->seq_len);
        }

        s2 = sqrt(((double) sb_in->num_seq * s2 - s1 * s1)/ ((double) sb_in->num_seq * ((double) sb_in->num_seq -1.0)));
        s1 = s1 / (double) sb_in->num_seq;



        LOG_MSG("Mean length: %f stdev: %f",s1,s2);
        for(j = 0; j < 10;j++){
                s1_t = 0.0;
                s2_t = 0.0;
                for(i = 0;i < sb_in->num_seq;i++){
                        r = rk_normal(&rndstate, s1,s2);
                        s1_t += r;
                        s2_t += r*r;
                        //fprintf(stdout,"%f ", r);
                }
                s2_t = sqrt(((double) sb_in->num_seq * s2_t - s1_t * s1_t)/ ((double) sb_in->num_seq * ((double) sb_in->num_seq -1.0)));
                s1_t = s1_t / (double) sb_in->num_seq;

                LOG_MSG("Mean length: %f stdev: %f",s1_t,s2_t);
        }

        MMALLOC(sb_out,sizeof(struct seq_buffer));
        sb_out->malloc_num = num ;
        sb_out->num_seq = num;
        sb_out->sequences = NULL;
        sb_out->max_len = 0.0;
        sb_out->L = 12;

        MMALLOC(sb_out->sequences, sizeof(struct chromosome*) *sb_out->malloc_num );
        for(i = 0; i < sb_out->num_seq;i++){
                sb_out->sequences[i] = NULL;
                RUNP(sb_out->sequences[i] = alloc_ihmm_seq());
                while(sb_out->sequences[i]->malloc_len <= sb_out->max_len){
                        RUN(realloc_ihmm_seq(sb_out->sequences[i]));
                }
        }
        free_ihmm_sequences(sb_in);
        free_ihmm_sequences(sb_out);
        return OK;
ERROR:
        return FAIL;
}

int main(const int argc, char * argv[])
{
        RUN(emit_sequences_from_random_model(argv[1], 100));
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}
