
#include "tldevel.h"
#include "tlrng.h"
#include "sequence_alloc.h"
#include "sequence_struct.h"

#define SEQUENCES_SIM_IMPORT
#include "sequences_sim.h"


int sim_sequences(int N, int L,int len,struct seq_buffer** seq_buf,struct rng_state* rng)
{
        struct seq_buffer* sb = NULL;
        struct ihmm_sequence* iseq = NULL;
        double* b = NULL;
        double r;
        int i,j,c;

        if(*seq_buf){
                sb = *seq_buf;
        }else{
                RUN(alloc_sequence_buffer(&sb, N));
        }

        RUN(get_null_model_emissions(&b, L));

        for(i = 1;i < L;i++){
                b[i] = b[i]+b[i-1];
                //fprintf(stdout,"%d %f\n",i,b[i]);
        }

        for(i = 0; i < N;i++){
                iseq = sb->sequences[i];
                if(iseq->malloc_len <= len){
                        realloc_ihmm_seq(iseq, len);
                }
                iseq->seq_len = len;
                for(j = 0; j < len;j++){
                        r = tl_random_double(rng);

                        for(c = 0; c < L;c++){
                                if(r < b[c]){
                                        iseq->seq[j] = c;
                                        break;
                                }
                        }
                        //iseq->seq[j] = tl_random_int(rng, L);
                        //fprintf(stdout,"%d,", iseq->seq[j]);
                }
                snprintf(iseq->name, MAX_SEQUENCE_NAME_LEN,"Seq%d",i+1);
                //fprintf(stdout,"\n");
        }
        sb->num_seq = N;
        sb->max_len = len;
        *seq_buf = sb;
        return OK;
ERROR:
        return FAIL;
}
