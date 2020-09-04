
#include "tldevel.h"
#include "tlrng.h"
#include "tlseqbuffer.h"
#include "sequence_alloc.h"
#include "sequence_struct.h"

#include "null_model_emission.h"

#define SEQUENCES_SIM_IMPORT

#include "sequences_sim.h"


int sim_sequences(int N, int L,int len,struct tl_seq_buffer** seq_buf,struct rng_state* rng)
{
        struct tl_seq_buffer* sb = NULL;
        struct tl_seq* s = NULL;
        //struct ihmm_sequence* iseq = NULL;
        double* b = NULL;
        double r;
        int i,j,c;

        if(*seq_buf){
                sb = *seq_buf;
        }else{
                RUN(alloc_tl_seq_buffer(&sb,N));
        }

        RUN(get_null_model_emissions(&b, L));

        for(i = 1;i < L;i++){
                b[i] = b[i]+b[i-1];
                //fprintf(stdout,"%d %f\n",i,b[i]);
        }

        for(i = 0; i < N;i++){
                s = sb->sequences[i];
                //iseq = sb->sequences[i];
                while(s->malloc_len <= len){
                        resize_tl_seq(s);
                }
                //if(iseq->malloc_len <= len){
                //realloc_ihmm_seq(iseq, len);
                //}
                s->len = len;
                //iseq->seq_len = len;
                for(j = 0; j < len;j++){
                        r = tl_random_double(rng);

                        for(c = 0; c < L;c++){
                                if(r < b[c]){
                                        s->seq[j] = c;
                                        break;
                                }
                        }
                        //iseq->seq[j] = tl_random_int(rng, L);
                        //fprintf(stdout,"%d,", iseq->seq[j]);
                }
                snprintf(s->name, TL_SEQ_MAX_NAME_LEN,"Seq%d",i+1);
                //fprintf(stdout,"\n");
        }
        sb->num_seq = N;
        sb->max_len = len;
        *seq_buf = sb;
        return OK;
ERROR:
        return FAIL;
}
