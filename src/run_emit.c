#include "emit_random.h"

struct seq_buffer* emit_sequences_from_random_model(struct seq_buffer* sb_in, int num)
{
        struct seq_buffer* sb_out = NULL;

        int i,j,c;
        double s1,s2;
        //double s1_t,s2_t;
        double r;

        rk_state rndstate;

        RUN(rk_randomseed(&rndstate));

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

        /*for(j = 0; j < 10;j++){
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
        */

        MMALLOC(sb_out,sizeof(struct seq_buffer));
        sb_out->malloc_num = num ;
        sb_out->num_seq = num;
        sb_out->sequences = NULL;
        sb_out->max_len = 0.0;
        sb_out->L = sb_in->L;

        sb_out->background = NULL;
        MMALLOC(sb_out->background,sizeof(float) * sb_out->L );
        for(i = 0;i < sb_out->L;i++){
                sb_out->background[i] = sb_in->background[i];
                //fprintf(stdout,"%d %f",i,sb_out->background[i]);
        }


        for(i = 1;i < sb_out->L;i++){
                sb_out->background[i] += sb_out->background[i-1];
        }


        MMALLOC(sb_out->sequences, sizeof(struct chromosome*) *sb_out->malloc_num );
        for(i = 0; i < sb_out->num_seq;i++){
                sb_out->sequences[i] = NULL;
                RUNP(sb_out->sequences[i] = alloc_ihmm_seq());


                snprintf(sb_out->sequences[i]->name, 256, "RANDOM%d", i+1);
                sb_out->sequences[i]->seq_len = rk_normal(&rndstate, s1,s2);
                if( sb_out->sequences[i]->seq_len > sb_out->max_len){
                        sb_out->max_len = sb_out->sequences[i]->seq_len;
                }
                while(sb_out->sequences[i]->malloc_len <= sb_out->sequences[i]->seq_len){
                        RUN(realloc_ihmm_seq(sb_out->sequences[i]));
                }

                for(j = 0;j < sb_out->sequences[i]->seq_len;j++){
                        //r = random_float_zero_to_x(1.0);
                         r = rk_double(&rndstate);
                         //fprintf(stdout,"%f\n",r);
                         for(c = 0 ; c < sb_out->L;c++){
                                 if(r <= sb_out->background[c]){
                                         sb_out->sequences[i]->seq[j] = c;
                                         break;
                                 }
                         }
                }
        }

        return sb_out;
ERROR:
        return NULL;
}
