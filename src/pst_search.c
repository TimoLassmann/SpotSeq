#include "tldevel.h"
#include "tlrng.h"
#include "tlseqio.h"
#include "tlmisc.h"
#include "tlalphabet.h"

#include "pst.h"
#include "pst_structs.h"

#define  PST_SEARCH_IMPORT
#include "pst_search.h"


int search_db(struct pst* p, char* filename, double thres)
{
        struct file_handler* f = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct rng_state* rng = NULL;
        struct alphabet* alphabet = NULL;

        double z_score;
        float P_M, P_R;
        int chunk,i,len,index;
        double a,b,v;

        a = p->a;
        b = p->b;
        v = p->var;

        //LOG_MSG("%s",infile);
        if(!my_file_exists(filename)){
                ERROR_MSG("File %s not found");
        }
        RUNP(rng = init_rng(0));

        RUN(open_fasta_fastq_file(&f, filename, TLSEQIO_READ));

        FILE* f_ptr = NULL;

        f_ptr = fopen("scores.txt", "w");
        chunk =1;
        while(1){
                RUN(read_fasta_fastq_file(f, &sb, 100000));
                if(chunk == 1){
                        if(sb->L == TL_SEQ_BUFFER_DNA){
                                RUN(create_alphabet(&alphabet, rng, TLALPHABET_NOAMBIGUOUS_DNA));
                        }else if(sb->L == TL_SEQ_BUFFER_PROTEIN){
                                RUN(create_alphabet(&alphabet, rng, TLALPHABET_NOAMBIGIOUS_PROTEIN ));
                        }
                }
                if(sb->num_seq == 0){
                        break;
                }
                LOG_MSG("Working on chunk: %d",chunk);
                for(i = 0; i < sb->num_seq;i++){
                        len = sb->sequences[i]->len;
                        RUN(convert_to_internal(alphabet, (uint8_t*)sb->sequences[i]->seq,len));
                        RUN(score_pst(p, sb->sequences[i]->seq, len, &P_M,&P_R));
                        P_M = P_M - P_R;
                        index = MACRO_MIN(p->max_observed_len, len);
                        index = p->fit_index[index];
                        a = p->fit[index][PST_FIT_A];
                        b = p->fit[index][PST_FIT_B];
                        v = p->fit[index][PST_FIT_V];

                        z_score = (P_M - (a + b * (double)len )) / v;
                        if(z_score >= thres){
                                fprintf(stdout,"Hit: %f\t%s\n",z_score,sb->sequences[i]->name);
                        }
                        fprintf(f_ptr,"%d,%f\n",len,z_score);

                }
                chunk++;
        }

        RUN(close_seq_file(&f));
        fclose(f_ptr);
        free_rng(rng);
        free_tl_seq_buffer(sb);
        free_alphabet(alphabet);
        return OK;
ERROR:
        return FAIL;
}
