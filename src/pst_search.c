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

        double m,z_score;
        float P_M, P_R;
        int chunk,i,len;
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

                        m = a + b * (double) len;
                        z_score = (P_M - m) / v;
                        if(z_score >= thres){
                                fprintf(stdout,"Hit: %f\t%s\n",z_score,sb->sequences[i]->name);
                        }


                }
                chunk++;
        }

        RUN(close_seq_file(&f));

        free_rng(rng);
        free_tl_seq_buffer(sb);
        free_alphabet(alphabet);
        return OK;
ERROR:
        return FAIL;
}
