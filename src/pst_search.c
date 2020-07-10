#include <omp.h>
#include "tldevel.h"
#include "tlrng.h"
#include "tlseqio.h"
#include "tlmisc.h"
#include "tlalphabet.h"

#include "pst.h"
#include "pst_structs.h"

#include "search_db.h"

#define  PST_SEARCH_IMPORT
#include "pst_search.h"

int search_db(struct pst* p, char* filename, double thres)
{
        struct file_handler* f = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct rng_state* rng = NULL;
        struct alphabet* alphabet = NULL;

        int chunk,i;
        int hits = 0;
        //LOG_MSG("%s",infile);
        if(!my_file_exists(filename)){
                ERROR_MSG("File %s not found");

        }

        RUNP(rng = init_rng(0));
        RUN(open_fasta_fastq_file(&f, filename, TLSEQIO_READ));
        chunk =1;
        while(1){
                RUN(read_fasta_fastq_file(f, &sb, 1000000));
                int alloc = 0;

                for(i = 0; i < sb->num_seq;i++){
                        alloc+= sb->sequences[i]->malloc_len;
                }
                LOG_MSG("CHUNK:%d %d",chunk, alloc);
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
#ifdef HAVE_OPENMP
                omp_set_num_threads(4);
#pragma omp parallel shared(sb,alphabet,p) private(i)
                {
#pragma omp for schedule(dynamic) nowait
#endif
                        for(i = 0; i < sb->num_seq;i++){
                                struct tl_seq* seq = sb->sequences[i];
                                double a,b,v;
                                int len = seq->len;
                                double z_score;
                                float P_M, P_R;

                                convert_to_internal(alphabet, (uint8_t*)seq->seq,len);
                                //len = MACRO_MIN(len, 100);
                                score_pst(p, seq->seq, len, &P_M,&P_R);

                                P_M = P_M - P_R;
                                int index = MACRO_MIN(p->max_observed_len, len);
                                index = p->fit_index[index];
                                a = p->fit[index][PST_FIT_A];
                                b = p->fit[index][PST_FIT_B];
                                v = p->fit[index][PST_FIT_V];
                                z_score = (P_M - (a + b * (double)len )) / v;
                                if(z_score >= thres){
                                        hits++;
                                        //fprintf(stdout,"Hit: %f\t%s\n",z_score,seq->name);
                                }
                        }
#ifdef HAVE_OPENMP
                }
#endif
                //free_tl_seq_buffer(sb);
                //sb = NULL;
                chunk++;
        }
        RUN(close_seq_file(&f));
        free_rng(rng);
        free_tl_seq_buffer(sb);
        free_alphabet(alphabet);
        LOG_MSG("Found %d hits", hits);
        return OK;
ERROR:
        return FAIL;
}



int search_db_hdf5(struct pst* p, char* filename, double thres)
{
        int chunk,i;

        //LOG_MSG("%s",infile);
        if(!my_file_exists(filename)){
                ERROR_MSG("File %s not found");

        }


        struct hdf_seq_store*h = NULL;
        int numseq;
        int hits = 0;
        while(1){
                RUN(read_hdf_seq_store_chunk(&h, filename));
                if(!h->num_seq){
                        break;
                }
#ifdef HAVE_OPENMP
                omp_set_num_threads(8);
#pragma omp parallel shared(h) private(i)
                {
#pragma omp for schedule(dynamic) nowait
#endif
                        for(i = 0; i < h->num_seq;i++){



                                uint8_t* seq = h->seq + h->len[i];
                                int len = h->len[i+1] - h->len[i];

                                double a,b,v;

                                double z_score;
                                float P_M, P_R;
                                score_pst(p, seq, len, &P_M,&P_R);

                                P_M = P_M - P_R;
                                int index = MACRO_MIN(p->max_observed_len, len);
                                index = p->fit_index[index];
                                a = p->fit[index][PST_FIT_A];
                                b = p->fit[index][PST_FIT_B];
                                v = p->fit[index][PST_FIT_V];
                                z_score = (P_M - (a + b * (double)len )) / v;
                                if(z_score >= thres){
                                        hits++;
                                        //fprintf(stdout,"Hit: %f\t%s\n",z_score,seq->name);
                                }
                        }
#ifdef HAVE_OPENMP
                }

#endif


                numseq+= h->num_seq;
                LOG_MSG("Scanned %d seqs",numseq);
        }
        free_hdf_seq_store(h);
        LOG_MSG("Found %d hits", hits);
        return OK;
ERROR:
        return FAIL;
}
