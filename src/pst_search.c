#include <omp.h>
#include <string.h>

#include "tldevel.h"
#include "tlrng.h"
#include "tlseqio.h"
#include "tlmisc.h"
#include "tlalphabet.h"

#include "pst.h"
#include "pst_structs.h"
#include "pst_hash.h"

#include "search_db.h"

#define  PST_SEARCH_IMPORT
#include "pst_search.h"

static int copy_sequences(struct tl_seq_buffer* sb, struct tl_seq* a);

int search_db(struct pst* p, char* filename, double thres,struct tl_seq_buffer** hits, uint64_t* db_size)
{

#ifdef HAVE_OPENMP
        omp_lock_t writelock;
        omp_set_num_threads(8);
        omp_init_lock(&writelock);
#endif
        struct file_handler* f = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct rng_state* rng = NULL;
        struct alphabet* alphabet = NULL;
        struct tl_seq_buffer* h = NULL;
        int chunk,i;
        int n_hits = 0;
        uint64_t total_nseq =0;

        if(*hits){
                h = *hits;
        }else{
                RUN(alloc_tl_seq_buffer(&h, 10000));
        }

        //LOG_MSG("%s",infile);
        if(!my_file_exists(filename)){
                ERROR_MSG("File %s not found");
        }

        RUNP(rng = init_rng(42));
        RUN(open_fasta_fastq_file(&f, filename, TLSEQIO_READ));
        chunk =1;
        total_nseq = 0;
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
                total_nseq+= sb->num_seq;
                if(sb->num_seq == 0){
                        break;
                }
                for(i = 0; i < sb->num_seq;i++){
                        struct tl_seq* seq = sb->sequences[i];
                        int len = seq->len;
                        convert_to_internal(alphabet, (uint8_t*)seq->seq,len);
                }
#ifdef HAVE_OPENMP


#pragma omp parallel default(shared)
#pragma omp for schedule(dynamic) nowait
#endif
                        for(i = 0; i < sb->num_seq;i++){
                                double z_score;
                                float score;

                                //len = MACRO_MIN(len, 100);
                                score_pst(p, sb->sequences[i]->seq, sb->sequences[i]->len, &score);

                                z_score_pst(p, sb->sequences[i]->len, score, &z_score);
                                if(z_score >= thres){
                                        omp_set_lock(&writelock);
                                        //int thread_id  = omp_get_thread_num();
                                        //LOG_MSG("Thread %d locking ", thread_id);
                                        copy_sequences(h, sb->sequences[i]);
                                        n_hits++;
                                        if(sb->sequences[i]->len > h->max_len){
                                                h->max_len = sb->sequences[i]->len;
                                        }
                                        omp_unset_lock(&writelock);

                                        //fprintf(stdout,"Hit: %f\t%s\n",z_score,seq->name);
                                }
                        }
//#ifdef HAVE_OPENMP
                        //              }
//#endif
                //free_tl_seq_buffer(sb);
                //sb = NULL;
                        //break;
                chunk++;
        }
        RUN(close_seq_file(&f));
        free_rng(rng);
        free_tl_seq_buffer(sb);
        free_alphabet(alphabet);
        LOG_MSG("Scanned %0.2f M sequences.", (double)total_nseq / 1000000.0);
        //LOG_MSG("Found %d hits", n_hits);
        *hits = h;
        *db_size = total_nseq;
#ifdef HAVE_OPENMP
        omp_destroy_lock(&writelock);
#endif
        return OK;
ERROR:
        return FAIL;
}

int copy_sequences(struct tl_seq_buffer* sb, struct tl_seq* a)
{
        int printed;
        struct tl_seq* target = NULL;
        if(sb->num_seq == sb->malloc_num){
                RUN(resize_tl_seq_buffer(sb));
        }

        target = sb->sequences[sb->num_seq];

        printed = snprintf(target->name, TL_SEQ_MAX_NAME_LEN, "%s", a->name);
        ASSERT(printed < TL_SEQ_MAX_NAME_LEN ,"characters printed entirely fills buffer");

        while(target->malloc_len <= a->len){
                RUN(resize_tl_seq(target));
        }
        memcpy(target->seq, a->seq, a->len);
        target->len = a->len;

        sb->num_seq++;
        return OK;
ERROR:
        return FAIL;
}

int search_db_hdf5(struct pst* p, char* filename, double thres)
{
        int chunk,i;
#ifdef HAVE_OPENMP
        omp_lock_t writelock;
        omp_set_num_threads(8);


        omp_init_lock(&writelock);
#endif


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

#pragma omp parallel shared(h) private(i)
                {
#pragma omp for schedule(dynamic) nowait
#endif
                        for(i = 0; i < h->num_seq;i++){



                                uint8_t* seq = h->seq + h->len[i];
                                int len = h->len[i+1] - h->len[i];

                                double z_score;
                                float score;
                                score_pst(p, seq, len, &score);
                                z_score_pst(p, len, score, &z_score);
                                if(z_score >= thres){
                                        omp_set_lock(&writelock);
                                        hits++;
                                        omp_unset_lock(&writelock);
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

#ifdef HAVE_OPENMP
        omp_destroy_lock(&writelock);
#endif

        return OK;
ERROR:
        return FAIL;
}
