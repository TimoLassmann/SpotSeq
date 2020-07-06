#include "tldevel.h"
#include "tllogsum.h"
#include "tlmisc.h"
#include "tlrng.h"
#include "tlseqio.h"
#include "tlalphabet.h"

#include "esl_stopwatch.h"

#define  PST_BUILD
#include "pst_build.h"

#include "pst_hash.h"
#include "pst_calibrate.h"
#include "pst_io.h"
#include "pst.h"

static int recursive_rename(char* name,  int i);

static int read_training_sequences(struct tl_seq_buffer** seq_buf,struct rng_state* rng,char* infile);

int create_pst_model(struct rng_state* rng, char* train_seq,char* seq_db,char* out_model, double p_min, double gamma, double z_thres)
{

        struct tl_seq_buffer* sb = NULL;
        struct pst* p = NULL;
        struct count_hash* h = NULL;
        /* sanity checks */

        /* 1: files  */

        if(!my_file_exists(train_seq)){
                ERROR_MSG("File %s not found.", train_seq);
        }
        if(!my_file_exists(seq_db)){
                ERROR_MSG("File %s not found.", seq_db);
        }
        /* if the outmodel exits rename to <out_model>_(i+1) and continue */
        RUN(recursive_rename(out_model, 0));

        /* set default  parameters  */
        if(TLSAFE_EQ(p_min, 0.0)){
                p_min = 0.0001;
        }
        if(TLSAFE_EQ(gamma, 0.0)){
                gamma = 0.01;
        }
        if(TLSAFE_EQ(z_thres, 0.0)){
                z_thres = 10.0;
        }
        LOG_MSG("Running with: %f %f %f", p_min, gamma, z_thres);

        RUN(read_training_sequences(&sb, rng, train_seq));
        DECLARE_TIMER(timer);
        START_TIMER(timer);
        RUN(fill_exact_hash(&h, sb));
        STOP_TIMER(timer);
        GET_TIMING(timer);

        free_tl_seq_buffer(sb);

        RUN(run_build_pst(&p,(float) p_min,(float)gamma,h));
        free_exact_hash(h);
        LOG_MSG("Lets calibrate");
        START_TIMER(timer);
        RUN(calibrate_pst(p, seq_db,z_thres));
        STOP_TIMER(timer);
        GET_TIMING(timer);
        DESTROY_TIMER(timer);


        RUN(write_pst_hdf5(p,   out_model));

        free_pst(p);
        return OK;
ERROR:
        return FAIL;
}

int read_training_sequences(struct tl_seq_buffer** seq_buf,struct rng_state* rng,char* infile)
{
        struct file_handler* f = NULL;
        struct tl_seq_buffer* sb = NULL;

        struct alphabet* a = NULL;
        int i;

        RUN(open_fasta_fastq_file(&f, infile, TLSEQIO_READ));

        RUN(read_fasta_fastq_file(f, &sb, 100000));
        if(sb->num_seq == 100000){
                WARNING_MSG("More than 100K training sequences.");
        }
        RUN(close_seq_file(&f));

        if(sb->L == TL_SEQ_BUFFER_DNA){
                RUN(create_alphabet(&a, rng, TLALPHABET_NOAMBIGUOUS_DNA));
        }else if(sb->L == TL_SEQ_BUFFER_PROTEIN){
                RUN(create_alphabet(&a, rng, TLALPHABET_NOAMBIGIOUS_PROTEIN));
        }

        for(i = 0; i < sb->num_seq;i++){
                RUN(convert_to_internal(a, (uint8_t*)sb->sequences[i]->seq, sb->sequences[i]->len));
        }
        free_alphabet(a);
        *seq_buf = sb;
        return OK;
ERROR:
        return NULL;
}

int recursive_rename(char* name,  int i)
{
        char* buffer1 = NULL;
        int b_len1 = 256;
        char* buffer2 = NULL;
        int b_len2 = 256;
        int rc;


        MMALLOC(buffer1, sizeof(char)* b_len1);
        MMALLOC(buffer2, sizeof(char)* b_len2);
        if(i){
                rc = snprintf(buffer1, b_len1, "%s_%d",name,i);
                while(rc >= b_len1){
                        b_len1 = b_len1 + b_len1 /2;
                        MREALLOC(buffer1, sizeof(char) * b_len1);
                        rc = snprintf(buffer1, b_len1, "%s_%d",name,i);

                }
        }else{
                rc = snprintf(buffer1, b_len1, "%s",name);
                while(rc >= b_len1){
                        b_len1 = b_len1 + b_len1 /2;
                        MREALLOC(buffer1, sizeof(char) * b_len1);
                        rc = snprintf(buffer1, b_len1, "%s",name);

                }
        }
        rc = snprintf(buffer2, b_len2, "%s_%d",name,i+1);
        while(rc >= b_len2){
                b_len2 = b_len2 + b_len2 /2;
                MREALLOC(buffer2, sizeof(char) * b_len2);
                rc = snprintf(buffer2, b_len2, "%s_%d",name,i+1);

        }

        //LOG_MSG("Searching for: %s %s", buffer1,buffer2);
        while (1){
                if(my_file_exists(buffer1) && my_file_exists(buffer2)){
                        /* both files exits - look deeper */
                        //LOG_MSG("both  %s %s", buffer1,buffer2);
                        RUN(recursive_rename(name, i+1));
                }else if(my_file_exists(buffer1) && ! my_file_exists(buffer2)){
                        /* buffer 2 is name of new file */
                        //LOG_MSG(" mv %s %s", buffer1,buffer2);
                        RUN(rename(buffer1, buffer2));
                }else if(!my_file_exists(buffer1)){
                        /*  I am done.  */
                        //LOG_MSG("done  %s %s", buffer1,buffer2);
                        break;
                }
        }

        MFREE(buffer1);
        MFREE(buffer2);
        return OK;
ERROR:
        return FAIL;

}
