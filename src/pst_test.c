
#include <string.h>
#include <zlib.h>

#include "tldevel.h"
#include "tllogsum.h"
#include "tlmisc.h"
#include "tlrng.h"
#include "tlseqio.h"
#include "tlalphabet.h"

#include "esl_stopwatch.h"
#include "sim_seq_lib.h"

#include "pst_io.h"

#define PST_IMPORT
#include "pst.h"



static int test_match_insert(char* infile);

int main(int argc, char *argv[])
{
        char alphabet[] = "ACGT";
        LOG_MSG("Hello World");
        char* filename = NULL;
        char* test_seq = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct pst* p = NULL;
        struct rng_state* rng = NULL;

        int i,j,c;
        float out;
        LOG_MSG("%d",argc);
        if(argc == 2){
                filename = argv[1];
                RUN(test_match_insert(filename));
        }else{

        }
        //GGTTTACT
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int test_match_insert(char* infile)
{
        struct file_handler* f = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct pst* p = NULL;
        struct pst* p_io = NULL;
        struct rng_state* rng = NULL;

        struct alphabet* a = NULL;
        struct count_hash* h = NULL;

        double s[5];
        char* test_seq = NULL;

        float P_M;
        float P_R;
        int test_len = 0;
        int i_point;
        int i,j;
        LOG_MSG("%s",infile);
        if(!my_file_exists(infile)){
                ERROR_MSG("File %s not found");
        }
        RUN(open_fasta_fastq_file(&f, infile, TLSEQIO_READ));


        RUN(read_fasta_fastq_file(f, &sb, 1000000));


        RUN(close_seq_file(&f));


        RUNP(rng = init_rng(0));
        if(sb->L == TL_SEQ_BUFFER_DNA){
                RUN(create_alphabet(&a, rng, TLALPHABET_DEFAULT_DNA));
        }else if(sb->L == TL_SEQ_BUFFER_PROTEIN){
                RUN(create_alphabet(&a, rng, TLALPHABET_DEFAULT_PROTEIN));
        }
        for(i = 0; i < sb->num_seq;i++){
                RUN(convert_to_internal(a, (uint8_t*)sb->sequences[i]->seq, sb->sequences[i]->len));

        }

        DECLARE_TIMER(timer);
        START_TIMER(timer);
        RUN(fill_exact_hash(&h, sb));
        STOP_TIMER(timer);
        GET_TIMING(timer);

        LOG_MSG("L: %d",h->L);
        RUN(run_build_pst(&p, 0.01,h));
        free_exact_hash(h);

        RUN(write_pst_hdf5(p, "test_pst.h5"));


        RUN(read_pst_hdf5(&p_io, "test_pst.h5"));



        START_TIMER(timer);
        for(j = 0;j <  5;j++){
                s[j] = 0.0;
        }
        for (j = 0; j < sb->num_seq;j++){
                RUN(score_pst(p, sb->sequences[j]->seq, sb->sequences[j]->len, &P_M,&P_R));
                //LOG_MSG("M:%f R:%f", P_M,P_R);
                P_M = P_M - P_R;
                s[0]++;
                s[1] += P_M;
                s[2] += P_M* P_M;
        }
        s[3] = s[1] / s[0];

        s[4] = sqrt ( (s[0] * s[2] -  pow(s[1], 2.0)) /  (s[0] * ( s[0] - 1.0)));
        STOP_TIMER(timer);
        GET_TIMING(timer);

        LOG_MSG("TRAIN: %f %f  (average p = %f)", s[3],s[4],   exp2f(s[3]) / (1.0 + exp2f(s[3])));
        //LOG_MSG("TRAIN %f %f", s[3],s[4]);
        //START_TIMER(timer);

        for(j = 0;j <  5;j++){
                s[j] = 0.0;
        }

        for (j = 0; j < sb->num_seq;j++){
                s[0]++;
                s[1] += sb->sequences[j]->len;
                s[2] += sb->sequences[j]->len * sb->sequences[j]->len;
        }
        s[3] = s[1] / s[0];

        s[4] = sqrt ( (s[0] * s[2] -  pow(s[1], 2.0)) /  (s[0] * ( s[0] - 1.0)));
        test_len = s[3];

        for(j = 0;j <  5;j++){
                s[j] = 0.0;
        }
        for (j = 0; j < sb->num_seq;j++){
                RUN(generate_random_seq(&test_seq, &test_len, rng));
                RUN(convert_to_internal(a, test_seq, test_len));
                //RUN(insert_seq(test_seq, test_len, sb->sequences[0]->seq, sb->sequences[0]->len, rng));
                RUN(score_pst(p, test_seq, test_len, &P_M,&P_R));
                //LOG_MSG("M:%f R:%f", P_M,P_R);
                P_M = P_M - P_R;
                //LOG_MSG("M:%f R:%f", P_M,P_R);
                s[0]++;
                s[1] += P_M;
                s[2] += P_M* P_M;

        }
        s[3] = s[1] / s[0];

        s[4] = sqrt ( (s[0] * s[2] -  pow(s[1], 2.0)) /  (s[0] * ( s[0] - 1.0)));
        START_TIMER(timer);
        for (j = 0; j <  sb->num_seq;j++){
                score_pst(p, test_seq, test_len, &P_M,&P_R);

        }
        STOP_TIMER(timer);
        GET_TIMING(timer);
//LOG_MSG("RANDOM: %f %f", s[3],s[4]);
        LOG_MSG("RANDOM: %f %f  (average p = %f)", s[3],s[4],   exp2f(s[3]) / (1.0 + exp2f(s[3])));

        //if()
        free_pst(p);
        p = p_io;
        //DECLARE_TIMER(timer);

        START_TIMER(timer);
        for(j = 0;j <  5;j++){
                s[j] = 0.0;
        }
        for (j = 0; j < sb->num_seq;j++){
                RUN(score_pst(p, sb->sequences[j]->seq, sb->sequences[j]->len, &P_M,&P_R));
                //LOG_MSG("M:%f R:%f", P_M,P_R);
                P_M = P_M - P_R;
                s[0]++;
                s[1] += P_M;
                s[2] += P_M* P_M;
        }
        s[3] = s[1] / s[0];

        s[4] = sqrt ( (s[0] * s[2] -  pow(s[1], 2.0)) /  (s[0] * ( s[0] - 1.0)));
        STOP_TIMER(timer);
        GET_TIMING(timer);

        LOG_MSG("TRAIN: %f %f  (average p = %f)", s[3],s[4],   exp2f(s[3]) / (1.0 + exp2f(s[3])));
        //LOG_MSG("TRAIN %f %f", s[3],s[4]);
        //START_TIMER(timer);

        for(j = 0;j <  5;j++){
                s[j] = 0.0;
        }

        for (j = 0; j < sb->num_seq;j++){
                s[0]++;
                s[1] += sb->sequences[j]->len;
                s[2] += sb->sequences[j]->len * sb->sequences[j]->len;
        }
        s[3] = s[1] / s[0];

        s[4] = sqrt ( (s[0] * s[2] -  pow(s[1], 2.0)) /  (s[0] * ( s[0] - 1.0)));
        test_len = s[3];

        for(j = 0;j <  5;j++){
                s[j] = 0.0;
        }
        for (j = 0; j < sb->num_seq;j++){
                RUN(generate_random_seq(&test_seq, &test_len, rng));
                RUN(convert_to_internal(a, test_seq, test_len));
                //RUN(insert_seq(test_seq, test_len, sb->sequences[0]->seq, sb->sequences[0]->len, rng));
                RUN(score_pst(p, test_seq, test_len, &P_M,&P_R));
                //LOG_MSG("M:%f R:%f", P_M,P_R);
                P_M = P_M - P_R;
                //LOG_MSG("M:%f R:%f", P_M,P_R);
                s[0]++;
                s[1] += P_M;
                s[2] += P_M* P_M;

        }
        s[3] = s[1] / s[0];

        s[4] = sqrt ( (s[0] * s[2] -  pow(s[1], 2.0)) /  (s[0] * ( s[0] - 1.0)));
        START_TIMER(timer);
        for (j = 0; j <  sb->num_seq;j++){
                score_pst(p, test_seq, test_len, &P_M,&P_R);

        }
        STOP_TIMER(timer);
        GET_TIMING(timer);
//LOG_MSG("RANDOM: %f %f", s[3],s[4]);
        LOG_MSG("RANDOM: %f %f  (average p = %f)", s[3],s[4],   exp2f(s[3]) / (1.0 + exp2f(s[3])));


        DESTROY_TIMER(timer);
        if(test_seq){
                MFREE(test_seq);
        }
        free_alphabet(a);

        free_pst(p);
        free_tl_seq_buffer(sb);
        free_rng(rng);
        return OK;
ERROR:
        return FAIL;
}
