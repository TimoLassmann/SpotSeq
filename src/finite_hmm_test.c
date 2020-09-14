#include "tldevel.h"
#include "tllogsum.h"
#include "tlrng.h"

#include "finite_hmm_struct.h"
#include "finite_hmm.h"
#include "finite_hmm_alloc.h"
#include "finite_hmm_stats.h"
#include "finite_hmm_plot.h"
#include "finite_hmm_score.h"

#include "global.h"

#include "null_model_emission.h"

static int generate_simple_fhmm(struct fhmm** f);

static int run_forward_diff_len(struct fhmm* fhmm,struct fhmm_dyn_mat*dm, uint8_t* seq, int len);

static int random_seq_test(struct fhmm* fhmm, struct fhmm_dyn_mat*dm);
/* Purpose: test fhmm search scoring */
int main(void)
{
        struct fhmm* fhmm = NULL;
        struct fhmm_dyn_mat* dm = NULL;
        uint8_t test_seq[100];

        init_logsum();

        RUN(generate_simple_fhmm(&fhmm));
        RUN(alloc_fhmm_dyn_mat(&dm, 1024, fhmm->K));

        //configure_target_len(fhmm, 10, 0);
        //fhmm_calibrate(fhmm, dm, 42);
        //exit(0);
        RUN(plot_finite_hmm_dot(fhmm, "fhmm_test_model.dot",0.01f));

        test_seq[0] = 0;
        test_seq[1] = 1;
        test_seq[2] = 2;
        test_seq[3] = 3;
        LOG_MSG("Single  hit test");
        RUN(run_forward_diff_len(fhmm,dm,  test_seq, 4));


        //exit(0);
        //exit(0);
        test_seq[0] = 0;
        test_seq[1] = 3;
        test_seq[2] = 1;
        test_seq[3] = 2;
        test_seq[4] = 3;
        LOG_MSG("Single  hit test");
        RUN(run_forward_diff_len(fhmm,dm,  test_seq, 5));
        //exit(0);
        test_seq[0] = 3;
        test_seq[1] = 1;
        test_seq[2] = 2;
        test_seq[3] = 3;
        test_seq[4] = 0;
        LOG_MSG("Single  hit test");
        RUN(run_forward_diff_len(fhmm,dm, test_seq, 5));

        test_seq[0] = 0;
        test_seq[1] = 0;
        test_seq[2] = 3;
        test_seq[3] = 1;
        test_seq[4] = 2;
        test_seq[5] = 3;
        test_seq[6] = 0;
        test_seq[7] = 0;
        LOG_MSG("Single  hit test");
        RUN(run_forward_diff_len(fhmm,dm,  test_seq, 8));
        exit(0);
        test_seq[0] = 0;
        test_seq[1] = 0;
        test_seq[2] = 0;
        test_seq[3] = 0;
        test_seq[4] = 3;
        test_seq[5] = 1;
        test_seq[6] = 2;
        test_seq[7] = 3;
        test_seq[8] = 0;
        test_seq[9] = 0;
        test_seq[10] = 0;
        test_seq[11] = 0;

        LOG_MSG("Single  hit test");
        RUN(run_forward_diff_len(fhmm, dm, test_seq, 12));

        test_seq[0] = 0;
        test_seq[1] = 3;
        test_seq[2] = 1;
        test_seq[3] = 2;
        test_seq[4] = 3;
        test_seq[5] = 0;
        test_seq[6] = 0;
        test_seq[7] = 0;
        test_seq[8] = 1;
        test_seq[9] = 2;
        test_seq[10] = 3;
        test_seq[11] = 0;

        LOG_MSG("Multi  hit test");
        RUN(run_forward_diff_len(fhmm,dm,  test_seq, 12));


        //random_seq_test(fhmm,dm);
        free_fhmm_dyn_mat(dm);
        free_fhmm(fhmm);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int random_seq_test(struct fhmm* fhmm,struct fhmm_dyn_mat* dm)
{
        struct rng_state* rng = NULL;
        uint8_t* seq = NULL;
        int malloc_len;
        int len;
        int i,j,c;
        double avg;

        RUNP(rng = init_rng(42));

        malloc_len = 10000;

        //realloc_dyn_matrices(fhmm, malloc_len + 2);
        FILE* f_ptr = NULL;


        RUNP(f_ptr = fopen("scores.csv", "w"));
        MMALLOC(seq, sizeof(uint8_t) * malloc_len);
        for(i = 400;i < 500;i+=100){
                len = i;
                avg = 0.0;
                if(len+2 >= dm->alloc_matrix_len){
                        resize_fhmm_dyn_mat(dm, len, fhmm->K);
                }
                for(c = 0; c < 1000000;c++){
                        for(j = 0; j < len;j++){

                                seq[j] = tl_random_int(rng, 4);
                                //fprintf(stdout,"%d",seq[j]);
                        }

                        //configure_target_len(fhmm, len, 1);
                        forward(fhmm, dm, &fhmm->f_score, seq, len,1);
                        random_model_score(len,&fhmm->r_score);// ,seq, len,len );
                        fprintf(f_ptr,"%f\n", (fhmm->f_score - fhmm->r_score) / logf(2.0));
                        avg += LOGISTIC_FLT(fhmm->f_score - fhmm->r_score);//mm->f_score - fhmm->r_score;//(fhmm->f_score- fhmm->r_score) / log(2.0);
                        //fprintf(stdout,"SCORE (uni) %f  %f  %f\t", fhmm->f_score, fhmm->r_score, LOGISTIC_FLT(fhmm->f_score - fhmm->r_score));
                }
                fprintf(stdout," %d %f\n",len, avg / 1000.0);
        }
        MFREE(seq);
        fclose(f_ptr);
        return OK;
ERROR:
        if(seq){
                MFREE(seq);
        }
        return FAIL;
}

int run_forward_diff_len(struct fhmm* fhmm, struct fhmm_dyn_mat*dm, uint8_t* seq, int len)
{
        int* path = NULL;
        int i;
        LOG_MSG("SeqLen: %d", len);
        for(i = MACRO_MAX(1, len -10); i < len+10   ;i++){

                //configure_target_len(fhmm, MACRO_MAX(1, i-2), 0);
                //configure_target_len(fhmm, i, 0);
                forward(fhmm, dm, &fhmm->f_score, seq, len,0);
                backward(fhmm, dm, &fhmm->b_score, seq, len,0);
                random_model_score(len,&fhmm->r_score);// ,seq, len,len );


                //random_model_score(fhmm->background , &fhmm->r_score ,seq, len,len);
                fprintf(stdout,"LEN: %d\t",i);
                fprintf(stdout,"SCORE (  uni) %f\n", fhmm->f_score);
                fprintf(stdout,"LEN: %d\t",i);
                fprintf(stdout,"SCORE (  uni) %f %f %f  %f  bit:%f\n", fhmm->b_score,scaledprob2prob(fhmm->f_score) - scaledprob2prob(fhmm->b_score),  fhmm->r_score, LOGISTIC_FLT(fhmm->f_score - fhmm->r_score),(fhmm->f_score - fhmm->r_score) / logf(2.0f));


                //posterior_decoding(fhmm, dm, fhmm->f_score, seq, len, path);

                //exit(0);
                //break;

                //configure_target_len(fhmm, i, 1);
                //configure_target_len(fhmm, MACRO_MAX(1, i-2), 1);
                forward(fhmm, dm, &fhmm->f_score, seq, len,1);
                backward(fhmm, dm, &fhmm->b_score, seq, len,1);
                random_model_score(len,&fhmm->r_score);// ,seq, len,len );
                fprintf(stdout,"LEN: %d\t",i);
                fprintf(stdout,"SCORE (multi) %f \n", fhmm->f_score);
                fprintf(stdout,"LEN: %d\t",i);
                fprintf(stdout,"SCORE (multi) %f %f %f  %f\n", fhmm->b_score, scaledprob2prob(fhmm->f_score) - scaledprob2prob(fhmm->b_score),fhmm->r_score, LOGISTIC_FLT(fhmm->f_score - fhmm->r_score));
                posterior_decoding(fhmm, dm, fhmm->f_score, seq, len, path);
                //posterior_decoding(fhmm, dm, fhmm->f_score, seq, len, path);
                //LOG_MSG("%f", scaledprob2prob(fhmm->f_score - fhmm->b_score));
                //exit(0);
        }

        return OK;
}

/* generate simple HMM A->C->G->T  */


int generate_simple_fhmm(struct fhmm** f)
{
        struct fhmm* fhmm = NULL;
        double* back = NULL;
        int i,j,o;
        RUNP(fhmm = alloc_fhmm());

        //fhmm->alloc_K = 1 + 5;
        fhmm->K = 4;//model->num_states;
        fhmm->L = 4;
        fhmm->alloc_K = 4;



        RUN(get_null_model_emissions(&back, fhmm->L));

        RUN(galloc(&fhmm->background,fhmm->L));
        for(i = 0;i < fhmm->L;i++){
                fhmm->background[i] = (float) back[i];
        }

        gfree(back);
        //RUN(alloc_dyn_matrices(fhmm));

        RUN(galloc(&fhmm->e, fhmm->K, fhmm->L));
        RUN(galloc(&fhmm->t, fhmm->K, fhmm->K));
        //RUN(galloc(&s2_t, fhmm->alloc_K, fhmm->alloc_K));
        for(i = 0; i < fhmm->K;i++){
                for(j = 0;j < fhmm->L;j++){
                        fhmm->e[i][j] = 0.0;
                }
                for(j = 0;j < fhmm->K;j++){
                        fhmm->t[i][j] = 0.0;
                }
        }

        /* set transition and emission in main model */

        for(i = 1;i < fhmm->K;i++){
                fhmm->t[i-1][i] = 1.0;
        }

        for(i = 0; i < fhmm->K;i++){
                for(j = 0; j < 4;j++){
                        fhmm->e[i][j] = (1.0 - 0.99) / 3.0;
                        if(i%4 == j){
                                fhmm->e[i][j] = 0.99;
                        }
                }
        }
        for(i = 0; i < fhmm->K;i++){
                for(j = 0; j < 4;j++){
                        //fhmm->e[i][j] = fhmm->background[j];
                        fprintf(stdout,"%f ",fhmm->e[i][j]);
                }
                fprintf(stdout,"\n");
        }

        //exit(0);


        /*fhmm->e[o+0][0] = 1.0;
        fhmm->e[o+1][1] = 1.0;
        fhmm->e[o+2][2] = 1.0;
        fhmm->e[o+3][3] = 1.0;*/

        /*fhmm->tSN = 0.0f;
        fhmm->tNN = 0.0f;
        fhmm->tNB = 0.0f;

        fhmm->tBX = 0.0f;
        fhmm->tXE = 0.0f;

        fhmm->tEC = 0.0f;
        fhmm->tCC = 0.0f;
        fhmm->tCT = 0.0f;

        fhmm->tEJ = 0.0f;
        fhmm->tJJ = 0.0f;
        fhmm->tJB = 0.0f;
        */

        /* convert probs into log space/ set tindex to allow for fast-ish dyn
         * programming in case there is a sparse transition matrix */
        RUN(setup_model(fhmm));


        *f = fhmm;
        return OK;
ERROR:
        return FAIL;
}
