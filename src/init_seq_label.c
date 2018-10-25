#include "init_seq_label.h"

int label_seq_based_on_random_fhmm(struct seq_buffer* sb, int k, double alpha)
{

        struct fhmm* fhmm = NULL;
        int i,j;
        rk_state rndstate;
        float sum;
        int expected_len;
        float s;                /* self transition */
        float e;                /* 1- r (exit) */
        float pseudocount = 0.0f; /* just to make sure everything is
                                   * connected; otherwise there is a
                                   * chance that the hmm model is
                                   * invalid (no path from start to
                                   * end */



        init_logsum();
        ASSERT(sb!= NULL, "No seq buffer");
        ASSERT(sb->num_seq != 0, "No sequences");

        rk_randomseed(&rndstate);

        RUNP(fhmm = alloc_fhmm());

        /* generate  HMM parameters - output hasto be a t and e matrix */

        fhmm->K = k + 2;
        fhmm->L = sb->L;
        MMALLOC(fhmm->background, sizeof(float)* sb->L);
        for(i = 0; i < sb->L;i++){
                fhmm->background[i] = sb->background[i];
        }


        fhmm->e = malloc_2d_float(fhmm->e, fhmm->K, fhmm->L, 0.0);
        fhmm->t = malloc_2d_float(fhmm->t, fhmm->K, fhmm->K, 0.0);
        /* Fill with random counts... */
        /* First emission (cos easy...) */

        for(i = 0; i < fhmm->K;i++){
                sum = 0.0f;
                for(j = 0 ; j < fhmm->L;j++){
                        fhmm->e[i][j] = rk_gamma(&rndstate,alpha, 1.0) + pseudocount;
                        fprintf(stdout,"%f %d %d %f\n",fhmm->e[i][j],i,j,rk_gamma(&rndstate,alpha, 1.0));
                        sum += fhmm->e[i][j];
                }
                ASSERT(sum != 0.0, "No sum - sampling gamma failed.");
                for(j = 0 ; j < fhmm->L;j++){
                        fhmm->e[i][j] /= sum;
                }
        }
        //exit(0);

        /* Second transitions */
        /* from start state */
        fhmm->t[IHMM_START_STATE][IHMM_START_STATE] = 0.0;
        fhmm->t[IHMM_START_STATE][IHMM_END_STATE] = 0.0f;
        sum = 0.0f;
        for(i = 2; i < fhmm->K;i++){
                fhmm->t[IHMM_START_STATE][i] = rk_gamma(&rndstate,alpha, 1.0) + pseudocount;
                sum += fhmm->t[IHMM_START_STATE][i];
        }
        ASSERT(sum != 0.0, "No sum - sampling gamma failed.");
        for(i = 2; i < fhmm->K;i++){
                fhmm->t[IHMM_START_STATE][i] /= sum;
        }


        /* stop state - all zero  */
        for(i = 0; i < fhmm->K;i++){
                fhmm->t[IHMM_END_STATE][i] = 0.0f;
        }

        /* remaining states  */
        /* first I calculate the optimal S -> end state transition,
         * then weave this into the other parameters */
        expected_len = 0;
        for(i = 0; i < sb->num_seq;i++){
                expected_len += sb->sequences[i]->seq_len;
        }
        expected_len = expected_len / sb->num_seq;

        /* initialise transitions  */
        s = (double)expected_len / ((double) expected_len + 1.0);
        e = 1.0 -s;
        for(i = 2; i < fhmm->K;i++){
                sum = 0;
                for(j = 2; j < fhmm->K;j++){
                        fhmm->t[i][j] = rk_gamma(&rndstate,alpha, 1.0)+  pseudocount;
                        sum += fhmm->t[i][j];
                }
                ASSERT(sum != 0.0, "No sum - sampling gamma failed.");
                for(j = 2; j < fhmm->K;j++){
                        fhmm->t[i][j] = fhmm->t[i][j] /sum * s;
                }
                fhmm->t[i][IHMM_END_STATE] = e;

        }
        fprintf(stdout,"Emissions\n");

        for(i = 0; i < fhmm->K;i++){
                sum = 0.0;
                for(j = 0; j < fhmm->L;j++){
                        fprintf(stdout,"%0.3f ", fhmm->e[i][j]);
                        sum +=  fhmm->e[i][j];
                }
                fprintf(stdout,"%f\n",sum);
        }
        fprintf(stdout,"Transitions\n");
        for(i = 0; i < fhmm->K;i++){
                sum = 0.0;
                for(j = 0; j < fhmm->K;j++){
                        fprintf(stdout,"%0.3f ", fhmm->t[i][j]);
                        sum +=  fhmm->t[i][j];
                }
                fprintf(stdout,"%f\n",sum);
        }


        /* alloc dyn matrices (now that I know how many states there are) */

        RUN(alloc_dyn_matrices(fhmm));
        RUN(realloc_dyn_matrices(fhmm, sb->max_len));

        /* convert probs into log space/ set tindex to allow for fast-ish dyn
         * programming in case there is a sparse transition matrix */
        RUN(setup_model(fhmm));

        /* Now I *should* be ready to run forward/backward on all sequences.. */

        RUN(run_label_sequences(fhmm,sb, 8));

        free_fhmm(fhmm);
        return OK;
ERROR:
        free_fhmm(fhmm);
        return FAIL;
}
