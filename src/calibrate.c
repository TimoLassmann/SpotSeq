/* The code for fitting Gumbel (type I extreme value) distributions */
/* was taken from the easel library: */
/* https://github.com/EddyRivasLab/easel */



#include "calibrate.h"

#include "emit_random.h"
#include "run_score.h"
#include "finite_hmm.h"

#define eslCONST_PI    3.14159265358979323846264338328

static void lawless416(double *x, int n, double lambda, double *ret_f, double *ret_df);
static int esl_stats_DMean(const double *x, int n, double *opt_mean, double *opt_var);
static double esl_gumbel_invcdf(double p, double mu, double lambda);

//int esl_gumbel_FitComplete(double *x, int n, double *ret_mu, double *ret_lambda);

int calibrate(char* model_file,int num_threads, double* mu, double* lambda)
{
        FILE* fout_ptr = NULL;

        struct seq_buffer* sb = NULL;
        struct seq_buffer* sb_in = NULL;
        struct fhmm* fhmm = NULL;

        int num_seq = 100000;
        double* scores = NULL;

        int i;

        /* Step 1: emit sequences with same background and length
         * distribution as training sequences */
        LOG_MSG("Generating random sequences.");
        RUNP(sb_in = get_sequences_from_hdf5_model(model_file, IHMM_SEQ_READ_ONLY_SEQ));

        RUNP(sb = emit_sequences_from_random_model(sb_in, num_seq,0));
        free_ihmm_sequences(sb_in);
        sb_in = NULL;

        //sb = emit_sequences_from_random_model(struct seq_buffer *sb_in, int num)

        /* Step 2: score random sequence, get scores into an array */

        LOG_MSG("Loading model.");
        RUNP(fhmm = init_fhmm(model_file));


        //RUN(run_score_sequences(fhmm, sb,num_threads));

        /* copy scores over */
        MMALLOC(scores, sizeof(double)* num_seq);
        for(i = 0; i < sb->num_seq;i++){
                scores[i] = sb->sequences[i]->score;
        }

        /* Step 3: fit gumbel distribution  */
        RUN(esl_gumbel_FitComplete(scores, num_seq, mu, lambda));
        fout_ptr = fopen("scores.csv", "w");
        fprintf(fout_ptr,"Name,Scores\n");
        for(i = 0; i < sb->num_seq;i++){
                fprintf(fout_ptr,"%d,%f\n",i,scores[i]);
        }
        fclose(fout_ptr);
        fprintf(stdout,"MU:%f\nLAMBDA:%f\n",*mu,*lambda);

        free_ihmm_sequences(sb_in);
        free_ihmm_sequences(sb);
        MFREE(scores);

        return OK;
ERROR:
        free_ihmm_sequences(sb);
        free_ihmm_sequences(sb_in);
        MFREE(scores);
        return FAIL;
}

//#ifdef ITEST

int main (int argc,char * argv[])
{
        //FILE* fout_ptr = NULL;
        double mu, lambda;

        //int i,j;
        //double  est_mu,est_lambda;
        //double* scores = NULL;
        //double r;
        rk_state rndstate;

        RUN(rk_randomseed(&rndstate));

        fprintf(stdout,"Running libhdf5glue sanity tests\n");
        RUN(calibrate("test.h5", 8, &mu, &lambda));

        /*lambda = 0.4;
        mu = -20;
        MMALLOC(scores, sizeof(double) * 10000);

        fout_ptr = fopen("scores.csv", "w");
        fprintf(fout_ptr,"Name,Scores\n");
        for(i = 0; i < 10000;i++){
                r = rk_double(&rndstate);
                scores[i] = esl_gumbel_invcdf(r,mu,lambda);
                fprintf(fout_ptr,"%d,%f\n",i,scores[i]);
        }
        fclose(fout_ptr);
        fprintf(stdout,"MU:%f\nLAMBDA:%f\n",mu,lambda);
        RUN(esl_gumbel_FitComplete(scores,10000,  &est_mu, &est_lambda));

        fprintf(stdout,"MU:%f\nLAMBDA:%f\n",est_mu,est_lambda);
        MFREE(scores);*/

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


//#endif


/* Function:  esl_gumbel_invcdf()
 *
 * Purpose:   Calculates the inverse CDF for a Gumbel distribution
 *            with parameters <mu> and <lambda>. That is, returns
 *            the quantile <x> at which the CDF is <p>.
 */
double esl_gumbel_invcdf(double p, double mu, double lambda)
{
    return mu - ( log(-1. * log(p)) / lambda);
}


/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 *
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *
 * SVN $Id: esl_stats.c 854 2013-02-25 22:00:19Z wheelert $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/hmmer/3.1/esl_stats.c $
 *****************************************************************/


/* Function: esl_gumbel_FitComplete()
 * Synopsis: Estimates $\mu$, $\lambda$ from complete data.
 * Date:     SRE, Fri Nov 14 07:56:29 1997 [St. Louis] - HMMER's EVDMaxLikelyFit()
 *
 * Purpose:  Given an array of Gumbel-distributed samples <x[0]..x[n-1]>,
 *           find maximum likelihood parameters <mu> and <lambda>.
 *
 * Algorithm: Uses approach described in [Lawless82]. Solves
 *            for lambda using Newton/Raphson iterations,
 *            then substitutes lambda into Lawless' equation 4.1.5
 *            to get mu.
 *
 * Args:     x          - list of Gumbel distributed samples
 *           n          - number of samples
 *           ret_mu     : RETURN: ML estimate of mu
 *           ret_lambda : RETURN: ML estimate of lambda
 *
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslENOHALT> if the fit doesn't converge.
 */
