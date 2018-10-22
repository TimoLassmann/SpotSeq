/* The code for fitting Gumbel (type I extreme value) distributions */
/* was taken from the easel library: */
/* https://github.com/EddyRivasLab/easel */



#include "calibrate.h"

#include "ihmm_seq.h"
#include "emit_random.h"
#include "run_score.h"
#include "finite_hmm.h"

#define eslCONST_PI    3.14159265358979323846264338328

static void lawless416(double *x, int n, double lambda, double *ret_f, double *ret_df);
static int esl_stats_DMean(const double *x, int n, double *opt_mean, double *opt_var);
static double esl_gumbel_invcdf(double p, double mu, double lambda);

static int esl_gumbel_FitComplete(double *x, int n, double *ret_mu, double *ret_lambda);

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
        RUNP(sb_in = get_sequences_from_hdf5_model(model_file));

        RUNP(sb = emit_sequences_from_random_model(sb_in, num_seq));
        free_ihmm_sequences(sb_in);
        sb_in = NULL;

        //sb = emit_sequences_from_random_model(struct seq_buffer *sb_in, int num)

        /* Step 2: score random sequence, get scores into an array */

        LOG_MSG("Loading model.");
        RUNP(fhmm = init_fhmm(model_file));


        RUN(run_score_sequences(fhmm, sb,num_threads));

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
        FILE* fout_ptr = NULL;
        double mu, lambda;

        int i,j;
        double  est_mu,est_lambda;
        double* scores = NULL;
        double r;
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
int esl_gumbel_FitComplete(double *x, int n, double *ret_mu, double *ret_lambda)
{

        double  variance;
        double  lambda, mu;
        double  fx;			/* f(x)  */
        double  dfx;			/* f'(x) */
        double  esum;                 /* \sum e^(-lambda xi) */
        double  tol = 1e-5;
        int     i;

        /* 1. Find an initial guess at lambda
         *    (Evans/Hastings/Peacock, Statistical Distributions, 2000, p.86)
         */
        esl_stats_DMean(x, n, NULL, &variance);
        lambda = eslCONST_PI / sqrt(6.*variance);

        /* 2. Use Newton/Raphson to solve Lawless 4.1.6 and find ML lambda
         */
        for (i = 0; i < 100; i++)
        {
                lawless416(x, n, lambda, &fx, &dfx);
                if (fabs(fx) < tol) break;             /* success */
                lambda = lambda - fx / dfx;	     /* Newton/Raphson is simple */
                if (lambda <= 0.) lambda = 0.001;      /* but be a little careful  */
        }

        /* 2.5: If we did 100 iterations but didn't converge, Newton/Raphson failed.
         *      Resort to a bisection search. Worse convergence speed
         *      but guaranteed to converge (unlike Newton/Raphson).
         *      We assume that fx is a monotonically decreasing function of x;
         *      i.e. fx > 0 if we are left of the root, fx < 0 if we
         *      are right of the root.
         */
        if (i == 100)
        {
                double left, right, mid;
                DPRINTF1(("esl_gumbel_FitComplete(): Newton/Raphson failed; switchover to bisection"));

                /* First bracket the root */
                left  = 0.;	                 	/* for sure */
                right = eslCONST_PI / sqrt(6.*variance);  /* an initial guess */
                lawless416(x, n, lambda, &fx, &dfx);
                while (fx > 0.)
                {
                        right *= 2.;		/* arbitrary leap to the right */
                        if (right > 100.) /* no reasonable lambda should be > 100, we assert */
                              ERROR_MSG("Failed to bracket root in esl_gumbel_FitComplete().");
                        lawless416(x, n, right, &fx, &dfx);
                }

                /* Now, bisection search in left/right interval */
                for (i = 0; i < 100; i++)
                {
                        mid = (left + right) / 2.;
                        lawless416(x, n, mid, &fx, &dfx);
                        if (fabs(fx) < tol) break;             /* success */
                        if (fx > 0.)	left = mid;
                        else          right = mid;
                }
                if (i == 100){
                        ERROR_MSG("Even bisection search failed in esl_gumbel_FitComplete().");
                }

                lambda = mid;
        }

        /* 3. Substitute into Lawless 4.1.5 to find mu
         */
        esum = 0.;
        for (i = 0; i < n; i++){
                esum  += exp(-lambda * x[i]);
        }
        mu = -log(esum / n) / lambda;

        *ret_lambda = lambda;
        *ret_mu     = mu;
        return OK;
ERROR:
        return FAIL;
}

/*****************************************************************
 * Complete data, maximum a posteriori parameters
 *****************************************************************/

/* lawless416()b
 * SRE, Thu Nov 13 11:48:50 1997 [St. Louis]
 *
 * Purpose:  Equation 4.1.6 from [Lawless82], pg. 143, and
 *           its first derivative with respect to lambda,
 *           for finding the ML fit to Gumbel lambda parameter.
 *           This equation gives a result of zero for the maximum
 *           likelihood lambda.
 *
 * Args:     x      - array of sample values
 *           n      - number of samples
 *           lambda - a lambda to test
 *           ret_f  - RETURN: 4.1.6 evaluated at lambda
 *           ret_df - RETURN: first derivative of 4.1.6 evaluated at lambda
 *
 * Return:   (void)
 */
static void lawless416(double *x, int n, double lambda, double *ret_f, double *ret_df)
{
        double esum;			/* \sum e^(-lambda xi)      */
        double xesum;			/* \sum xi e^(-lambda xi)   */
        double xxesum;		/* \sum xi^2 e^(-lambda xi) */
        double xsum;			/* \sum xi                  */
        int i;

        esum = xesum = xsum  = xxesum = 0.;
        for (i = 0; i < n; i++)
        {
                xsum   += x[i];
                xesum  += x[i] * exp(-1. * lambda * x[i]);
                xxesum += x[i] * x[i] * exp(-1. * lambda * x[i]);
                esum   += exp(-1. * lambda * x[i]);
        }
        *ret_f  = (1./lambda) - (xsum / n)  + (xesum / esum);
        *ret_df = ((xesum / esum) * (xesum / esum))
                - (xxesum / esum)
                - (1. / (lambda * lambda));
}

/*****************************************************************
 * 1. Summary statistics calculations (means, variances)
 *****************************************************************/

/* Function:  esl_stats_DMean()
 * Synopsis:  Calculates mean and $\sigma^2$ for samples $x_i$.
 *
 * Purpose:   Calculates the sample mean and $s^2$, the unbiased
 *            estimator of the population variance, for a
 *            sample of <n> numbers <x[0]..x[n-1]>, and optionally
 *            returns either or both through <ret_mean> and
 *            <ret_var>.
 *
 *            <esl_stats_FMean()> and <esl_stats_IMean()> do the same,
 *            for float and integer vectors.
 *
 * Args:      x        - samples x[0]..x[n-1]
 *            n        - number of samples
 *            opt_mean - optRETURN: mean
 *            opt_var  - optRETURN: estimate of population variance
 *
 * Returns:   <eslOK> on success.
 */
int esl_stats_DMean(const double *x, int n, double *opt_mean, double *opt_var)
{
        double sum   = 0.;
        double sqsum = 0.;
        int i;

        for (i = 0; i < n; i++)
        {
                sum   += x[i];
                sqsum += x[i]*x[i];
        }
        if (opt_mean != NULL)  *opt_mean = sum / (double) n;
        if (opt_var  != NULL)  *opt_var  = (sqsum - sum*sum/(double)n) / ((double)n-1);
        return OK;
}
