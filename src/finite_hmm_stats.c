/* These functions are from HMMER / EASEL  */

#include "tldevel.h"
#include "tlrng.h"
#include "finite_hmm_struct.h"
#include "finite_hmm_alloc.h"
#include "finite_hmm.h"

#define FINITE_HMM_STATS_IMPORT
#include "finite_hmm_stats.h"


#define eslCONST_PI    3.14159265358979323846264338328


static int esl_gumbel_FitComplete(double *x, int n, double *ret_mu, double *ret_lambda);

static int esl_stats_DMean(const double *x, int n, double *opt_mean, double *opt_var);
static int lawless416(double *x, int n, double lambda, double *ret_f, double *ret_df);
static double esl_gumbel_invcdf(double p, double mu, double lambda);

int fhmm_calibrate(struct fhmm* fhmm,struct fhmm_dyn_mat* dm, int seed)
{
        double* f_scores = NULL;
        struct rng_state* rng = NULL;
        uint8_t* seq = NULL;
        int malloc_len = 10000;


        int i,j;//c;
        //double avg;

        double tailp = 0.04;
        int sim_len = 200;
        int sim_N = 200;
        double score;
        double   gmu, glam;
        double P;
        FILE* f_ptr = NULL;





        RUNP(rng = init_rng(seed));

        MMALLOC(seq, sizeof(uint8_t) * malloc_len);

        MMALLOC(f_scores, sizeof(double) * sim_N);
        malloc_len = sim_len;
        /* rubbish! ... */
        if(malloc_len >= dm->alloc_matrix_len){
                resize_fhmm_dyn_mat(dm, malloc_len, fhmm->K);
        }

        for(i = 0;i < sim_N;i++){
                //sim_len = ;
                //sim_len = tl_random_gaussian(rng, 300,150);
                //sim_len = MACRO_MAX(50, sim_len);
                //sim_len = MACRO_MIN(5000, sim_len);

                for(j = 0; j < sim_len;j++){
                        seq[j] = tl_random_int(rng,fhmm->L);
                }
                score_seq_fwd(fhmm, dm, seq, sim_len, 1, &score, &P);
                //LOG_MSG("%f %f", score,P);
                //configure_target_len(fhmm, sim_len, 1);
                //forward(fhmm, dm, &fhmm->f_score, seq, sim_len);
                //random_model_score(sim_len,&fhmm->r_score);// ,seq, len,len );
                f_scores[i] = score;
                //LOG_MSG("%f %f -> %f ", fhmm->f_score,fhmm->r_score, f_scores[i]);


        }
        RUN(esl_gumbel_FitComplete(f_scores, sim_N, &gmu, &glam));

        //fprintf(stdout,"lambda:%f  tau: %f\n", gmu,glam);


        /* Explanation of the eqn below: first find the x at which the Gumbel tail
         * mass is predicted to be equal to tailp. Then back up from that x
         * by log(tailp)/lambda to set the origin of the exponential tail to 1.0
         * instead of tailp.
         */

        fhmm->tau = esl_gumbel_invcdf(1.0-tailp, gmu, glam) + (log(tailp) / fhmm->lambda);

        //fprintf(stdout,"lambda:%f  tau: %f\n", fhmm->lambda,fhmm->tau);
        /*RUNP(f_ptr = fopen("scores.csv", "w"));
        sim_N = 10000;
        for(i = 0;i < sim_N;i++){
                sim_len = 100;
                for(j = 0; j < sim_len;j++){
                        seq[j] = tl_random_int(rng,fhmm->L);
                }
                score_seq_fwd(fhmm, dm, seq, sim_len, 1, &score, &P);
                P = exp(P);
                fprintf(f_ptr,"%f,%f,%f,",score ,P, P* (double) sim_N);
                sim_len = 400;
                for(j = 0; j < sim_len;j++){
                        seq[j] = tl_random_int(rng,fhmm->L);
                }
                score_seq_fwd(fhmm, dm, seq, sim_len, 1, &score, &P);
                P = exp(P);
                fprintf(f_ptr,"%f,%f,%f,",score ,P, P* (double) sim_N);
                sim_len = 1600;
                for(j = 0; j < sim_len;j++){
                        seq[j] = tl_random_int(rng,fhmm->L);
                }
                score_seq_fwd(fhmm, dm, seq, sim_len, 1, &score, &P);
                P = exp(P);
                fprintf(f_ptr,"%f,%f,%f\n",score ,P, P* (double) sim_N);
//fprintf(f_ptr,"%f,%f\n",score,P);

                //if(P <= 0.05){
                        //LOG_MSG("%f %f -> %f p:%f", fhmm->f_score,fhmm->r_score, score,P);
                //}
        }
        fclose(f_ptr);*/
        MFREE(f_scores);
        MFREE(seq);

        return OK;
ERROR:
        if(seq){
                MFREE(seq);
        }
        return FAIL;
}





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
                //ESL_DPRINTF1(("esl_gumbel_FitComplete(): Newton/Raphson failed; switchover to bisection"));

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
                if (i == 100)
                        ERROR_MSG( "Even bisection search failed in esl_gumbel_FitComplete().");

                lambda = mid;
        }

        /* 3. Substitute into Lawless 4.1.5 to find mu
         */
        esum = 0.;
        for (i = 0; i < n; i++)
                esum  += exp(-lambda * x[i]);
        mu = -log(esum / n) / lambda;

        *ret_lambda = lambda;
        *ret_mu     = mu;
        return OK;
ERROR:
        return FAIL;
}

int lawless416(double *x, int n, double lambda, double *ret_f, double *ret_df)
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
        return OK;
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
int
esl_stats_DMean(const double *x, int n, double *opt_mean, double *opt_var)
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

/* Function:  esl_gumbel_invcdf()
 *
 * Purpose:   Calculates the inverse CDF for a Gumbel distribution
 *            with parameters <mu> and <lambda>. That is, returns
 *            the quantile <x> at which the CDF is <p>.
 */
double
esl_gumbel_invcdf(double p, double mu, double lambda)
{
    return mu - ( log(-1. * log(p)) / lambda);
}

