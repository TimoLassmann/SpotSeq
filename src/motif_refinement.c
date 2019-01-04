
#include "motif_refinement.h"
struct motif* init_motif(int W, int L);
int copy_motif(struct motif* source,struct motif* target);
void free_motif(struct motif* m);
int print_motif(struct motif* m);


int em_algorithm(double** counts,int W, int L, struct seq_buffer* sb)
{

        struct motif* tmp_motif = NULL;
        struct motif* best_motif = NULL;
        struct motif* org_motif = NULL;
        struct ihmm_sequence* s = NULL;
        uint8_t* tmp_seq;
        /* initialize - initial lambda is sqtr(N number oif sequences ) / n - number of overlapping kmers  */
        /* up to 1/ 2W  (2 times width of motif) */

        int N,n;



        double lambda = 0.0;
        double lambda_u = 0.0;
        double start_lambda = 0.0;

        double lambda_b = 0.0;
        double lambda_b_u = 0.0;
        double start_lambda_b = 0.0;


        double sum;
        double score_m;
        double score_b;
        double local_p;
        double likelihood = 0.0;
        double old_likelihood = 0.0;

        int i,j,c, offset, iter;
        RUNP(org_motif = init_motif(W, L));

        for(i = 0; i < W;i++){
                for(j = 0; j < L; j++){
                        org_motif->count_matrix[i][j] = counts[i][j];
                }
        }
        RUNP(best_motif = init_motif(W, L));
        RUNP(tmp_motif = init_motif(W, L));

        N = sb->num_seq;
        n = 0;
        sum = 0.0;
        for(i = 0; i < N;i++){
                s = sb->sequences[i];
                for(j = 0; j < s->seq_len; j++){
                        org_motif->background_counts[s->seq[j]]+= 1.0;
                        sum++;

                }
                n += s->seq_len - W;
        }
        //print_motif(org_motif);
        for(j = 0; j < L;j++){
                org_motif->background_freq[j] = org_motif->background_counts[j] / sum;
        }

        for(i = 0; i < W;i++){
                sum = 0.0;
                for(j = 0; j < L;j++){
                        sum += org_motif->count_matrix[i][j] + (double)N/ 1000000.0 * org_motif->background_freq[j];
                }
                for(j = 0; j < L;j++){
                        org_motif->freq_matrix[i][j] = (org_motif->count_matrix[i][j] + (double)N/ 1000000.0 * org_motif->background_freq[j]) / sum;
                }

        }


        print_motif(org_motif);

        start_lambda = sqrt((double)N) / (double)n;



        while (start_lambda <= (1.0 /( 2.0 *(double) W))){
                LOG_MSG("Testing lambda %f.",  start_lambda);

                lambda = start_lambda;
                likelihood = 0.0;
                sum = 0.0;
                for(i = 0; i < N;i++){
                        s = sb->sequences[i];

                        for(j = 0; j < s->seq_len-W; j++){
                                s->u[j] = lambda ;
                                sum += s->u[j];
                        }

                        for(j = 0; j < s->seq_len-W; j++){
                                //        s->u[j] /= sum;
                                // s->u[j] *= lambda ;
                        }


                }
                fprintf(stdout,"sum:%f n:%d\n", sum, n );


                RUN(copy_motif(org_motif, tmp_motif));


                /* iterations here */
                for(iter = 0; iter < 1000; iter++){
                        likelihood = 0.0;
                        lambda_u = 0.0;
                        lambda_b_u  = 0.0;
                        /* set counts to 0 */
                        for(j = 0; j < L;j++){
                                tmp_motif->background_counts[j] = 0.0;
                        }

                        for(i = 0; i < W;i++){
                                for(j = 0; j < L;j++){
                                        tmp_motif->count_matrix[i][j] = 0.0;
                                }

                        }

                        /* step 1 - estimate z_ij - store in sb->u  */
                        for(i = 0; i < N;i++){
                                s = sb->sequences[i];


                                sum = 0.0;
                                for(j = 0; j < s->seq_len- W; j++){
                                        /* calculate P(X_j | M) and P(X_j | B) */
                                        score_m = 1.0;
                                        score_b = 1.0;
                                        //score_m = prob2scaledprob(1.0);
                                        //score_b = prob2scaledprob(1.0);
                                        tmp_seq = s->seq +j;
                                        for(c = 0; c < W;c++){
                                                //score_m  += prob2scaledprob(tmp_motif->freq_matrix[c][tmp_seq[c]]);
                                                //score_b += prob2scaledprob(tmp_motif->background_freq[tmp_seq[c]]);
                                                score_m = score_m * tmp_motif->freq_matrix[c][tmp_seq[c]];
                                                score_b = score_b * tmp_motif->background_freq[tmp_seq[c]];

                                        }
                                        //score_m += prob2scaledprob(lambda);
                                        //score_b += prob2scaledprob(1.0-lambda);


                                        s->u[j] =  (score_m *  lambda) / ((score_m *  lambda)+(score_b* (1.0-lambda)));
                                        likelihood += s->u[j] * log(score_m *  lambda);
                                        likelihood += (1.0- s->u[j]) * log(score_b* (1.0-lambda));
                                        //likelihood += (1.0- s->u[j]) * (score_b);

                                        //if(!i && j < 20){
                                        //fprintf(stdout,"%f %f %f\n",log(score_m *  lambda),  log(score_b* (1.0-lambda)n),s->u[j] );
                                        //}

                                        sum += s->u[j];

                                        lambda_u += s->u[j];


                                }
                                //fprintf(stdout,"sum:%f %f lambda:%f\n", sum, sum / (double)n,lambda);
                                for(j = 0; j < s->seq_len - W; j++){
                                        //fprintf(stdout,"%d %f ", j,s->u[j]);
                                        //s->u[j] /= sum;
                                        //fprintf(stdout,"%f\n",s->u[j]);
                                        //s->u[j] *= lambda;

                                        ///fprintf(stdout,"%f %f\n",s->u[j], (1.0 - s->u[j]));
                                        //  lambda_u += s->u[j];
                                        // lambda_b_u += (1.0 - s->u[j]);



                                }
                                for(j = 0; j < s->seq_len- W; j++){
                                        tmp_seq = s->seq +j;
                                        for(c = 0; c < W;c++){
                                               tmp_motif->count_matrix[c][tmp_seq[c]] += s->u[j];
                                               tmp_motif->background_counts[tmp_seq[c]] += (1.0 - s->u[j]);
                                        }
                                }

                                //print_motif(tmp_motif);
                                //exit(0);

                                //if(iter == 10){
                                //       exit(0);
                                //}
                                /* smooth z... */

                                /*for(offset = 0;offset < W;offset++){
                                        for(j = offset; j < s->seq_len - (W*2);j += W){
                                                local_p = 0.0;
                                                for(c = 0; c < W; c++){
                                                        local_p += s->u[j+c];
                                                }

                                                if(local_p > 1.0f){
                                                        for(c = 0; c < W; c++){
                                                                if(s->u[j+c]){
                                                                        s->u[j+c] /= local_p;
                                                                }
                                                        }
                                                }
                                        }
                                }*/


                        }


                        if(fabs(likelihood - old_likelihood) < 10e-6){
                                break;

                        }
                        old_likelihood = likelihood;
                        /* update counts */
                        sum = 0.0;
                        for(j = 0; j < L;j++){
                                sum += tmp_motif->background_counts[j];
                        }

                        for(j = 0; j < L;j++){
                                tmp_motif->background_freq[j] = tmp_motif->background_counts[j] / sum;
                        }

                        for(i = 0; i < W;i++){
                                sum = 0.0;
                                for(j = 0; j < L;j++){
                                        sum += tmp_motif->count_matrix[i][j]  + (double)N/ 1000000.0 *  tmp_motif->background_freq[j];
                                }
                                for(j = 0; j < L;j++){
                                        tmp_motif->freq_matrix[i][j] = (tmp_motif->count_matrix[i][j]+ (double)N/ 1000000.0 *  tmp_motif->background_freq[j]) / sum;
                                        //tmp_motif->freq_matrix[i][j] = (tmp_motif->count_matrix[i][j]+ tmp_motif->background_freq[j])/ sum;//+ tmp_motif->background_freq[j]) / sum;
                                }


                        }


                        //print_motif(tmp_motif);

                        //lambda_u = lambda_u / (double)n;
                        //lambda_b_u  = lambda_b_u / (double)n;

                        lambda = lambda_u /  (double)n;//(lambda_u +lambda_b_u   );
                        LOG_MSG("iteration %d ll:%f lambde_u %f %f  new lambda = %f ",iter,likelihood, lambda_u ,lambda_b_u,lambda);
                        //print_motiflambda_u(tmp_motif);
                        //LOG_MSG("NEWL = %f",lambda );
                        //fprintf(st)
                }
                tmp_motif->log_likelihood = likelihood;
                LOG_MSG("ll:%f (%d)",likelihood, iter);
                if(best_motif->log_likelihood  <  tmp_motif->log_likelihood){

                        copy_motif(tmp_motif, best_motif);
                }

                start_lambda *= 2.0;


        }



        //print_motif(org_motif);

        print_motif(best_motif);


        for(i = 0; i < W;i++){
                for(j = 0; j < L; j++){
                        counts[i][j] = best_motif->count_matrix[i][j];
                }
        }

        free_motif(org_motif);
        free_motif(best_motif);
        free_motif(tmp_motif);

        /* step 0 - estimate z_ij - store in sb->u  */

        /* 2nd step estimate lambda (mixing factor) */

        /* 3rd step  */
        return OK;
ERROR:
        return FAIL;
}


struct motif* init_motif(int W, int L)
{

        struct motif* m = NULL;
        int i;
        MMALLOC(m, sizeof(struct motif));
        m->L = L;
        m->W = W;
        m->log_likelihood = -INFINITY;
        m->background_counts = NULL;
        m->background_freq = NULL;
        m->count_matrix = NULL;
        m->freq_matrix =  NULL;
        m->count_matrix = galloc(m->count_matrix,W,L,0.0);
        m->freq_matrix = galloc(m->freq_matrix,W,L,0.0);

        MMALLOC(m->background_counts, sizeof(double) * L);
        MMALLOC(m->background_freq, sizeof(double) * L);

        for(i = 0; i < L;i++){
                m->background_counts[i] = 0.0;
                m->background_freq[i] = 0.0;

        }

        return m;
ERROR:
        free_motif(m);
        return NULL;
}

int print_motif(struct motif* m)
{
        int i,j;
        double sum = 0.0;
        fprintf(stdout,"Log-likelihood:%f\n", m->log_likelihood);
        fprintf(stdout,"Background\n");

        for(j = 0; j < m->L;j++){
                fprintf(stdout,"%3.3f\t", m->background_counts[j]);

        }
        for(j = 0; j < m->L;j++){

                fprintf(stdout,"%3.3f\t",m->background_freq[j]);
        }
        fprintf(stdout,"\n");
        fprintf(stdout,"\n");

        for(i = 0; i < m->W;i++){
                sum = 0.0;
                for(j = 0; j < m->L;j++){
                        fprintf(stdout,"%3.3f\t",m->count_matrix[i][j]);
                        sum += m->count_matrix[i][j];

                }
                fprintf(stdout,"SUM:%3.3f\t",sum);
                for(j = 0; j < m->L;j++){

                        fprintf(stdout,"%3.3f\t",m->freq_matrix[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        return OK;

}

int copy_motif(struct motif* source,struct motif* target)
{
        int i,j,L,W;
        ASSERT(source->L == target->L , "Different alphabet len!");
        ASSERT(source->W == target->W , "Different motif len!");
        L = source->L;
        W = source->W;
        target->log_likelihood = source->log_likelihood;
        for(j = 0; j < L;j++){
                target->background_counts[j] =  source->background_counts[j];
                target->background_freq[j] = source->background_freq[j];

         }

        for(i = 0; i < W;i++){
                for(j = 0; j < L;j++){
                        target->count_matrix[i][j] = source->count_matrix[i][j];
                        target->freq_matrix[i][j] = source->freq_matrix[i][j];
                }
        }

        return OK;
ERROR:
        return FAIL;
}

void free_motif(struct motif* m)
{
        if(m){
                if(m->count_matrix){
                        gfree(m->count_matrix);
                }
                if(m->freq_matrix){
                        gfree(m->freq_matrix);
                }
                MFREE(m);
        }
}
