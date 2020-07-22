
#include "finite_hmm.h"
#include "finite_hmm_alloc.h"

#include "tllogsum.h"


#define eslCONST_LOG2  0.69314718055994529

double esl_exp_surv(double x, double mu, double lambda);
double esl_exp_logsurv(double x, double mu, double lambda);


static float rel_entropy(float* vec, float* back,int N);

int score_seq_fwd(struct fhmm* fhmm , struct fhmm_dyn_mat* m, uint8_t* a, int len,int mode, double* ret_score,double*P)
{
        float f,r;
        double s;
        //configure_target_len(fhmm,len, mode);
        forward(fhmm, m, &f, a, len,mode);
        random_model_score(len,&r);// ,seq, len,len );
        s =  (f - r) / eslCONST_LOG2;
        *P = esl_exp_logsurv(s, fhmm->tau,fhmm->lambda);
        *ret_score = s;
        //fprintf(f_ptr,"%f,%f,%f\n",score,P, P* (double) sim_N);
        return OK;
}


/* Calculate the bayesian information criteria score for a finite HMM */
/* ML is the maximum likelihood (i.e. product of p(x^k | M ))  */
/* data is the number of residues in the training dataset */
int calculate_BIC( struct fhmm* fhmm, double ML, double data,double* BIC)
{
        int i,j,c;
        double num_param;

        ASSERT(fhmm != NULL, "No model");

        ASSERT(data != 0, "No data");

        num_param = 0.0;
        for(i = 0; i < fhmm->K;i++){
                c = 0;
                for(j = 0; j < fhmm->K;j++){
                        if(fhmm->t[i][j] != prob2scaledprob(0.0)){
                                c++;
                        }
                }
                c = c - 1;       /* we need to subtract one because the last transition is not a free parameter - the sumhas to be 1 */
                num_param += c;
        }
        for(i = 0; i < fhmm->K;i++){
                c = 0;
                for(j = 0; j < fhmm->L;j++){
                        if( fhmm->e[i][j] != prob2scaledprob(0.0)){
                                c++;
                        }
                }
                 c = c - 1;       /* we need to subtract one because the last emission is not a free parameter - the sumhas to be 1 */
                 num_param += c;
        }

        *BIC = -2.0 * ML + num_param * log(data);

        return OK;
ERROR:
        return FAIL;


}

/* Function:  esl_exp_surv()
 *
 * Purpose:   Calculates the survivor function, $P(X>x)$ (that is, 1-CDF,
 *            the right tail probability mass) for an exponential distribution,
 *            given value <x>, offset <mu>, and decay parameter <lambda>.
 */
double
esl_exp_surv(double x, double mu, double lambda)
{
  if (x < mu) return 1.0;
  return exp(-lambda * (x-mu));
}

/* Function:  esl_exp_logsurv()
 *
 * Purpose:   Calculates the log survivor function, $\log P(X>x)$ (that is,
 *            log(1-CDF), the log of the right tail probability mass) for an
 *            exponential distribution, given value <x>, offset <mu>, and
 *            decay parameter <lambda>.
 */
double
esl_exp_logsurv(double x, double mu, double lambda)
{
  if (x < mu) return 0.0;
  return -lambda * (x-mu);
}


int random_model_score(int len,float* ret_score)
{

        float  r;                /* self transition */
        float  e;                /* 1- r (exit) */

        ASSERT(len > 0, "Seq is of length 0");

        /* initialise transitions  */
        r = (float )len / ((float ) len + 1.0);
        e = 1.0 -r;
        r = prob2scaledprob(r);
        e = prob2scaledprob(e);

        *ret_score = (float) len * r + e;

        return OK;
ERROR:
        return FAIL;
}

int forward(struct fhmm* fhmm , struct fhmm_dyn_mat* m, float* ret_score, uint8_t* a, int len, int mode)
{
        int i,j,c,f;

        float** NBECJ = NULL;
        float** matrix = NULL;

        //const float* trans = 0;
        float tSN;
        float tNN;
        float tNB;

        float tBX;
        float tXE;

        float tEC;
        float tCC;
        float tCT;

        float tEJ;
        float tJJ;
        float tJB;
        float p,q;
        //float tmp = 0;




        ASSERT(fhmm != NULL, "No model");
        ASSERT(m != NULL, "No dyn programming  matrix");
        ASSERT(a != NULL, "No sequence");
        ASSERT(len > 0, "Seq is of length 0");

        //ASSERT(len == fhmm->config_len, "Model configured for len %d but seqlen is %d. ", fhmm->config_len,len);

        matrix = m->F_matrix;
        NBECJ = m->F_NBECJ;



        if(mode){
                q = 0.5f;
                p = (float) len / ((float)len + 3.0f);
        }else{
                q = 0.0f;
                p = (float) len / ((float)len + 2.0f);

        }

        tSN = prob2scaledprob(1.0f);
        tNN = prob2scaledprob(p);
        tNB = prob2scaledprob(1.0f - p);
        tBX = prob2scaledprob(2.0f / (float) (fhmm->K * ( fhmm->K + 1.0f)));
        tXE = prob2scaledprob(1.0f);
        //LOG_MSG("%f", scaledprob2prob(fhmm->tBX));
        tEC = prob2scaledprob(1.0f - q);
        tCC = prob2scaledprob(p);
        tCT = prob2scaledprob(1.0f - p);

        tEJ = prob2scaledprob(q);
        tJJ = prob2scaledprob(p);
        tJB = prob2scaledprob(1.0f - p);

        for(j = 0; j < fhmm->K;j++){
                matrix[0][j] = -INFINITY;
        }


        NBECJ[0][N_STATE] = prob2scaledprob(1.0f) + tSN; /* start -> N state has 1.0 prob */
        NBECJ[0][B_STATE] = tNB;
        NBECJ[0][E_STATE] = -INFINITY;
        NBECJ[0][C_STATE] = -INFINITY;
        NBECJ[0][J_STATE] = -INFINITY;

        /* fprintf(stdout,"\n"); */
        /* fprintf(stdout,"%0.2f\t%0.2f: \t", NBECJ[0][N_STATE],  NBECJ[0][B_STATE]); */
        /* for(j = 0; j < fhmm->K;j++){ */
        /*         fprintf(stdout,"%0.4f ",matrix[0][j]); */
        /* } */
        /* fprintf(stdout,": %0.2f\t%0.2f\t%0.2f\n", NBECJ[0][E_STATE],  NBECJ[0][J_STATE], NBECJ[0][C_STATE]); */


        for(i = 1; i < len+1;i++){
                for(j = 0; j < fhmm->K;j++){
                        matrix[i][j] = -INFINITY;
                }
                for(j = 0; j < 5;j++){
                        NBECJ[i][j] = -INFINITY;
                }

                for(j = 0; j < fhmm->K;j++){
                        for(c = 1; c < fhmm->tindex[j][0];c++){
                                f = fhmm->tindex[j][c];
                                matrix[i][f] = logsum(matrix[i][f], matrix[i-1][j] + fhmm->t[j][f]);
                        }
                }
                for(j = 0;j < fhmm->K;j++){
                        /* add transition from B state */
                        matrix[i][j] = logsum(matrix[i][j], NBECJ[i-1][B_STATE] + tBX);
                        matrix[i][j] += fhmm->e[j][a[i-1]];
                        //cur[c] += fhmm->e[c][a[i-1]];
                        NBECJ[i][E_STATE] = logsum(NBECJ[i][E_STATE], matrix[i][j] + tXE);
                }
                /* J */
                NBECJ[i][J_STATE] = logsum(NBECJ[i-1][J_STATE] + tJJ, NBECJ[i][E_STATE] + tEJ);
                /* C */
                NBECJ[i][C_STATE] = logsum(NBECJ[i-1][C_STATE] + tCC, NBECJ[i][E_STATE] + tEC);
                /* N */
                NBECJ[i][N_STATE] = NBECJ[i-1][N_STATE] + tNN;
                /* B */
                NBECJ[i][B_STATE] = logsum(NBECJ[i][N_STATE] + tNB, NBECJ[i][J_STATE]+ tJB);


                /* fprintf(stdout,"%0.2f\t%0.2f: \t", NBECJ[i][N_STATE],  NBECJ[i][B_STATE]); */
                /* for(j = 0; j < fhmm->K;j++){ */
                /*         fprintf(stdout,"%0.4f\t",matrix[i][j]); */
                /* } */
                /* fprintf(stdout,": %0.2f\t%0.2f\t%0.2f\n", NBECJ[i][E_STATE],  NBECJ[i][J_STATE], NBECJ[i][C_STATE]); */
        }
        /* LOG_MSG("%f %f %f ",NBECJ[len][C_STATE] , fhmm->tCT,NBECJ[len][C_STATE] + fhmm->tCT); */
        *ret_score = (NBECJ[len][C_STATE] + tCT);
        return OK;
ERROR:
        return FAIL;
}

int backward(struct fhmm* fhmm,struct fhmm_dyn_mat* m , float* ret_score, uint8_t* a, int len,int mode)
{
        float tSN;
        float tNN;
        float tNB;

        float tBX;
        float tXE;

        float tEC;
        float tCC;
        float tCT;

        float tEJ;
        float tJJ;
        float tJB;
        float p,q;

        int i,j,c,f;
        float** matrix = NULL;
        float** NBECJ = NULL;
        ASSERT(fhmm != NULL, "No model");
        ASSERT(m != NULL, "No dyn programming  matrix");
        ASSERT(a != NULL, "No sequence");

        ASSERT(len > 0, "Seq is of length 0");
        //ASSERT(len == fhmm->config_len, "Model configured for len %d but seqlen is %d. ", fhmm->config_len,len);


        if(mode){
                q = 0.5f;
                p = (float) len / ((float)len + 3.0f);
        }else{
                q = 0.0f;
                p = (float) len / ((float)len + 2.0f);

        }

        tSN = prob2scaledprob(1.0f);
        tNN = prob2scaledprob(p);
        tNB = prob2scaledprob(1.0f - p);
        tBX = prob2scaledprob(2.0f / (float) (fhmm->K * ( fhmm->K + 1.0f)));
        tXE = prob2scaledprob(1.0f);
        //LOG_MSG("%f", scaledprob2prob(fhmm->tBX));
        tEC = prob2scaledprob(1.0f - q);
        tCC = prob2scaledprob(p);
        tCT = prob2scaledprob(1.0f - p);

        tEJ = prob2scaledprob(q);
        tJJ = prob2scaledprob(p);
        tJB = prob2scaledprob(1.0f - p);


        //LOG_MSG("%d len",len);
        matrix = m->B_matrix;
        NBECJ = m->B_NBECJ;

        NBECJ[len][J_STATE] = -INFINITY;
        NBECJ[len][B_STATE] = -INFINITY;
        NBECJ[len][N_STATE] = -INFINITY;

        NBECJ[len][C_STATE] = tCT;
        NBECJ[len][E_STATE] = NBECJ[len][C_STATE] +   tEC;

        for(j = 0; j < fhmm->K;j++){
                matrix[len][j] = NBECJ[len][E_STATE] + tXE + fhmm->e[j][a[len-1]];
                //LOG_MSG("EMIT: %d (%d)  %f", a[len-1],len,matrix[len][j]);
        }

        //fprintf(stdout,"\n");
        for(i = len-1; i >= 1; i-- ){
                //LOG_MSG("Looking at: %d %d", i, a[i-1]);
                NBECJ[i][B_STATE] = -INFINITY;
                for(j = 0; j < fhmm->K;j++){
                        NBECJ[i][B_STATE] = logsum(NBECJ[i][B_STATE], matrix[i+1][j] + tBX);
                }
                /* not sure about this ...  */
                NBECJ[i][J_STATE] = logsum(NBECJ[i+1][J_STATE] + tJJ, NBECJ[i][B_STATE] + tJB);
                /* C state */
                NBECJ[i][C_STATE] = NBECJ[i+1][C_STATE] + tCC;
                /* E state */
                NBECJ[i][E_STATE] = logsum(NBECJ[i][J_STATE] + tEJ, NBECJ[i][C_STATE] + tEC);
                /* N state  */
                NBECJ[i][N_STATE] = logsum(NBECJ[i+1][N_STATE] + tNN, NBECJ[i][B_STATE] + tNB);

                for(j = 0; j < fhmm->K;j++){
                        matrix[i][j] = NBECJ[i][E_STATE] + tXE;
                }
                for(j = 0; j < fhmm->K;j++){

                        for(c = 1; c < fhmm->tindex[j][0];c++){
                                f = fhmm->tindex[j][c];
                                matrix[i][j] = logsum(matrix[i][j],fhmm->t[j][f] + matrix[i+1][f]);
                        }
                }
                for(j = 0; j < fhmm->K;j++){
                        //matrix[i][j] = logsum(matrix[i][j], NBECJ[i][E_STATE] + fhmm->tXE);
                        matrix[i][j] += fhmm->e[j][(int)a[i-1]];
                        //LOG_MSG("EMIT: %d (%d)  %f", a[i-1],i,matrix[i][j]);
                }
        }
        NBECJ[0][B_STATE] = -INFINITY;
        for(j = 0; j < fhmm->K;j++){
                NBECJ[0][B_STATE] = logsum(NBECJ[0][B_STATE], matrix[1][j] + tBX);
                //LOG_MSG("B_STATE: %f %f %f", NBECJ[0][B_STATE], fhmm->tBX, fhmm->tNB);
        }
        NBECJ[0][J_STATE] = -INFINITY;
        NBECJ[0][C_STATE] = -INFINITY;
        NBECJ[0][E_STATE] = -INFINITY;

        NBECJ[0][N_STATE] = logsum(NBECJ[1][N_STATE] + tNN, NBECJ[0][B_STATE] + tNB);
        /*for(j = 0; j < fhmm->K;j++){
                matrix[0][j] = -INFINITY;
                }*/

        *ret_score = NBECJ[0][N_STATE] + tSN;
        return OK;
ERROR:
        return FAIL;
}


int posterior_decoding(struct fhmm* fhmm,double** Fmatrix, double** Bmatrix,double score,uint8_t* a, int len,int* path)
{
        int i,j,c,f,best;
        int state;

        double* this_F = 0;
        double* this_B = 0;

        ASSERT(fhmm != NULL, "No model");
        ASSERT(a != NULL, "No sequence");
        ASSERT(len > 0, "Seq is of length 0");


        //const double* trans = 0;
        double max = prob2scaledprob(0.0);
        double total = score;

        for(i = 1; i <= len;i++){
                //last_F = Fmatrix[i-1];
                this_F = Fmatrix[i];
                this_B = Bmatrix[i];
                for(j = 0; j < fhmm->K;j++){
                        this_F[j] = scaledprob2prob((this_F[j]  +( this_B[j] - fhmm->e[j][a[i-1]] )) - total);
                }
        }
        /* need to do separately because sequence[len] is undefined */
        i = len+1;

        this_F = Fmatrix[i];
        this_B = Bmatrix[i];
        for(j = 0; j < fhmm->K;j++){
                this_F[j] = scaledprob2prob((this_F[j]  +( this_B[j] )) - total);
        }

        for(j = 0; j < fhmm->K ;j++){
                Fmatrix[len][j] = Fmatrix[len][j]  + Fmatrix[len+1][END_STATE] ;
                Bmatrix[len][j] = END_STATE;
        }
        best = -1;
        /* Fill B with best transition pointers... */
        for(i = len-1; i >= 0; i-- ){
                for(j = 0; j < fhmm->K;j++){
                        //trans = hmm->transitions[j];
                        max = -INFINITY;
                        for(c = 1; c < fhmm->tindex[j][0];c++){
                                f = fhmm->tindex[j][c];
                                if( Fmatrix[i+1][f] > max){
                                        max  = Fmatrix[i+1][f] ;
                                        best = f;
                                }
                        }
                        Fmatrix[i][j] = max + Fmatrix[i][j];
                        Bmatrix[i][j] = best;
                }
        }
        state = START_STATE;


        /* traceback */
        i = 0;
        while (i <= len){
                //going to
                if(i){
                        path[i-1] = state;
                }

                state =  Bmatrix[i][state];
                i++;
        }


        /*for(i =0; i < len;i++){
                fprintf(stdout,"%2d %2d %2d ", i,(int) a[i], path[i] );
                for(j = 0; j < fhmm->L;j++){
                        fprintf(stdout,"%2f ",scaledprob2prob(fhmm->e[path[i]][j]));;
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");*/
        return OK;
ERROR:
        return FAIL;
}

/* reads finite model from hdf5 file  */
struct fhmm* init_fhmm(char* filename)
{
        struct fhmm* fhmm = NULL;
        //char model_name = NULL;
        ASSERT(filename!= NULL, "No filename");
        /* allocate finite hmm */
        RUNP(fhmm = alloc_fhmm());

        /* get HMM parameters  */
        //RUN(read_fhmm_parameters(fhmm,filename, model_name));

        /* alloc dyn matrices (now that I know how many states there are) */
        //RUN(alloc_dyn_matrices(fhmm));

        /* convert probs into log space/ set tindex to allow for fast-ish dyn
         * programming in case there is a sparse transition matrix */
        //RUN(setup_model(fhmm));

        return fhmm;
ERROR:
        free_fhmm(fhmm);
        return NULL;
}

int setup_model(struct fhmm* fhmm)
{
        int i,j,c;

        ASSERT(fhmm != NULL, "no model");
        /* I need to allocate tindex & convert probs to log space */
        init_logsum();

        /* calculate relative entropy  */
        fhmm->H = 0.0;
        for(i = 0; i < fhmm->K;i++){
                fhmm->H += rel_entropy(fhmm->e[i], fhmm->background, fhmm->L);
        }
        fhmm->H /= (double) fhmm->K;
        fhmm->lambda =  0.69314718055994529 + 1.44 / ((double) fhmm->K * fhmm->H);

        //RUNP(fhmm->tindex = galloc(fhmm->tindex, fhmm->K , fhmm->K+1, 0));
        RUN(galloc(&fhmm->tindex, fhmm->alloc_K , fhmm->alloc_K+1));

        for(i = 0; i < fhmm->K;i++){
                for(j = 0; j < fhmm->K+1;j++){
                        fhmm->tindex[i][j] = 0;
                }
        }

        for(i = 0; i < fhmm->K;i++){
                for(j = 0 ; j < fhmm->L;j++){
                        //LOG_MSG("%d: %f %f   %f",j, fhmm->e[i][j], fhmm->background[j],prob2scaledprob(fhmm->e[i][j] / fhmm->background[j]));
                        fhmm->e[i][j] = prob2scaledprob(fhmm->e[i][j] / fhmm->background[j]);
                }
        }

        for(i = 0; i < fhmm->K;i++){
                for(j = 0; j < fhmm->K;j++){
                        fhmm->t[i][j] = prob2scaledprob(fhmm->t[i][j]);
                }
        }

        for(i = 0; i < fhmm->K;i++){
                c = 0;
                for(j = 0; j < fhmm->K;j++){
                        if(fhmm->t[i][j]  != -INFINITY){
                                //fprintf(stdout,"%d-> %d\n",i, j);
                                fhmm->tindex[i][c+1] = j;
                                c++;
                        }
                }
                fhmm->tindex[i][0] = c+1;
        }
        //exit(0);
        /* background */
        for(i = 0; i < fhmm->L;i++){
                fhmm->background[i] = prob2scaledprob(fhmm->background[i]);
        }
        return OK;
ERROR:
        return FAIL;
}

float rel_entropy(float* vec, float* back,int N)
{
        float kl = 0.0f;
        int i;


        for(i = 0; i < N; i++)
                if (vec[i] > 0.0f) {
                        if (back[i] == 0.0f){
                                kl = INFINITY;

                                return kl;
                        }else{

                                kl += vec[i] * log(vec[i]/back[i]);
                        }
                }
        //return(1.44269504 * kl); /* converts to bits */

        return (kl * 1.44269504f);   /* log2(e)  */
}

int remove_state_for_ploting(struct fhmm*fhmm, int state)
{
        int i,j;
        int a,b;
        ASSERT(fhmm != NULL, "No model");

        float** tmp_trans = NULL;
        float** tmp_emit = NULL;

        //RUNP(tmp_trans = galloc(tmp_trans,fhmm->K-1, fhmm->K-1, 0.0));
        //RUNP(tmp_emit = galloc(tmp_emit,fhmm->K-1, fhmm->L, 0.0));

        RUN(galloc(&tmp_trans,fhmm->K-1, fhmm->K-1));
        RUN(galloc(&tmp_emit,fhmm->K-1, fhmm->L));

        for(i = 0; i < fhmm->K-1;i++){
                for(j = 0; j < fhmm->K-1;j++){
                        tmp_trans[i][j] = 0.0;
                }
        }
        for(i = 0; i < fhmm->K-1;i++){
                for(j = 0; j < fhmm->L;j++){
                        tmp_emit[i][j] = 0.0;
                }
        }

        a = 0;
        for(i = 0; i < fhmm->K;i++){
                if(i != state){
                        b = 0;
                        for(j = 0; j < fhmm->K;j++){
                                if(j != state){
                                        tmp_trans[a][b] = fhmm->t[i][j];
                                        b++;
                                }
                        }
                        a++;
                }
        }
        a = 0;
        for(i = 0; i < fhmm->K;i++){
                if(i != state){
                        for(j = 0; j < fhmm->L;j++){
                                tmp_emit[a][j] = fhmm->e[i][j];
                        }
                        a++;
                }

        }
        gfree(fhmm->e);
        gfree(fhmm->t);
        fhmm->e = tmp_emit;
        fhmm->t= tmp_trans;
        fhmm->K = fhmm->K -1;

        return OK;
ERROR:
        return FAIL;
}

int convert_fhmm_scaled_to_prob(struct fhmm* fhmm)
{
        int i,j;
        ASSERT(fhmm != NULL, "no model");
        /* I need to allocate tindex & convert probs to log space */
        init_logsum();
        //LOG_MSG("states:%d L:%d", fhmm->K, fhmm->L);
        for(i = 0; i < fhmm->K;i++){
                for(j = 0 ; j < fhmm->L;j++){
                        //LOG_MSG("Emission %d %d %f",i,j,fhmm->e[i][j]);
                        fhmm->e[i][j] = scaledprob2prob(fhmm->e[i][j]);
                }
        }

        for(i = 0; i < fhmm->K;i++){
                for(j = 0; j < fhmm->K;j++){
                        //LOG_MSG("Trans %d %d %f",i,j,fhmm->t[i][j]);
                        fhmm->t[i][j] = scaledprob2prob(fhmm->t[i][j]);
                }
        }

        for(i = 0; i < fhmm->L;i++){
                fhmm->background[i] = scaledprob2prob(fhmm->background[i]);
        }
        return OK;
ERROR:
        return FAIL;

}

int convert_fhmm_log_to_prob_for_sampling(struct fhmm* fhmm)
{
        int i,j;
        double sum;
        ASSERT(fhmm != NULL, "no model");
        /* I need to allocate tindex & convert probs to log space */
        init_logsum();


        sum = 0.0;

        for(i = 0; i < fhmm->L;i++){
                fhmm->background[i] = scaledprob2prob(fhmm->background[i]);
                sum += fhmm->background[i];
        }
        for(i = 0; i < fhmm->L;i++){
                fhmm->background[i] /= sum;
        }


        for(i = 1; i < fhmm->L;i++){
                fhmm->background[i] += fhmm->background[i-1];
        }

        /* Step one convert log probabilities */
        for(i = 0;i < fhmm->K;i++){
                sum = 0.0;
                //fprintf(stdout,"K:%d\t",i);
                for(j = 0; j < fhmm->K;j++){
                        //fprintf(stdout,"%f ",fhmm->t[i][j]);
                        fhmm->t[i][j] = scaledprob2prob(fhmm->t[i][j]);
                        sum += fhmm->t[i][j];
                }
                //fprintf(stdout,"\n");
                for(j = 0; j < fhmm->K;j++){
                        fhmm->t[i][j] /= sum;
                }

                for(j = 1;j < fhmm->K;j++){
                        fhmm->t[i][j] += fhmm->t[i][j-1];
                }
                //fprintf(stdout,"K:%d\t",i);
                //for(j = 0; j < fhmm->K;j++){
                //        fprintf(stdout,"%f ",fhmm->t[i][j]);
                //}
                //fprintf(stdout,"\n");

        }
        //fprintf(stdout,"\n");
        /* now emission probabilities  */
        for(i = 0; i < fhmm->K;i++){
                sum = 0.0;
                for(j = 0 ; j < fhmm->L;j++){
                        fhmm->e[i][j] = scaledprob2prob(fhmm->e[i][j]);
                        sum += fhmm->e[i][j];
                }
                //ASSERT(sum != 0.0 , "Sum cannot be zero!");
                for(j = 0 ; j < fhmm->L;j++){
                        fhmm->e[i][j] /= sum;
                }
                for(j = 1; j < fhmm->L;j++){
                        fhmm->e[i][j] += fhmm->e[i][j-1];
                }
                /*fprintf(stdout,"K:%d\t",i);
                for(j = 0; j < fhmm->L;j++){
                        fprintf(stdout,"%f ",fhmm->e[i][j]);
                }
                fprintf(stdout,"\n");*/
        }
        return OK;
ERROR:
        return FAIL;
}


/*int configure_target_len(struct fhmm* fhmm,int len,  int multihit)
{
        ASSERT(fhmm != NULL, "No model");


        float p,q;

        if(multihit){
                q = 0.5f;
                p = (float) len / ((float)len + 3.0f);
        }else{
                q = 0.0f;
                p = (float) len / ((float)len + 2.0f);

        }

        fhmm->tSN = prob2scaledprob(1.0f);
        fhmm->tNN = prob2scaledprob(p);
        fhmm->tNB = prob2scaledprob(1.0f - p);

        fhmm->tBX = prob2scaledprob(2.0f / (float) (fhmm->K * ( fhmm->K + 1.0f)));
        fhmm->tXE = prob2scaledprob(1.0f);
        //LOG_MSG("%f", scaledprob2prob(fhmm->tBX));
        fhmm->tEC = prob2scaledprob(1.0f - q);
        fhmm->tCC = prob2scaledprob(p);
        fhmm->tCT = prob2scaledprob(1.0f - p);

        fhmm->tEJ = prob2scaledprob(q);
        fhmm->tJJ = prob2scaledprob(p);
        fhmm->tJB = prob2scaledprob(1.0f - p);

        fhmm->config_len = len;
        return OK;
ERROR:
        return FAIL;
        }*/
