
#include "finite_hmm.h"

#include "tllogsum.h"


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

int forward(struct fhmm* fhmm , struct fhmm_dyn_mat* m, float* ret_score, uint8_t* a, int len)
{
        int i,j,c,f;

        float** NBECJ = NULL;
        float** matrix = NULL;

        const float* trans = 0;

        float tmp = 0;




        ASSERT(fhmm != NULL, "No model");
        ASSERT(m != NULL, "No dyn programming  matrix");
        ASSERT(a != NULL, "No sequence");
        ASSERT(len > 0, "Seq is of length 0");


        matrix = m->F_matrix;
        NBECJ = m->F_NBECJ;




        for(j = 0; j < fhmm->K;j++){
                matrix[0][j] = -INFINITY;
        }


        NBECJ[0][N_STATE] = prob2scaledprob(1.0f); /* start -> N state has 1.0 prob */
        NBECJ[0][B_STATE] = fhmm->tNB;
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
                        matrix[i][j] = logsum(matrix[i][j], NBECJ[i-1][B_STATE] + fhmm->tBX);
                        matrix[i][j] += fhmm->e[j][a[i-1]];
                        //cur[c] += fhmm->e[c][a[i-1]];
                        NBECJ[i][E_STATE] = logsum(NBECJ[i][E_STATE], matrix[i][j] + fhmm->tXE);
                }
                /* J */
                NBECJ[i][J_STATE] = logsum(NBECJ[i-1][J_STATE] + fhmm->tJJ, NBECJ[i][E_STATE] + fhmm->tEJ);
                /* C */
                NBECJ[i][C_STATE] = logsum(NBECJ[i-1][C_STATE] + fhmm->tCC, NBECJ[i][E_STATE] + fhmm->tEC);
                /* N */
                NBECJ[i][N_STATE] = NBECJ[i-1][N_STATE] + fhmm->tNN;
                /* B */
                NBECJ[i][B_STATE] = logsum(NBECJ[i][N_STATE] + fhmm->tNB, NBECJ[i][J_STATE]+ fhmm->tJB);


                /* fprintf(stdout,"%0.2f\t%0.2f: \t", NBECJ[i][N_STATE],  NBECJ[i][B_STATE]); */
                /* for(j = 0; j < fhmm->K;j++){ */
                /*         fprintf(stdout,"%0.4f\t",matrix[i][j]); */
                /* } */
                /* fprintf(stdout,": %0.2f\t%0.2f\t%0.2f\n", NBECJ[i][E_STATE],  NBECJ[i][J_STATE], NBECJ[i][C_STATE]); */
        }
        /* LOG_MSG("%f %f %f ",NBECJ[len][C_STATE] , fhmm->tCT,NBECJ[len][C_STATE] + fhmm->tCT); */
        *ret_score = (NBECJ[len][C_STATE] + fhmm->tCT);
        return OK;
ERROR:
        return FAIL;
}

int backward(struct fhmm* fhmm,struct fhmm_dyn_mat* m , float* ret_score, uint8_t* a, int len)
{
        int i,j,c,f;
        float** matrix = NULL;
        float** NBECJ = NULL;
        ASSERT(fhmm != NULL, "No model");
        ASSERT(m != NULL, "No dyn programming  matrix");
        ASSERT(a != NULL, "No sequence");

        ASSERT(len > 0, "Seq is of length 0");

        matrix = m->B_matrix;
        NBECJ = m->B_NBECJ;

        NBECJ[len][J_STATE] = -INFINITY;
        NBECJ[len][B_STATE] = -INFINITY;
        NBECJ[len][N_STATE] = -INFINITY;

        NBECJ[len][C_STATE] = fhmm->tCT;
        NBECJ[len][E_STATE] = NBECJ[len][C_STATE] +   fhmm->tEC;

        for(j = 0; j < fhmm->K;j++){
                matrix[len][j] = NBECJ[len][E_STATE] + fhmm->tXE + fhmm->e[j][a[len-1]];
        }

        //fprintf(stdout,"\n");
        for(i = len-1; i >= 1; i-- ){
                //LOG_MSG("Looking at: %d %d", i, a[i-1]);
                NBECJ[i][B_STATE] = -INFINITY;
                for(j = 0; j < fhmm->K;j++){
                        NBECJ[i][B_STATE] = logsum(NBECJ[i][B_STATE], matrix[i+1][j] + fhmm->tBX);
                }
                /* not sure about this ...  */
                NBECJ[i][J_STATE] = logsum(NBECJ[i+1][J_STATE] + fhmm->tJJ, NBECJ[i][B_STATE] + fhmm->tJB);
                /* C state */
                NBECJ[i][C_STATE] = NBECJ[i+1][C_STATE] + fhmm->tCC;
                /* E state */
                NBECJ[i][E_STATE] = logsum(NBECJ[i][J_STATE] + fhmm->tEJ, NBECJ[i][C_STATE] + fhmm->tEC);
                /* N state  */
                NBECJ[i][N_STATE] = logsum(NBECJ[i+1][N_STATE] + fhmm->tNN, NBECJ[i][B_STATE] + fhmm->tNB);

                for(j = 0; j < fhmm->K;j++){
                        matrix[i][j] = NBECJ[i][E_STATE] + fhmm->tXE;
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
                NBECJ[0][B_STATE] = logsum(NBECJ[0][B_STATE], matrix[1][j] + fhmm->tBX);
                //LOG_MSG("B_STATE: %f %f %f", NBECJ[0][B_STATE], fhmm->tBX, fhmm->tNB);
        }
        NBECJ[0][J_STATE] = -INFINITY;
        NBECJ[0][C_STATE] = -INFINITY;
        NBECJ[0][E_STATE] = -INFINITY;

        NBECJ[0][N_STATE] = logsum(NBECJ[1][N_STATE] + fhmm->tNN, NBECJ[0][B_STATE] + fhmm->tNB);
        for(j = 0; j < fhmm->K;j++){
                matrix[0][j] = -INFINITY;
        }

        *ret_score = NBECJ[0][N_STATE];
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

        //RUNP(fhmm->tindex = galloc(fhmm->tindex, fhmm->K , fhmm->K+1, 0));
        RUN(galloc(&fhmm->tindex, fhmm->alloc_K , fhmm->alloc_K));

        for(i = 0; i < fhmm->alloc_K;i++){
                for(j = 0; j < fhmm->alloc_K;j++){
                        fhmm->tindex[i][j] = 0;
                }
        }

        for(i = 0; i < fhmm->K;i++){
                for(j = 0 ; j < fhmm->L;j++){
                        LOG_MSG("%d: %f %f   %f",j, fhmm->e[i][j], fhmm->background[j],prob2scaledprob(fhmm->e[i][j] / fhmm->background[j]));
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

int remove_state_for_ploting(struct fhmm*fhmm, int state)
{
        int i,j;
        int a,b;
        ASSERT(fhmm != NULL, "No model");

        double** tmp_trans = NULL;
        double** tmp_emit = NULL;

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


