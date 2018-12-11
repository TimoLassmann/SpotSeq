
#include "finite_hmm.h"

int random_model_score(double* b, double* ret_score,  uint8_t* a, int len, int expected_len)
{
        int i;
        double score;
        double r;                /* self transition */
        double e;                /* 1- r (exit) */

        ASSERT(b != NULL, "No background probabilities");
        ASSERT(a != NULL, "No sequence");
        ASSERT(len > 0, "Seq is of length 0");

        /* initialise transitions  */
        r = (double)expected_len / ((double) expected_len + 1.0);
        e = 1.0 -r;
        r = prob2scaledprob(r);
        e = prob2scaledprob(e);

        /* initialise scores */
        score = prob2scaledprob(1.0);

        for(i = 0; i < len; i++){
                score += b[a[i]];
                score += r;
        }
        score += e;

        *ret_score = score;

        return OK;
ERROR:
        return FAIL;
}

int forward(struct fhmm* fhmm,double** matrix, double* ret_score, uint8_t* a, int len)
{
        int i,j,c,f;

        double* last= 0;
        double* cur = 0;
        const double* trans = 0;

        double tmp = 0;

        ASSERT(fhmm != NULL, "No model");
        ASSERT(matrix != NULL, "No dyn programming  matrix");
        ASSERT(a != NULL, "No sequence");

        ASSERT(len > 0, "Seq is of length 0");

        cur = matrix[0];

        for(j = 0; j < fhmm->K;j++){
                cur[j]  = -INFINITY;
        }
        cur[IHMM_START_STATE] = 0.0;

        for(i = 1; i < len+1;i++){
                last = cur;
                cur = matrix[i];
                for(j = 0; j < fhmm->K;j++){
                        //LOG_MSG("writing to state: %d",j);
                        cur[j] = -INFINITY;
                }
                for(j = 0; j < fhmm->K;j++){
                        tmp = last[j];
                        trans = fhmm->t[j];
                        for(c = 1; c < fhmm->tindex[j][0];c++){
                                f = fhmm->tindex[j][c];
                                cur[f] = logsum(cur[f], tmp + trans[f] );//+ hmm->emissions[c][(int)a[i-1]]);
                        }
                }
                for(c = 2;c < fhmm->K;c++){
                        cur[c] += fhmm->e[c][a[i-1]];
                }
        }
        //All goes to 1.
        last = cur;//matrix[len];
        cur = matrix[len+1];

        for(j = 0; j < fhmm->K;j++){
                cur[j] = -INFINITY;// prob2scaledprob(0.0);
        }

        for(j = 2; j < fhmm->K;j++){
                cur[IHMM_END_STATE] = logsum(cur[IHMM_END_STATE],last[j] + fhmm->t[j][IHMM_END_STATE]);
        }
        *ret_score = cur[IHMM_END_STATE];
        //fhmm->f_score = cur[IHMM_END_STATE];// matrix[ENDSTATE][i];
        return OK;
ERROR:
        return FAIL;
}

int backward(struct fhmm* fhmm,double** matrix, double* ret_score, uint8_t* a, int len)
{
        int i,j,c,f;

        double* next= 0;
        double* cur = 0;
        const double* trans = 0;

        ASSERT(fhmm != NULL, "No model");
        ASSERT(matrix != NULL, "No dyn programming  matrix");
        ASSERT(a != NULL, "No sequence");

        ASSERT(len > 0, "Seq is of length 0");

        cur = matrix[len+1];

        for(j = 0; j < fhmm->K;j++){
                cur[j] = -INFINITY;
        }

        cur[IHMM_END_STATE] = 0.0f;

        next = cur;

        cur = matrix[len];
        for(j = 0; j < fhmm->K;j++){
                cur[j] = fhmm->t[j][IHMM_END_STATE] + next[IHMM_END_STATE];
        }
        for(c = 2;c < fhmm->K;c++){
                cur[c] += fhmm->e[c][a[len-1]];
        }

        // backward recursion...
        for(i = len-1; i > 0; i-- ){
                next = cur;
                cur = matrix[i];
                for(j = 0; j < fhmm->K;j++){
                        trans = fhmm->t[j];
                        cur[j] = -INFINITY;
                        for(c = 1; c < fhmm->tindex[j][0];c++){
                                f = fhmm->tindex[j][c];
                                cur[j] = logsum(cur[j],trans[f] + next[f]);// hmm->emissions[c][(int)a[i]]);// + next[c]);
                        }
                }
                for(j = 2; j < fhmm->K;j++){
                        cur[j] += fhmm->e[j][(int)a[i-1]];
                }
        }

        cur = matrix[0];
        next = matrix[1];

        for(j = 0; j < fhmm->K;j++){
                cur[j] = -INFINITY;// prob2scaledprob(0.0f);
        }
        for(i = 0; i < fhmm->K;i++){
                cur[IHMM_START_STATE] = logsum(cur[IHMM_START_STATE], fhmm->t[IHMM_START_STATE][i] + next[i]);//  + hmm->emissions[i][(int)a[0]]  );
        }
        *ret_score = cur[IHMM_START_STATE];
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
                Fmatrix[len][j] = Fmatrix[len][j]  + Fmatrix[len+1][IHMM_END_STATE] ;
                Bmatrix[len][j] = IHMM_END_STATE;
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
        state = IHMM_START_STATE;


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
        RUN(alloc_dyn_matrices(fhmm));

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

        RUNP(fhmm->tindex = galloc(fhmm->tindex, fhmm->K , fhmm->K+1, 0));

        for(i = 0; i < fhmm->K;i++){
                for(j = 0 ; j < fhmm->L;j++){
                        fhmm->e[i][j] = prob2scaledprob(fhmm->e[i][j]);
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
                                fhmm->tindex[i][c+1] = j;
                                c++;
                        }
                }
                fhmm->tindex[i][0] = c+1;
        }

        /* background */
        for(i = 0; i < fhmm->L;i++){
                fhmm->background[i] = prob2scaledprob(fhmm->background[i]);
        }
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
        LOG_MSG("states:%d L:%d", fhmm->K, fhmm->L);
        for(i = 0; i < fhmm->K;i++){
                for(j = 0 ; j < fhmm->L;j++){
                        fhmm->e[i][j] = scaledprob2prob(fhmm->e[i][j]);
                }
        }

        for(i = 0; i < fhmm->K;i++){
                for(j = 0; j < fhmm->K;j++){
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

int alloc_dyn_matrices(struct fhmm* fhmm)
{
        ASSERT(fhmm!= NULL, "No model");

        fhmm->alloc_matrix_len = 1024;

        RUNP(fhmm->F_matrix = galloc(fhmm->F_matrix, fhmm->alloc_matrix_len, fhmm->K, 0.0));
        RUNP(fhmm->B_matrix = galloc(fhmm->B_matrix, fhmm->alloc_matrix_len, fhmm->K, 0.0));
        return OK;
ERROR:
        return FAIL;
}

int realloc_dyn_matrices(struct fhmm* fhmm,int new_len)
{
        ASSERT(fhmm != NULL, "No model");
        ASSERT(new_len > 0, "newlen has to be > 0");
        ASSERT(fhmm->alloc_matrix_len > 0, "No matrix allocated yet...");

        if(fhmm->alloc_matrix_len < new_len){
                while(fhmm->alloc_matrix_len < new_len){
                        fhmm->alloc_matrix_len = fhmm->alloc_matrix_len << 1;
                }
                RUNP(fhmm->F_matrix = galloc(fhmm->F_matrix, fhmm->alloc_matrix_len, fhmm->K, 0.0));
                RUNP(fhmm->B_matrix = galloc(fhmm->B_matrix, fhmm->alloc_matrix_len, fhmm->K, 0.0));
        }
        return OK;
ERROR:
        return FAIL;
}

struct fhmm*  read_fhmm_parameters(struct hdf5_data* hdf5_data, char* group)
{
        struct fhmm* fhmm = NULL;
        double sum;
        int i,j;


        ASSERT(hdf5_data != NULL, "No filename");

        RUNP(fhmm = alloc_fhmm());



        hdf5_open_group(group,hdf5_data);

        if((hdf5_data->dataset = H5Dopen(hdf5_data->group, "K",H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        //printf ("H5Dopen returns: %d\n", hdf5_data->dataset);
        hdf5_read_attributes(hdf5_data,hdf5_data->dataset);
        hdf5_data->datatype  = H5Dget_type(hdf5_data->dataset );     /* datatype handle */
        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        hdf5_data->rank      = H5Sget_simple_extent_ndims(hdf5_data->dataspace);
        hdf5_data->status  = H5Sget_simple_extent_dims(hdf5_data->dataspace,hdf5_data->dim , NULL);
        hdf5_data->status = H5Dread(hdf5_data->dataset, hdf5_data->datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fhmm->K);
        if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");

        if((hdf5_data->dataset = H5Dopen(hdf5_data->group, "L",H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        //printf ("H5Dopen returns: %d\n", hdf5_data->dataset);
        hdf5_read_attributes(hdf5_data,hdf5_data->dataset);
        hdf5_data->datatype  = H5Dget_type(hdf5_data->dataset );     /* datatype handle */
        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        hdf5_data->rank      = H5Sget_simple_extent_ndims(hdf5_data->dataspace);
        hdf5_data->status  = H5Sget_simple_extent_dims(hdf5_data->dataspace,hdf5_data->dim , NULL);
        hdf5_data->status = H5Dread(hdf5_data->dataset, hdf5_data->datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fhmm->L);
        if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");


        hdf5_read_dataset("emission",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read emission_counts");
        fhmm->e = (double**) hdf5_data->data;

        hdf5_read_dataset("transition",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read transition_counts");
        fhmm->t = (double**) hdf5_data->data;

        hdf5_read_dataset("transition_index",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read transition_index");
        fhmm->tindex = (int**) hdf5_data->data;


        hdf5_read_dataset("background",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 1, "Could not read transition_background");
        fhmm->background = (double*) hdf5_data->data;


        hdf5_close_group(hdf5_data);
        return fhmm;
ERROR:
        return NULL;

}


struct fhmm* alloc_fhmm(void)
{

        struct fhmm* fhmm = NULL;

        MMALLOC(fhmm, sizeof(struct fhmm));
        fhmm->F_matrix = NULL;
        fhmm->B_matrix = NULL;
        fhmm->e = NULL;
        fhmm->t = NULL;
        fhmm->tindex = NULL;
        fhmm->background = NULL;
        fhmm->K = 0;
        fhmm->L = 0;
        fhmm->f_score = 0.0;
        fhmm->b_score = 0.0;

        fhmm->alloc_matrix_len = 0;


        return fhmm;
ERROR:
        free_fhmm(fhmm);
        return NULL;
}

void free_fhmm(struct fhmm* fhmm)
{
        if(fhmm){
                if(fhmm->F_matrix){
                        gfree(fhmm->F_matrix);
                        //free_2d((void**) fhmm->F_matrix);
                }
                if(fhmm->B_matrix){
                        gfree(fhmm->B_matrix);
                        //free_2d((void**) fhmm->B_matrix);
                }

                if(fhmm->e){
                        gfree(fhmm->e);
                        //free_2d((void**) fhmm->e);
                }
                if(fhmm->t){
                        gfree(fhmm->t);
                        //free_2d((void**) fhmm->t);
                }
                if(fhmm->background){
                        MFREE(fhmm->background);
                }
                if(fhmm->tindex){
                        gfree(fhmm->tindex);
                        //free_2d((void**)fhmm->tindex);
                }

                MFREE(fhmm);
        }
}
