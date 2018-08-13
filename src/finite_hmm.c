
#include "finite_hmm.h"

static struct fhmm* alloc_fhmm(void);
static int setup_model(struct fhmm* fhmm);
static int alloc_dyn_matrices(struct fhmm* fhmm);
static int read_hmm_parameters(struct fhmm* fhmm, char* filename);


int random_model_score(struct fhmm* fhmm, uint8_t* a, int len, int expected_len)
{
        int i;

        float score;
        float r;                /* self transition */
        float e;                /* 1- r (exit) */
        float* b;
        ASSERT(fhmm != NULL, "No model");
        ASSERT(a != NULL, "No sequence");
        ASSERT(len > 0, "Seq is of length 0");

        b = fhmm->background;
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


        fhmm->r_score = score;


        return OK;
ERROR:
        return FAIL;

}


int forward(struct fhmm* fhmm, uint8_t* a, int len)
{
        int i,j,c,f;

        float** matrix = NULL;
        float* last= 0;
        float* cur = 0;
        const float* trans = 0;

        float tmp = 0;

        ASSERT(fhmm != NULL, "No model");
        ASSERT(a != NULL, "No sequence");

        ASSERT(len > 0, "Seq is of length 0");

        matrix = fhmm->F_matrix;
        cur = matrix[0];

        for(j = 0; j < fhmm->K;j++){
                cur[j]  = -INFINITY;
        }
        cur[IHMM_START_STATE] = 0.0f;

        for(i = 1; i < len+1;i++){
                last = cur;
                cur = matrix[i];
                for(j = 0; j < fhmm->K;j++){
                        cur[j] = -INFINITY;
                }
                for(j = 0; j < fhmm->K;j++){
                        tmp = last[j];
                        trans = fhmm->t[j];
                        for(c = 1; c <   fhmm->tindex[j][0];c++){
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
        fhmm->f_score = cur[IHMM_END_STATE];// matrix[ENDSTATE][i];
        return OK;
ERROR:
        return FAIL;
}

/* reads finite model from hdf5 file  */
struct fhmm* init_fhmm(char* filename)
{
        struct fhmm* fhmm = NULL;
        ASSERT(filename!= NULL, "No filename");



        RUNP(fhmm = alloc_fhmm());

        /* get HMM parameters  */
        RUN(read_hmm_parameters(fhmm,filename));

        /* alloc dyn matrices (now that I know how many states there are) */
        RUN(alloc_dyn_matrices(fhmm));

        /* convert probs into log space/ set tindex to allow for fast-ish dyn
         * programming in case there is a sparse transition matrix */
        RUN(setup_model(fhmm));

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

        RUNP(fhmm->tindex = malloc_2d_int(fhmm->tindex, fhmm->K , fhmm->K+1, 0));

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






int alloc_dyn_matrices(struct fhmm* fhmm)
{
        ASSERT(fhmm!= NULL, "No model");

        fhmm->alloc_matrix_len = 1024;

        RUNP(fhmm->F_matrix = malloc_2d_float(fhmm->F_matrix, fhmm->alloc_matrix_len, fhmm->K, 0.0));
        RUNP(fhmm->B_matrix = malloc_2d_float(fhmm->B_matrix, fhmm->alloc_matrix_len, fhmm->K, 0.0));
        return OK;
ERROR:
        return FAIL;
}

int realloc_dyn_matrices(struct fhmm* fhmm,int new_len)
{
        ASSERT(fhmm != NULL, "No model");
        ASSERT(new_len > 0, "newlen has to be > 0");
        if(fhmm->alloc_matrix_len < new_len){
                while(fhmm->alloc_matrix_len < new_len){
                        fhmm->alloc_matrix_len = fhmm->alloc_matrix_len << 1;
                }
                RUNP(fhmm->F_matrix = malloc_2d_float(fhmm->F_matrix, fhmm->alloc_matrix_len, fhmm->K, 0.0));
                RUNP(fhmm->B_matrix = malloc_2d_float(fhmm->B_matrix, fhmm->alloc_matrix_len, fhmm->K, 0.0));
        }
        return OK;
ERROR:
        return FAIL;
}

int read_hmm_parameters(struct fhmm* fhmm, char* filename)
{

        struct hdf5_data* hdf5_data = NULL;

        float sum;
        int i,j;

        ASSERT(fhmm!= NULL, "No model");
        ASSERT(filename != NULL, "No filename");
        ASSERT(my_file_exists(filename) != 0,"File %s does not exist.",filename);



        /* read in hdf5 file and get emission and transition matrix */



        hdf5_data = hdf5_create();

        hdf5_open_file(filename,hdf5_data);

        hdf5_read_attributes(hdf5_data,hdf5_data->file);
        ASSERT(hdf5_data->num_attr != 0 , "Could not find attributes");
        print_attributes(hdf5_data);
        get_group_names(hdf5_data);
        fprintf(stdout,"Groups:\n");
        for(i = 0; i < hdf5_data->grp_names->num_names;i++){
                fprintf(stdout,"%d %s\n",i,hdf5_data->grp_names->names[i]);
        }

        hdf5_open_group("imodel",hdf5_data);
        hdf5_read_attributes(hdf5_data,hdf5_data->group);
        ASSERT(hdf5_data->num_attr != 0 , "Could not find attributes");
        print_attributes(hdf5_data);
        for(i = 0; i < hdf5_data->num_attr;i++){
                if(!strncmp("Number of states", hdf5_data->attr[i]->attr_name, 16)){
                        fhmm->K = hdf5_data->attr[i]->int_val;
                }

                if(!strncmp("Number of letters", hdf5_data->attr[i]->attr_name, 17)){
                        fhmm->L = hdf5_data->attr[i]->int_val;
                }

        }
        hdf5_close_group(hdf5_data);

        hdf5_open_group("SequenceInformation",hdf5_data);

        hdf5_read_dataset("background",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 1, "Could not read transition_counts");
        fhmm->background = (float*) hdf5_data->data;
        hdf5_close_group(hdf5_data);


        hdf5_open_group("fmodel",hdf5_data);

        hdf5_read_dataset("emission",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read transition_counts");
        fhmm->e = (float**) hdf5_data->data;


        hdf5_read_dataset("transition",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read transition_counts");
        fhmm->t = (float**) hdf5_data->data;

        hdf5_close_group(hdf5_data);

        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);


        for(i = 0; i < fhmm->K;i++){
                sum = 0.0;
                for(j =0; j < fhmm->L;j++){
                        fprintf(stdout,"%f ",fhmm->e[i][j]);
                        sum += fhmm->e[i][j];
                }
                fprintf(stdout,"sum:%f\n",sum);
        }


        for(i = 0; i < fhmm->K;i++){
                sum = 0.0;
                for(j =0; j < fhmm->K;j++){
                        fprintf(stdout,"%f ",fhmm->t[i][j]);
                        sum += fhmm->t[i][j];
                }
                fprintf(stdout,"sum:%f\n",sum);
        }



        return OK;
ERROR:
        return FAIL;

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
                        free_2d((void**) fhmm->F_matrix);
                }
                if(fhmm->B_matrix){
                        free_2d((void**) fhmm->B_matrix);
                }

                if(fhmm->e){
                        free_2d((void**) fhmm->e);
                }
                if(fhmm->t){
                        free_2d((void**) fhmm->t);
                }
                if(fhmm->background){
                        MFREE(fhmm->background);
                }
                if(fhmm->tindex){
                        free_2d((void**)fhmm->tindex);
                }

                MFREE(fhmm);
        }
}


