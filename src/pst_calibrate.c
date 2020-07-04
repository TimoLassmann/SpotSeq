#include "tldevel.h"
#include "tlrng.h"
#include "tlseqio.h"
#include "tlmisc.h"
#include "tlalphabet.h"

#include "khash.h"

#include "pst_structs.h"
#include "pst.h"

#define PST_CALIBRATE_IMPORT
#include "pst_calibrate.h"

typedef struct window_score{
        double len;
        double s0;
        double s1;
} wscore;

KHASH_MAP_INIT_INT(whash, wscore)

static int score_all(struct pst* p, char* filename, double** sa, int** la, int* n);

static int merge_small_bins(struct window_score** arr, int num);
static int average_scores_in_bins(struct window_score** arr, int num);
static int linear_regression(struct window_score** arr, int num, double* a_var,double* b_var);
static int return_y(double a,double b, double x, double* y);

static int standard_deviation_from_fitted(double* score_arr, int* len_arr, int n_score, double a, double b, double* var);


/* run pst score */
/* record Log Odds scores */
/* summarize per length. (in hash -)  */
/* sort hash entries by length */
/* merge nearby entries if few sequences  */
/* linear regression  */
/* run pst again (?) to calculate standard deviation from fitted curve  */

int calibrate_pst(struct pst* p, char* filename,double threshold)
{
        struct window_score** ws_arr = NULL;
        int* len_arr = NULL;
        double* score_arr = NULL;

        int n_score;

        //double s0,s1,s2;
        double mean;
        double z_score;

        int i;

        FILE* f_ptr = NULL;

        khash_t(whash) *hash = kh_init(whash);
        khiter_t k;
        int num_hash_items;
        int ret;
        int outliers = 100;

        int iter = 1;

        /* score everything  */
        RUN(score_all(p, filename, &score_arr, &len_arr, &n_score));

        for(i = 0; i < n_score;i++){
                k = kh_put(whash, hash, len_arr[i], &ret);
                if (!ret){
                        kh_value(hash, k).s0++;
                        kh_value(hash, k).s1 += score_arr[i];
                }else{
                        kh_value(hash, k).len = (double) len_arr[i];
                        kh_value(hash, k).s0 = 1.0;
                        kh_value(hash, k).s1 = score_arr[i];
                }
        }
        /* putting these scores, or better array of structs, makes it easier for my brain.. */
        num_hash_items = kh_size(hash);

        MMALLOC(ws_arr, sizeof(struct window_score*) * num_hash_items);
        i = 0;
        for (k = kh_begin(hash); k != kh_end(hash); ++k){
                if (kh_exist(hash, k)){
                        ws_arr[i] = &kh_value(hash, k);
                        i++;
                }
        }

        /* merge small bins */
        RUN(merge_small_bins(ws_arr, num_hash_items));

        /* average scores  */
        RUN(average_scores_in_bins(ws_arr, num_hash_items));

        /* linear regression  */
        RUN(linear_regression(ws_arr, num_hash_items, &p->a, &p->b));
        MFREE(ws_arr);

        /* stdev from fitted curve  */
        /* calculate standard deviation of scores from fitted curve  */
        RUN(standard_deviation_from_fitted(score_arr, len_arr, n_score, p->a, p->b, &p->var));

        /* start outlier detection */
        //exit(0);
        outliers = 1;
        while (outliers){
                kh_destroy(whash,hash);

                khash_t(whash) *hash = kh_init(whash);
                ws_arr = NULL;
                outliers = 0;
                for(i = 0; i < n_score;i++){
                        if(len_arr[i] > 0){
                                RUN(return_y(p->a, p->b, len_arr[i] , &mean));
                                //mean = score_arr[i] - mean;
                                z_score = (score_arr[i]- mean) / p->var;
                                if(z_score >= threshold){
                                        outliers++;
                                        len_arr[i] = len_arr[i] * -1;
                                }else{

                                        k = kh_put(whash, hash, len_arr[i], &ret);
                                        if (!ret){
                                                kh_value(hash, k).s0++;
                                                kh_value(hash, k).s1 += score_arr[i];
                                        }else{
                                                kh_value(hash, k).len = (double) len_arr[i];
                                                kh_value(hash, k).s0 = 1.0;
                                                kh_value(hash, k).s1 = score_arr[i];
                                        }
                                }
                        }
                }
                num_hash_items = kh_size(hash);

                MMALLOC(ws_arr, sizeof(struct window_score*) * num_hash_items);
                i = 0;
                for (k = kh_begin(hash); k != kh_end(hash); ++k){
                        if (kh_exist(hash, k)){
                                ws_arr[i] = &kh_value(hash, k);
                                i++;
                        }
                }

                /* merge small bins */
                RUN(merge_small_bins(ws_arr, num_hash_items));

                /* average scores  */
                RUN(average_scores_in_bins(ws_arr, num_hash_items));

                /* linear regression  */
                RUN(linear_regression(ws_arr, num_hash_items, &p->a, &p->b));
                MFREE(ws_arr);
                LOG_MSG("Iter: %d outliers: %d", iter, outliers);
                iter++;

        }


        f_ptr = fopen("logscores.txt", "w");
        for(i = 0; i < n_score;i++){
                if(len_arr[i] < 0){
                        len_arr[i] = len_arr[i] * -1;
                }
                RUN(return_y(p->a, p->b, len_arr[i] , &mean));
                //mean = score_arr[i] - mean;
                z_score = (score_arr[i]- mean) / p->var;
                fprintf(f_ptr,"%d,%f,%f\n", len_arr[i] , score_arr[i], z_score);
        }
        fclose(f_ptr);


        MFREE(score_arr);
        MFREE(len_arr);

        return OK;
ERROR:
        return FAIL;
}


int score_all(struct pst* p, char* filename, double** sa, int** la, int* n)
{
        struct file_handler* f = NULL;
        struct tl_seq_buffer* sb = NULL;


        struct rng_state* rng = NULL;
        struct alphabet* a = NULL;

        float P_M, P_R;
        int* len_arr = NULL;
        double* score_arr = NULL;
        int n_score_alloc;
        int n_score;

        int chunk,i,len;
        n_score_alloc = 1000000;

        n_score = 0;
        MMALLOC(len_arr, sizeof(int) * n_score_alloc);
        MMALLOC(score_arr, sizeof(double) * n_score_alloc);

        //LOG_MSG("%s",infile);
        if(!my_file_exists(filename)){
                ERROR_MSG("File %s not found");
        }
        RUNP(rng = init_rng(0));

        RUN(open_fasta_fastq_file(&f, filename, TLSEQIO_READ));


        chunk =1;
        while(1){
                RUN(read_fasta_fastq_file(f, &sb, 100000));
                if(chunk == 1){
                        if(sb->L == TL_SEQ_BUFFER_DNA){
                                RUN(create_alphabet(&a, rng, TLALPHABET_DEFAULT_DNA));
                        }else if(sb->L == TL_SEQ_BUFFER_PROTEIN){
                                RUN(create_alphabet(&a, rng, TLALPHABET_DEFAULT_PROTEIN));
                        }
                }
                if(sb->num_seq == 0){
                        break;
                }
                LOG_MSG("Working on chunk: %d",chunk);
                for(i = 0; i < sb->num_seq;i++){
                        len = sb->sequences[i]->len;
                        RUN(convert_to_internal(a, (uint8_t*)sb->sequences[i]->seq,len));
                        RUN(score_pst(p, sb->sequences[i]->seq, len, &P_M,&P_R));
                        P_M = P_M - P_R;
                        score_arr[n_score] = P_M;
                        len_arr[n_score] = len;
                        n_score++;
                        if(n_score > n_score_alloc){
                                n_score_alloc = n_score_alloc + n_score_alloc / 2;
                                MREALLOC(len_arr, sizeof(int) * n_score_alloc);
                                MREALLOC(score_arr, sizeof(double) * n_score_alloc);
                        }
                }
                chunk++;
        }


        RUN(close_seq_file(&f));




        *sa = score_arr;
        *la = len_arr;
        *n = n_score;
        return OK;
ERROR:
        return FAIL;
}


int return_y(double a,double b, double x, double* y)
{

        *y =  a + b * x;
        return OK;
}

int merge_small_bins(struct window_score** arr, int num)
{
        int i,j;
        double cur_count;
        double cur_score;
        double cur_len;
        double n;
        ASSERT(arr != NULL,"No arr");
        ASSERT(num > 1, "No items");

        for(i = 0; i < num;i++){
                cur_len = arr[i]->len;
                cur_count = arr[i]->s0;
                cur_score = arr[i]->s1;
                n = 1;
                j = i+1;
                while(cur_count < 500.0 && j != num){
                        cur_len += arr[j]->len;
                        cur_count += arr[j]->s0;
                        cur_score += arr[j]->s1;
                        arr[j]->len = 0.0;
                        arr[j]->s0 = 0.0;
                        arr[j]->s1 = 0.0;
                        n = n + 1.0;
                        j++;
                }
                arr[i]->len = cur_len / n;
                arr[i]->s0 = cur_count;
                arr[i]->s1 = cur_score;
                i = j-1;

        }
        /* The last bin *may* have fewer than 500 sequences - but I think I don't care. */
        /*for(i = num-1; i > num-1000;i--){
                if(arr[i]->len){
                fprintf(stdout," %f %f %f\n", arr[i]->len,arr[i]->s0, arr[i]->s1);
                }
                }*/
        return OK;
ERROR:
        return FAIL;
}

int average_scores_in_bins(struct window_score** arr, int num)
{
        int i;
        ASSERT(arr != NULL,"No arr");
        ASSERT(num > 1, "No items");

        for(i = 0; i < num;i++){
                arr[i]->s1 = arr[i]->s1 / arr[i]->s0;
        }

        return OK;
ERROR:
        return FAIL;
}

/* y = a+bx
 */
int linear_regression(struct window_score** arr, int num, double* a_var,double* b_var)
{
        double a;
        double b;

        double sumX;
        double sumX2;
        double sumY;
        double sumXY;
        double n;

        int i;
        sumX = 0.0;
        sumX2 = 0.0;
        sumY = 0.0;
        sumXY = 0.0;
        ASSERT(arr != NULL,"No arr");
        ASSERT(num > 1, "No items");

        for(i = 0; i < num;i++){
                if(arr[i]->len > 0){
                        sumX = sumX + arr[i]->len;
                        sumX2 = sumX2 + arr[i]->len * arr[i]->len;
                        sumY = sumY + arr[i]->s1;
                        sumXY = sumXY + arr[i]->len* arr[i]->s1;
                }
        }
        n = (double) num;
        b = (n*sumXY-sumX*sumY)/(n*sumX2-sumX*sumX);
        a = (sumY - b*sumX)/n;
        *a_var = a;
        *b_var = b;
        printf("Values are: a=%0.2f and b = %0.2f\n",a,b);
        printf("\nEquation of best fit is: y = %0.2f + %0.2fx\n",a,b);


        fprintf(stdout,"my_fit <- function(x) { y <- %f + %f * x;return (y)}\n",a,b);
        return OK;
ERROR:
        return FAIL;
}


int standard_deviation_from_fitted(double* score_arr, int* len_arr, int n_score, double a, double b, double* var)
{

        /* calculate standard deviation of scores from fitted curve  */
        double s0 = 0.0;
        double s1 = 0.0;
        double s2 = 0.0;

        double mean;
        double v;
        int  i;

        for(i = 0; i < n_score;i++){
                if(len_arr[i] > 0){
                        mean = score_arr[i] - ( a + b * (double)len_arr[i]);
                        s0 += 1.0;
                        s1 += mean;
                        s2 += mean * mean;
                }
        }

        v =  sqrt ( (s0 * s2 -  pow(s1, 2.0)) /  (s0 * ( s0 - 1.0)));
        *var = v;
        return OK;
}
