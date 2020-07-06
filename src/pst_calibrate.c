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
        int min_len;
        int max_len;
        double len;
        double s0;
        double s1;
} wscore;

struct scorelen{
        double s;
        double l;
};

struct slen_store{
        struct scorelen** sl;
        int n_alloc;
        int n;
};

static int alloc_sl(struct slen_store** slen_store);
static int resize_sl(struct slen_store* sls);
static void free_sl(struct slen_store* sls);
static int sort_sl(struct slen_store* sls);

static int sort_sl_by_len(const void *a, const void *b);


KHASH_MAP_INIT_INT(whash, wscore)

//static int score_all(struct pst* p, char* filename, double** sa, int** la, int* n);
static int score_all(struct pst* p, char* filename , struct slen_store** sls);
static int merge_bins(struct window_score** arr, int num,int bin_size);

static int parcel_out_linear_regression(struct window_score** arr, int num, struct slen_store* sls);
static int average_scores_in_bins(struct window_score** arr, int num);
static int linear_regression(struct scorelen** arr, int s, int e, double* a_var,double* b_var);
//static int linear_regression(struct window_score** arr, int num, double* a_var,double* b_var);
static int return_y(double a,double b, double x, double* y);

static int standard_deviation_from_fitted(double* score_arr, int* len_arr, int n_score, double a, double b, double* var);

static int sort_wscore_by_len(const void *a, const void *b);

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
        struct slen_store* sls = NULL;
        struct scorelen* sl = NULL;
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
        RUN(score_all(p, filename,&sls));//   &score_arr, &len_arr, &n_score));

        RUN(sort_sl(sls));

        for(i = 0; i < sls->n;i++){
                sl = sls->sl[i];
                k = kh_put(whash, hash, sl->l, &ret);
                if (!ret){
                        kh_value(hash, k).s0++;
                        kh_value(hash, k).s1 += sl->s;
                }else{
                        kh_value(hash, k).len = sl->l;
                        kh_value(hash, k).min_len = sl->l;
                        kh_value(hash, k).max_len = sl->l;
                        kh_value(hash, k).s0 = 1.0;
                        kh_value(hash, k).s1 = sl->s;
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
        RUN(merge_bins(ws_arr, num_hash_items,1000));

        /* average scores  */
        //RUN(average_scores_in_bins(ws_arr, num_hash_items));
        parcel_out_linear_regression(ws_arr, num_hash_items, sls);
        exit(0);
        /* linear regression  */
        //RUN(linear_regression(ws_arr, num_hash_items, &p->a, &p->b));
        MFREE(ws_arr);

        /* stdev from fitted curve  */
        /* calculate standard deviation of scores from fitted curve  */
        RUN(standard_deviation_from_fitted(score_arr, len_arr, n_score, p->a, p->b, &p->var));

        /* start outlier detection */
        //exit(0);
        outliers = 1;
        while (outliers){
                kh_clear(whash,hash);

                ws_arr = NULL;
                outliers = 0;
                for(i = 0; i < n_score;i++){
                        if(len_arr[i] > 0){
                                RUN(return_y(p->a, p->b, len_arr[i] , &mean));
                                //mean = score_arr[i] - mean;
                                z_score = (score_arr[i]- mean) / p->var;
                                if(fabs(z_score) >= threshold){
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
                LOG_MSG("Iter: %d outliers: %d", iter, outliers);
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
                RUN(merge_bins(ws_arr, num_hash_items,1000));

                /* average scores  */
                RUN(average_scores_in_bins(ws_arr, num_hash_items));

                /* linear regression  */
                //RUN(linear_regression(ws_arr, num_hash_items, &p->a, &p->b));

                RUN(standard_deviation_from_fitted(score_arr, len_arr, n_score, p->a, p->b, &p->var));
                MFREE(ws_arr);

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

        kh_destroy(whash,hash);
        MFREE(score_arr);
        MFREE(len_arr);

        return OK;
ERROR:
        return FAIL;
}


int score_all(struct pst* p, char* filename , struct slen_store** sls)
{
        struct slen_store* sl_store = NULL;
        struct file_handler* f = NULL;
        struct tl_seq_buffer* sb = NULL;


        struct rng_state* rng = NULL;
        struct alphabet* a = NULL;

        float P_M, P_R;

        int chunk,i,len;

        RUN(alloc_sl(&sl_store));


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
                                RUN(create_alphabet(&a, rng, TLALPHABET_NOAMBIGUOUS_DNA));
                        }else if(sb->L == TL_SEQ_BUFFER_PROTEIN){
                                RUN(create_alphabet(&a, rng, TLALPHABET_NOAMBIGIOUS_PROTEIN ));
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
                        sl_store->sl[sl_store->n]->l = (double) len;
                        sl_store->sl[sl_store->n]->s = P_M;
                        sl_store->n++;
                        if(sl_store->n == sl_store->n_alloc){
                                RUN(resize_sl(sl_store));
                        }
                }
                chunk++;
        }

        RUN(close_seq_file(&f));

        free_rng(rng);
        free_tl_seq_buffer(sb);
        free_alphabet(a);
        *sls = sl_store;

        return OK;
ERROR:
        return FAIL;
}


int return_y(double a,double b, double x, double* y)
{

        *y =  a + b * x;
        return OK;
}

int merge_bins(struct window_score** arr, int num,int bin_size)
{
        int i,j,c;
        double cur_count;
        double cur_score;
        double cur_len;
        double n;

        int cur_min;
        int cur_max;

        ASSERT(arr != NULL,"No arr");
        ASSERT(num > 1, "No items");


        qsort(arr , num,  sizeof( wscore*),sort_wscore_by_len);
        cur_max = arr[0]->s0;
        for(i = 1; i < num ;i++){
                cur_max = MACRO_MAX(cur_max, arr[i]->s0);
        }
        bin_size = MACRO_MAX(bin_size, cur_max*2);
        //LOG_MSG("%d %d", cur_max, bin_size);
        //exit(0);
        for(i = 0; i < num;i++){
                cur_len = arr[i]->len;
                cur_min = arr[i]->min_len;
                cur_max = arr[i]->max_len;
                cur_count = arr[i]->s0;
                cur_score = arr[i]->s1;
                n = 1;
                j = i+1;
                while(cur_count < bin_size  && j != num){
                        cur_len += arr[j]->len;
                        cur_min = MACRO_MIN(cur_min, arr[j]->min_len);
                        cur_max = MACRO_MAX(cur_max, arr[j]->max_len);
                        cur_count += arr[j]->s0;
                        cur_score += arr[j]->s1;
                        arr[j]->len = 0.0;
                        arr[j]->s0 = 0.0;
                        arr[j]->s1 = 0.0;
                        n = n + 1.0;
                        j++;
                }
                arr[i]->len = cur_len / n;
                arr[i]->min_len = cur_min;
                arr[i]->max_len = cur_max;
                arr[i]->s0 = cur_count;
                arr[i]->s1 = cur_score;
                i = j-1;

        }
        /* check last bucket  */
        /* The last bin *may* have fewer than 500 sequences - but I think I don't care. */
        c = -1;
        for(i = num-1; i >= 0;i--){
                if(arr[i]->len){
                        c = i;
                        break;
                }
        }
        //LOG_MSG("C:%d",c);
        if(arr[c]->s0 < bin_size){
                cur_len = arr[c]->len;
                cur_min = arr[c]->min_len;
                cur_max = arr[c]->max_len;
                cur_count = arr[c]->s0;
                cur_score = arr[c]->s1;

                n = 1;
                j = c-1;
                while(cur_count < bin_size && j != 0){
                        if(arr[j]->len){

                                cur_len += arr[j]->len;
                                cur_min = MACRO_MIN(cur_min, arr[j]->min_len);
                                cur_max = MACRO_MAX(cur_max, arr[j]->max_len);
                                cur_count += arr[j]->s0;
                                cur_score += arr[j]->s1;
                                arr[j]->len = 0.0;
                                arr[j]->s0 = 0.0;
                                arr[j]->s1 = 0.0;
                                n = n + 1.0;
                        }
                        j--;

                }
                arr[c]->len = cur_len / n;
                arr[c]->min_len = cur_min;
                arr[c]->max_len = cur_max;
                arr[c]->s0 = cur_count;
                arr[c]->s1 = cur_score;
        }
        qsort(arr , num,  sizeof( wscore*),sort_wscore_by_len);

        /* The last bin *may* have fewer than 500 sequences - but I think I don't care. */
        /* for(i = 0; i < num ;i++){ */
        /*         fprintf(stdout,"%d %d %d:  %f %f %f\n",i,arr[i]->min_len,arr[i]->max_len,   arr[i]->len,arr[i]->s0, arr[i]->s1); */
        /* } */

        return OK;
ERROR:
        return FAIL;
}

int parcel_out_linear_regression(struct window_score** arr, int num, struct slen_store* sls)
{
        double s0 = 0.0;
        double s1 = 0.0;
        double s2 = 0.0;
        double threshold = 10.0;
        double mean;
        double v;
        double z_score;
        double a,b;
        int start_sl = 0;
        int stop_sl = 0;

        int i,j,c,index;
        int sanity_score_count;
        /* forward across empty buckets
           these come from the last qsort in the merge small bins function
        */
        for(c = 0; c < num;c++){
                if(arr[c]->len){
                        break;
                }
        }

        LOG_MSG("C: %d - %d",c,num);
        FILE* f_ptr = NULL;
        index= 0;
        f_ptr = fopen("piece.csv", "w");
        for(i = c;i < num;i++){
                sanity_score_count = 0;
                for(j = start_sl; j < sls->n;j++){
                        if(sls->sl[j]->l > arr[i]->max_len){

                                break;
                        }else{
                                sanity_score_count++;
                        }
                }

                //LOG_MSG("%d %f", sanity_score_count, arr[i]->s0);
                //LOG_MSG("%d %d", start_sl , j);
                //LOG_MSG("LEN: %d %d", arr[i]->min_len, arr[i]->max_len);
                stop_sl = j;
                int outliers = 1;

                while(outliers){
                        outliers = 0;
                        RUN(linear_regression(sls->sl, start_sl, stop_sl, &a, &b));
                        /* calc var  */
                        s0 = 0.0;
                        s1 = 0.0;
                        s2 = 0.0;
                        for(j = start_sl;j < stop_sl;j++){
                                ///fprintf(f_ptr,"%f,%f,%f\n", sls->sl[j]->l,sls->sl[j]->s, a + b * sls->sl[j]->l);
                                if(sls->sl[j]->l > 0.0){
                                        mean = sls->sl[j]->s - ( a + b * sls->sl[j]->l);
                                        s0 += 1.0;
                                        s1 += mean;
                                        s2 += mean * mean;
                                }

                        }
                        v =  sqrt ( (s0 * s2 -  pow(s1, 2.0)) /  (s0 * ( s0 - 1.0)));
                        for(j = start_sl;j < stop_sl;j++){
                                if(sls->sl[j]->l > 0.0){
                                        z_score = (sls->sl[j]->s - ( a + b * sls->sl[j]->l)) / v;
                                        if(fabs(z_score) >= threshold){
                                                outliers++;
                                                sls->sl[j]->l = sls->sl[j]->l * -1.0;
                                        }
                                }
                        }
                }
                LOG_MSG("%d a:%f\tb:%f\tv:%f",index,a,b,v);
                index++;
                for(j = start_sl;j < stop_sl;j++){
                        if(sls->sl[j]->l < 0.0){
                                sls->sl[j]->l *= -1.0;
                        }

                        z_score = (sls->sl[j]->s -  ( a + b * sls->sl[j]->l)) / v;
                        fprintf(f_ptr,"%f,%f,%f,%f\n", sls->sl[j]->l,sls->sl[j]->s, a + b * sls->sl[j]->l, z_score);
                }
                start_sl = stop_sl;


        }
        fclose(f_ptr);
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

int linear_regression(struct scorelen** arr, int s, int e, double* a_var,double* b_var)
{
        struct scorelen* sl = NULL;
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
        ASSERT(s < e, "No items");

        for(i = s; i < e;i++){
                sl = arr[i];
                if(sl->l > 0.0){
                        //LOG_MSG("%f %f", sl->l,sl->s);
                        //if(arr[i]->len > 0){
                        sumX = sumX + sl->l;
                        sumX2 = sumX2 + sl->l * sl->l;
                        sumY = sumY +  sl->s;
                        sumXY = sumXY + sl->l* sl->s;
                }
        }
        n = (double) (e -s);
        b = (n*sumXY-sumX*sumY)/(n*sumX2-sumX*sumX);
        a = (sumY - b*sumX)/n;
        *a_var = a;
        *b_var = b;
        //printf("Values are: a=%0.2f and b = %0.2f\n",a,b);
        //printf("\nEquation of best fit is: y = %0.2f + %0.2fx\n",a,b);


        //LOG_MSG("my_fit <- function(x) { y <- %f + %f * x;return (y)}",a,b);

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

        /*double* len_s0 = NULL;
        double* len_s1 = NULL;
        double* len_s2 = NULL;

        int len = 1001;
        int index;

        MMALLOC(len_s0, sizeof(double) * len);
        MMALLOC(len_s1, sizeof(double) * len);
        MMALLOC(len_s2, sizeof(double) * len);
        for(i = 0; i < 1001;i++){

                len_s0[i] = 0.0;
                len_s1[i] = 0.0;
                len_s2[i] = 0.0;

                }*/
        for(i = 0; i < n_score;i++){
                if(len_arr[i] > 0){
                        //index = len_arr[i] / 10;
                        //index = MACRO_MIN(index, 1000);

                        mean = score_arr[i] - ( a + b * (double)len_arr[i]);
                        s0 += 1.0;
                        s1 += mean;
                        s2 += mean * mean;
                        //len_s0[index] += 1.0;
                        //len_s1[index] += mean;
                        //len_s2[index] += mean* mean;
                }
        }
        /*for(i = 0; i < 1001;i++){
                v =  sqrt ( (len_s0[i] * len_s2[i] -  pow(len_s1[i], 2.0)) /  (len_s0[i] * ( len_s0[i] - 1.0)));
                LOG_MSG("%d %f\t%f %f %f",i,v,len_s0[i],len_s1[i],len_s2[i]);
                }*/
        v =  sqrt ( (s0 * s2 -  pow(s1, 2.0)) /  (s0 * ( s0 - 1.0)));
        //LOG_MSG("mean %f",v);
        *var = v;
        //exit(0);
        return OK;
ERROR:
        return FAIL;
}


int sort_wscore_by_len(const void *a, const void *b)
{












        wscore* const *one = a;
        wscore* const *two = b;

        if((*one)->len > (*two)->len){
                return 1;
        }else{
                return -1;
        }
}


int alloc_sl(struct slen_store** slen_store)
{
        struct slen_store* sls = NULL;
        int i;
        MMALLOC(sls, sizeof(struct slen_store));
        sls->n_alloc = 100000;
        sls->n = 0;
        sls->sl = NULL;
        MMALLOC(sls->sl, sizeof(struct scorelen*) * sls->n_alloc);
        for(i = 0; i < sls->n_alloc;i++){
                sls->sl[i] = NULL;
                MMALLOC(sls->sl[i],sizeof(struct scorelen));
        }
        *slen_store = sls;
        return OK;
ERROR:
        return FAIL;
}

int resize_sl(struct slen_store* sls)
{
        int old,i;

        old = sls->n_alloc;
        sls->n_alloc = sls->n_alloc + sls->n_alloc /2;

        MREALLOC(sls->sl, sizeof(struct scorelen*) * sls->n_alloc);
        for(i = old; i < sls->n_alloc;i++){
                sls->sl[i] = NULL;
                MMALLOC(sls->sl[i],sizeof(struct scorelen));
        }
        return OK;
ERROR:
        return FAIL;
}

void free_sl(struct slen_store* sls)
{
        int i;
        if(sls){
                for(i = 0; i < sls->n_alloc;i++){
                        MFREE(sls->sl[i]);
                }
                MFREE(sls->sl);
                MFREE(sls);
        }
}

int sort_sl(struct slen_store* sls)
{
        ASSERT(sls != NULL, "No sls");
        ASSERT(sls->n != 0, "No entries");
        qsort(sls->sl, sls->n, sizeof(struct scorelen*), sort_sl_by_len);
        return OK;
ERROR:
        return FAIL;
}

int sort_sl_by_len(const void *a, const void *b)
{
        struct scorelen * const *one = a;
        struct scorelen * const *two = b;

        if((*one)->l > (*two)->l){
                return 1;
        }else{
                return -1;
        }
}
