
#include <omp.h>
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
        double min_len;
        double max_len;
        double len;
        double s0;
        double s1;
        double s2;
} wscore;

struct scorelen{
        double s;
        double l;
};

struct slen_store{
        struct scorelen** sl;
        int n_alloc;
        int n;
        int offset;
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

static int parcel_out_linear_regression(double*** fit, struct window_score** arr, int num, struct slen_store* sls,double threshold);
static int average_scores_in_bins(struct window_score** arr, int num);
static int linear_regression_bins(struct window_score** arr, int num, double* a_var,double* b_var);
static int linear_regression(struct scorelen** arr, int s, int e, double* a_var,double* b_var);
static int linear_regression_var(struct scorelen** arr, int s, int e,double a, double b, double* a_var,double* b_var);
static int standard_deviation_from_fitted(double* score_arr, int* len_arr, int n_score, double a, double b, double* var);

static int sort_wscore_by_len(const void *a, const void *b);

//static int make_intervals(struct  scorelen** arr, int l, double min_ivlen, double min_items,struct window_score*** ret);
static int make_intervals(struct  scorelen** arr, int l, double min_ivlen, double min_items,struct window_score*** ws_ret, int* ws_num);

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

        int i,j,len,last;

        khash_t(whash) *hash = kh_init(whash);
        khiter_t k;
        int num_items;
        int ret;




        /* score everything  */
        RUN(score_all(p, filename,&sls));//   &score_arr, &len_arr, &n_score));

        RUN(sort_sl(sls));

        /*FILE* f_ptr = NULL;

        RUNP(f_ptr = fopen("pstscores.csv", "w"));
        for(i = 0; i < sls->n;i++){
                fprintf(f_ptr ,"%f\n", sls->sl[i]->s);
        }
        fclose(f_ptr);*/
        //exit(0);
        RUN(make_intervals( sls->sl,sls->n, 2.0,500.0,&ws_arr,&num_items));
        /*double s0 = 0.0;
        for(i = 0; i < num_items;i++){

                LOG_MSG("%d %f %f %f %f ",i, ws_arr[i]->len, ws_arr[i]->min_len,ws_arr[i]->max_len,ws_arr[i]->s0);
                s0 += ws_arr[i]->s0;
        }
        LOG_MSG("%d %f",sls->n,s0);
        */
        //exit(0);
        /*
        for(i = 0; i < sls->n;i++){
                sl = sls->sl[i];
                k = kh_put(whash, hash, sl->l, &ret);
                if (!ret){
                        kh_value(hash, k).s0++;
                        kh_value(hash, k).s1 += sl->s;
                        kh_value(hash, k).s2 += sl->s * sl->s;
                }else{
                        kh_value(hash, k).len = sl->l;
                        kh_value(hash, k).min_len = sl->l;
                        kh_value(hash, k).max_len = sl->l;
                        kh_value(hash, k).s0 = 1.0;
                        kh_value(hash, k).s1 = sl->s;
                        kh_value(hash, k).s2 = sl->s * sl->s;
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
        */
        /*
        double a,b;
        double av,bv;
        double z_score;
        int outlier = 1;
        while(outlier){
                outlier = 0;
                a = 0.0;
                b = 0.0;
                av = 0.0;
                bv = 0.0;

                RUN(linear_regression(sls->sl, 0 , sls->n, &a, &b));

                RUN(linear_regression_var(sls->sl,0 , sls->n,  a,b, &av,&bv));


                //f_ptr = fopen("GAGA.txt","w");
                for(i = 0; i < sls->n;i++){
                        if(sls->sl[i]->l > 0.0){
                        z_score =   (sls->sl[i]->s - (a + b * sls->sl[i]->l)) / (av + bv * sls->sl[i]->l);
                        if(z_score >= threshold){
                                sls->sl[i]->l = sls->sl[i]->l * -1.0;
                                outlier++;
                        }
                        }
                        //      fprintf(f_ptr,"%f,%f,%f,%f\n",sls->sl[i]->l,sls->sl[i]->s, sls->sl[i]->s - (a + b * sls->sl[i]->l),z_score);
                }
                LOG_MSG("%d %f %f %f %f",outlier,a,b,av,bv);


                //fclose(f_ptr);

        }
        FILE* f_ptr = NULL;
        f_ptr = fopen("GAGA.txt","w");
        for(i = 0; i < sls->n;i++){
                if(sls->sl[i]->l < 0.0){
                        sls->sl[i]->l *= -1.0;
                }
                z_score =   (sls->sl[i]->s - (a + b * sls->sl[i]->l)) / (av + bv * sls->sl[i]->l);
                fprintf(f_ptr,"%f,%f,%f,%f\n",sls->sl[i]->l,sls->sl[i]->s, sls->sl[i]->s - (a + b * sls->sl[i]->l),z_score);
        }
        fclose(f_ptr);

        exit(0);*/
/* merge small bins */
        //RUN(merge_bins(ws_arr, num_hash_items,5000));

        /* average scores  */
        //RUN(average_scores_in_bins(ws_arr, num_hash_items));

        RUN(parcel_out_linear_regression(&p->fit, ws_arr, num_items, sls, threshold));

        RUN(get_dim1(p->fit, &len));


        p->fit_index = NULL;

        p->max_observed_len = (int)( p->fit[len-1][PST_FIT_MAX]);
        RUN(galloc(&p->fit_index, p->max_observed_len+1));


        last = 0;
        for(i = 0; i < len; i++){
                for(j = last;j <= p->fit[i][PST_FIT_MAX];j++){
                        p->fit_index[j] = i;
                }
                last = p->fit[i][PST_FIT_MAX]+1;
        }

        FILE* f_ptr = NULL;
        f_ptr = fopen("GAGA.txt","w");
        for(i = 0; i < sls->n;i++){
                if(sls->sl[i]->l < 0.0){
                        sls->sl[i]->l *= -1.0;
                }
                int len = sls->sl[i]->l;
                int index = MACRO_MIN(p->max_observed_len, len);
                double a,b,v,z_score;
                index = p->fit_index[index];
                a = p->fit[index][PST_FIT_A];
                b = p->fit[index][PST_FIT_B];
                v = p->fit[index][PST_FIT_V];
                z_score = (sls->sl[i]->s - (a + b * (double)len )) / v;

                //z_score =   (sls->sl[i]->s - (a + b * sls->sl[i]->l)) / (av + bv * sls->sl[i]->l);
                fprintf(f_ptr,"%f,%f,%f,%f\n",sls->sl[i]->l,sls->sl[i]->s, sls->sl[i]->s - (a + b * sls->sl[i]->l),z_score);
        }
        fclose(f_ptr);


        for(i = 0; i < num_items;i++){
                MFREE(ws_arr[i]);
        }
        MFREE(ws_arr);
        free_sl(sls);


        kh_destroy(whash,hash);

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

        //float P_M, P_R;

        int chunk,i;//,len;

        RUN(alloc_sl(&sl_store));


        //LOG_MSG("%s",infile);
        if(!my_file_exists(filename)){
                ERROR_MSG("File %s not found");
        }

        RUNP(rng = init_rng(0));

        RUN(open_fasta_fastq_file(&f, filename, TLSEQIO_READ));

        sl_store->offset = 0;
        chunk =1;
        while(1){
                RUN(read_fasta_fastq_file(f, &sb, 1000000));
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
                //LOG_MSG("ALPA: %d %d", sb->L, p->L);
                while(sl_store->n_alloc < sl_store->offset + sb->num_seq){
                        RUN(resize_sl(sl_store));
                        //LOG_MSG("New: %d", sl_store->n_alloc);
                }
                //LOG_MSG("Working on chunk: %d",chunk);
#ifdef HAVE_OPENMP
                omp_set_num_threads(8);
#pragma omp parallel shared(sb,a) private(i)
                {
#pragma omp for schedule(dynamic) nowait
#endif

                        for(i = 0; i < sb->num_seq;i++){
                                int len = sb->sequences[i]->len;
                                float score;
                                convert_to_internal(a, (uint8_t*)sb->sequences[i]->seq,len);
                                score_pst(p, sb->sequences[i]->seq, len, &score);
                                //LOG_MSG("score: %f", P_M-P_R);
                                //P_M = (P_M - P_R)  / 0.69314718055994529;
                                //P_M = MACRO_MAX(P_M, -50.0f);
                                /*if(P_M < -500){
                                        score_pst(p, sb->sequences[i]->seq, len, &P_M,&P_R);
                                        LOG_MSG("weird: %d %f %f", len, P_M,P_R);
                                        }*/
                                sl_store->sl[i + sl_store->offset]->l = (double) len;
                                sl_store->sl[i + sl_store->offset]->s = score;
                                //sl_store->n++;
                                //if(sl_store->n == sl_store->n_alloc){
                                //RUN(resize_sl(sl_store));
                                //}

                        }
#ifdef HAVE_OPENMP
                }
#endif

                sl_store->offset += sb->num_seq;
                sl_store->n = sl_store->offset;

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

int make_intervals(struct  scorelen** arr, int l, double min_ivlen, double min_items,struct window_score*** ws_ret, int* ws_num)
{
        struct window_score** ws_arr = NULL;
//        khash_t(whash) *hash = kh_init(whash);
        //      khiter_t k;
        struct scorelen* sl = NULL;

        double  min_len;
        double max_len;
        double avg_len;
        double old_len;
        double s0;
        double m;
        double av,bv;
        int i,c,ret,num;
        int num_alloc;


        num_alloc = 1000;
        num = 0;

        MMALLOC(ws_arr, sizeof(struct window_score*) * num_alloc);

        s0 = 0.0;
        min_len = 100000000.0;
        max_len = 0.0;
        avg_len = 0.0;

        old_len = 0.0;

        for(i = 0; i < l;i++){
                sl = arr[i];
                if(sl->l > 0.0){

                        if(sl->l !=old_len){
                                if(s0 > min_items && fabs(max_len - min_len) >= min_ivlen){
                                        ws_arr[num] = NULL;
                                        MMALLOC(ws_arr[num], sizeof(struct window_score));


                                        ws_arr[num]->len = avg_len / s0;
                                        ws_arr[num]->min_len = min_len;
                                        ws_arr[num]->max_len = max_len;
                                        ws_arr[num]->s0 = s0;
                                        num++;

                                        if(num == num_alloc){
                                                num_alloc = num_alloc + num_alloc / 2;
                                                MREALLOC(ws_arr, sizeof(struct window_score*) * num_alloc);
                                                for(c = num ; c < num_alloc;c++){
                                                        ws_arr[c] = NULL;
                                                        MMALLOC(ws_arr[c], sizeof(struct window_score));
                                                }
                                        }

                                        /*LOG_MSG("Min: %f max= %f  num:%f",min_len,max_len, s0);
                                        k = kh_put(whash, hash, min_len, &ret);
                                        if (!ret){
                                                LOG_MSG("Weird - already exists... ");
                                        }else{
                                                kh_value(hash, k).len = avg_len /s0;
                                                kh_value(hash, k).min_len = min_len;
                                                kh_value(hash, k).max_len = max_len;
                                                kh_value(hash, k).s0 = s0;
                                                }*/
                                        s0 = 0.0;
                                        min_len = max_len+1;
                                        max_len = 0.0;
                                        avg_len = 0.0;
                                }
                                old_len = sl->l;
                        }
                        min_len = MACRO_MIN(min_len,sl->l);
                        max_len = MACRO_MAX(max_len,sl->l);
                        avg_len += sl->l;
                        s0++;

                }
        }
        ws_arr[num] = NULL;
        MMALLOC(ws_arr[num], sizeof(struct window_score));

        ws_arr[num]->len = avg_len / s0;
        ws_arr[num]->min_len = min_len;
        ws_arr[num]->max_len = max_len;
        ws_arr[num]->s0 = s0;
        num++;

        qsort(ws_arr , num,  sizeof( wscore*),sort_wscore_by_len);

         /* clean up last interval  */
        if(num > 1){

                ws_arr[num-2]->min_len = MACRO_MIN(ws_arr[num-1]->min_len,ws_arr[num-2]->min_len);
                ws_arr[num-2]->max_len = MACRO_MAX(ws_arr[num-1]->max_len,ws_arr[num-2]->max_len);
                ws_arr[num-2]->s0 = ws_arr[num-1]->s0+ws_arr[num-2]->s0;
        }
        MFREE(ws_arr[num-1]);
        num--;

        //LOG_MSG("%f - %f: %f %f",min_len,max_len,avg_len,m);
        //kh_destroy(whash,hash);
        *ws_ret = ws_arr;
        *ws_num = num;
        return OK;
ERROR:
        return FAIL;
}

int merge_bins(struct window_score** arr, int num,int bin_size)
{
        int i,j,c;
        double cur_count;
        double cur_score;
        double cur_score2;
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
        //bin_size = MACRO_MAX(bin_size, cur_max);
        LOG_MSG("%d %d", cur_max, bin_size);
        for(i = 0; i < num;i++){
                LOG_MSG("%d %f",i, arr[i]->len);
        }
        //exit(0);
        for(i = 0; i < num;i++){
                cur_len = arr[i]->len;
                cur_min = arr[i]->min_len;
                cur_max = arr[i]->max_len;
                cur_count = arr[i]->s0;
                cur_score = arr[i]->s1;
                cur_score2 = arr[i]->s2;
                n = 1;
                j = i+1;
                while(cur_count < bin_size  && j != num){
                        cur_len += arr[j]->len;
                        cur_min = MACRO_MIN(cur_min, arr[j]->min_len);
                        cur_max = MACRO_MAX(cur_max, arr[j]->max_len);
                        cur_count += arr[j]->s0;
                        cur_score += arr[j]->s1;
                        cur_score2 += arr[j]->s2;
                        arr[j]->len = 0.0;
                        arr[j]->s0 = 0.0;
                        arr[j]->s1 = 0.0;
                        arr[j]->s2 = 0.0;
                        n = n + 1.0;
                        j++;
                }
                arr[i]->len = cur_len / n;
                arr[i]->min_len = cur_min;
                arr[i]->max_len = cur_max;
                arr[i]->s0 = cur_count;
                arr[i]->s1 = cur_score;
                arr[i]->s2 = cur_score2;
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
                cur_score2  =arr[c]->s2;
                n = 1;
                j = c-1;
                while(cur_count < bin_size && j != 0){
                        if(arr[j]->len){

                                cur_len += arr[j]->len;
                                cur_min = MACRO_MIN(cur_min, arr[j]->min_len);
                                cur_max = MACRO_MAX(cur_max, arr[j]->max_len);
                                cur_count += arr[j]->s0;
                                cur_score += arr[j]->s1;
                                cur_score2 += arr[j]->s2;
                                arr[j]->len = 0.0;
                                arr[j]->s0 = 0.0;
                                arr[j]->s1 = 0.0;
                                arr[j]->s2 = 0.0;
                                n = n + 1.0;
                        }
                        j--;

                }
                arr[c]->len = cur_len / n;
                arr[c]->min_len = cur_min;
                arr[c]->max_len = cur_max;
                arr[c]->s0 = cur_count;
                arr[c]->s1 = cur_score;
                arr[c]->s2 = cur_score2;
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

int parcel_out_linear_regression(double*** fit, struct window_score** arr, int num, struct slen_store* sls,double threshold)
{

        double** fit_param = NULL;
        double s0 = 0.0;
        double s1 = 0.0;
        double s2 = 0.0;
        //double threshold = 10.0;
        double mean;
        double v;
        double z_score;
        double a,b;
        int start_sl = 0;
        int stop_sl = 0;
        int outliers = 1;
        int i,j,c,index;
        int sanity_score_count;
        /* forward across empty buckets
           these come from the last qsort in the merge small bins function
        */

        //LOG_MSG("C: %d - %d",c,num);
        RUN(galloc(&fit_param, num, 5));
        index= 0;
        //f_ptr = fopen("piece.csv", "w");
        for(i = 0;i < num;i++){
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
                //LOG_MSG("LEN: %f %f", arr[i]->min_len, arr[i]->max_len);
                stop_sl = j;
                outliers = 1;

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
                fit_param[index][PST_FIT_MIN] = arr[i]->min_len;
                fit_param[index][PST_FIT_MAX] = arr[i]->max_len;
                fit_param[index][PST_FIT_A] = a;
                fit_param[index][PST_FIT_B] = b;
                fit_param[index][PST_FIT_V] = v;
                index++;

                start_sl = stop_sl;


        }

        *fit = fit_param;
        return OK;
ERROR:
        return FAIL;
}

int average_scores_in_bins(struct window_score** arr, int num)
{
        int i;
        double s0,s1,s2;
        double v,m;
        ASSERT(arr != NULL,"No arr");
        ASSERT(num > 1, "No items");

        for(i = 0; i < num;i++){
                s0 = arr[i]->s0;
                s1 = arr[i]->s1;
                s2 = arr[i]->s2;
                v =  sqrt ( (s0 * s2 -  pow(s1, 2.0)) /  (s0 * ( s0 - 1.0)));
                m = s1 / s0;
                fprintf(stdout,"%d,%f,%f,%f,%f\n",i,s0,arr[i]->len, m,v);
                arr[i]->s1 = arr[i]->s1 / arr[i]->s0;
        }

        linear_regression_bins(arr, num, &s0, &s1);
        for(i = 0; i < num;i++){

                fprintf(stdout,"%d,%f,%f,%f,%f\n",i,arr[i]->len, arr[i]->s1, arr[i]->s1 - (s0 + s1 * arr[i]->len), v);
        }
        LOG_MSG("%f %f",s0,s1);
        exit(0);

        return OK;
ERROR:
        return FAIL;
}

int linear_regression_bins(struct window_score** arr, int num, double* a_var,double* b_var)
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
        n = 0.0;
        for(i = 0; i < num;i++){
                if(arr[i]->len> 0.0){
                        sumX = sumX + arr[i]->len;
                        sumX2 = sumX2 + arr[i]->len * arr[i]->len;
                        sumY = sumY +  arr[i]->s1;
                        sumXY = sumXY + arr[i]->len * arr[i]->s1;
                        n = n + 1.0;
                }
        }

        //n = (double) (num);
        b = (n*sumXY-sumX*sumY)/(n*sumX2-sumX*sumX);
        a = (sumY - b*sumX)/n;
        *a_var = a;
        *b_var = b;
        //LOG_MSG("my_fit <- function(x) { y <- %f + %f * x;return (y)}",a,b);
        return OK;
}

int linear_regression_var(struct scorelen** arr, int s, int e,double a, double b, double* a_var,double* b_var)
{
        struct scorelen* sl = NULL;


        double sumX;
        double sumX2;
        double sumY;
        double sumXY;

        double min_len;
        double max_len;
        double avg_len;
        double s0,s1,s2;
        double m,n;
        double av,bv;
        int i;

        sumX = 0.0;
        sumX2 = 0.0;
        sumY = 0.0;
        sumXY = 0.0;
        n = 0.0;

        s0 = 0.0;
        s1 = 0.0;
        s2 = 0.0;
        min_len = 100000000.0;
        max_len = 0.0;
        avg_len = 0.0;

        for(i = s; i < e;i++){
                sl = arr[i];
                if(sl->l > 0.0){
                        min_len = MACRO_MIN(min_len,sl->l);
                        max_len = MACRO_MAX(max_len,sl->l);
                        avg_len += sl->l;
                        s0++;
                        m = (a +b * sl->l);
                        m = fabs( sl->s - m);
                        s1 += m;
                        s2 += m * m;
                        if(s0 > 5000 && min_len !=  max_len){
                                m =  sqrt ( (s0 * s2 -  pow(s1, 2.0)) /  (s0 * ( s0 - 1.0)));
                                avg_len = avg_len / s0;
                                fprintf(stdout,"%f,%f\n",avg_len,m);
                                //LOG_MSG("%f - %f:%f %f %f",min_len,max_len,s0,avg_len,m);
                                sumX = sumX + avg_len;
                                sumX2 = sumX2 + avg_len * avg_len;
                                sumY = sumY +  m;
                                sumXY = sumXY + avg_len * m;
                                n = n + 1.0;
                                s0 = 0.0;
                                s1 = 0.0;
                                s2 = 0.0;
                                min_len = 100000000.0;
                                max_len = 0.0;
                                avg_len = 0.0;
                                i -= 2500;

                        }
                }
        }
        m =  sqrt ( (s0 * s2 -  pow(s1, 2.0)) /  (s0 * ( s0 - 1.0)));
        avg_len = avg_len / s0;
        fprintf(stdout,"%f,%f\n",avg_len,m);
        //LOG_MSG("%f - %f: %f %f",min_len,max_len,avg_len,m);
        sumX = sumX + avg_len;
        sumX2 = sumX2 + avg_len * avg_len;
        sumY = sumY +  m;
        sumXY = sumXY + avg_len * m;
        n = n + 1.0;
        bv = (n*sumXY-sumX*sumY)/(n*sumX2-sumX*sumX);
        av = (sumY - bv*sumX)/n;
        *a_var = av;
        *b_var = bv;
        LOG_MSG("my_fit <- function(x) { y <- %f + %f * x;return (y)}",av,bv);
        exit(0);
        return OK;

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
        n = 0;
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
                        n = n + 1.0;
                }
        }
        //n = (double) (e -s);
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
        ASSERT(score_arr != NULL, "No scores");
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
        sls->offset = 0;
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
