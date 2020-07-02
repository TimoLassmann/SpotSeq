
#include "tldevel.h"
#include "pst.h"

#define PST_HASH_IMPORT
#include "pst_hash.h"

static int sanity_check_hash(struct count_hash* hash);

int fill_exact_hash(struct count_hash** hash, struct tl_seq_buffer* sb)
{

        struct count_hash* h_struct  = NULL;
        khash_t(exact) *h = NULL;
        khiter_t k;
        char* seq;
        int len;
        uint64_t* mask = NULL;
        uint64_t key;
        uint64_t c;
        uint64_t x;
        int i,j;
        int ret;
        uint64_t max_len,l;

        MMALLOC(h_struct,sizeof(struct count_hash));

        h_struct->hash = NULL;
        h_struct->mask = NULL;
        h_struct->len = MAX_PST_LEN;
        h_struct->L = 0;
        h_struct->counts_l = NULL;

        if(sb->L == TL_SEQ_BUFFER_DNA){
                h_struct->L = 4;
        }else if(sb->L == TL_SEQ_BUFFER_PROTEIN){
                h_struct->L = 20;
        }

        MMALLOC(h_struct->mask, sizeof(uint64_t) * MAX_PST_LEN);

        MMALLOC(h_struct->counts_l, sizeof(float) * (MAX_PST_LEN+1));

        for(i = 0;i <= h_struct->len;i++){
                h_struct->counts_l[i] = 0.0f;
        }
        /* count  how many strings of length i there are in the set of sequences */
        for(i = 1;i <= h_struct->len;i++){
                for(j = 0; j < sb->num_seq;j++){
                        if(sb->sequences[j]->len >= i){
                                h_struct->counts_l[i] += (sb->sequences[j]->len  - (i -1));
                        }
                }
        }

        h_struct->min_prop = 0.0f;
        for(i = 1;i <= h_struct->len;i++){
                if(2.0 / h_struct->counts_l[i] > h_struct->min_prop){
                        h_struct->min_prop = 2.0 / h_struct->counts_l[i];
                }
                //LOG_MSG("Min prob: %f %f %f ", 2.0 / h_struct->counts_l[i], h_struct->counts_l[i], h_struct->min_prop);
        }

        h_struct->mask[0] = 0x1FULL;
        for(i = 1;i < h_struct->len;i++){
                h_struct->mask[i] = (h_struct->mask[i-1] << 5ULL) | 0x1FULL;
        }


        h = h_struct->hash;
        if(!h){
                h = kh_init(exact);
        }

        mask = h_struct->mask;

        max_len = h_struct->len;
        //LOG_MSG("HASH L: %d", h_struct->L);
        for(i = 0; i < sb->num_seq;i++){
                x = 0ULL;
                seq = sb->sequences[i]->seq;
                len = sb->sequences[i]->len;
                //fprintf(stdout,"%s\n", seq);
                l = 0;
                for(j = 0; j < len;j++){
                        x = (x << 5ULL) | (uint64_t)seq[j];
                        for(c = 0; c <= l;c++){
                                key = ((c+1) << 60ULL) | (x & mask[c]); /* first 4 bits is length  */
                                //LOG_MSG("%lud", key);
                                k = kh_put(exact, h, key, &ret);
                                if(!ret){
                                        //k = kh_get(32, h, r);
                                        //fprintf(stdout,"Exists - increment\n");
                                        kh_value(h, k)++;
                                }else{
                                        kh_value(h, k) =1;
                                        //fprintf(stdout,"set to 1 \n");
                                }
                        }
                        l++;
                        l = MACRO_MIN(l, max_len-1);
                }
        }
        h_struct->hash = h;

        *hash = h_struct;
        /*i = 0;
        for (k = kh_begin(h); k != kh_end(h); ++k){
                if (kh_exist(h, k)){
                        //dkey = kh_key(h, k);
                        //fprintf(stdout, "%lud %d\n",kh_key(h,k),kh_value(h, k));
                        i += kh_value(h, k);
                        //kh_value(h, k) = 1;
                }
        }
        */

        //fprintf(stdout," hash size: %d, count: %d L:%d\n", kh_size(h),i, h_struct->L);

        //exit(0);


        RUN(sanity_check_hash(h_struct));
        *hash = h_struct;
        return OK;
ERROR:
        return FAIL;
}


int sanity_check_hash(struct count_hash* hash)
{
        khiter_t hi;
        uint64_t key;

        int i,j,c;
        int sum;



        sum = 0;
        for(i = 0; i < hash->L;i++){
                for(j = 0; j < hash->L;j++){
                        key = (2ULL << 60ULL) | ((uint64_t) i << 5ULL) | (uint64_t)j ;
                        c = 0;
                        hi= kh_get(exact, hash->hash, key);
                        if(hi != kh_end(hash->hash)){
                                //fprintf(stdout, "%lud %d\n",kh_key(hash->hash,hi),kh_value(hash->hash, hi));
                                //LOG_MSG("%c: %f", "ACGT"[i],kh_value(h->hash, hi));
                                c  = kh_value(hash->hash, hi);
                        }
                        sum += c;
                        /*if(hash->L == 4){
                                fprintf(stdout,"%c","ACGT"[i]);
                                fprintf(stdout,"%c","ACGT"[j]);
                                fprintf(stdout,"\t%d (sum: %d)\n",c,sum);

                        }else{
                                fprintf(stdout,"%c","ACDEFGHIKLMNPQRSTVWY"[i]);
                                fprintf(stdout,"%c","ACDEFGHIKLMNPQRSTVWY"[j]);
                                fprintf(stdout,"\t%d (sum: %d)\n",c,sum);

                                }*/
                }
        }

        //fprintf(stdout,"counted %d <-> %f",sum, hash->counts_l[2]);
        //ASSERT(fabsf((float)sum - hash->counts_l[2]) < 0.5f,"counts don't add up: %d %f",sum,hash->counts_l[2]);
        return OK;

ERROR:
        return FAIL;
}


void free_exact_hash(struct count_hash* hash)
{
        if(hash){
                kh_destroy(exact, hash->hash);
                MFREE(hash->counts_l);
                MFREE(hash->mask);
                MFREE(hash);
        }
}
