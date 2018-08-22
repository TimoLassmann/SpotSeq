
#include "pst.h"
#include "ihmm_seq.h"


static int analyze_label_sequences_with_pst(struct seq_buffer* sb);
static int qsort_ints_cmp(const void *a, const void *b);

int count_int_string(int*p,int** suffix,int h,int len);
int binsearch_down(int*p,int** suffix,int h,int len);
int binsearch_up(int*p,int** suffix,int h,int len);
int cmp_int_array(int* a, int*b, int len);

int analyze_label_sequences_with_pst(struct seq_buffer* sb)
{
        struct pst* pst = NULL;
        int i,j,c = 0;
        int** suffix_array = NULL;
        int* alphabet = NULL;
        int* tmp;
        int L;
        int min_L;
        int suffix_len;
        int total_sequence_len;
        /* first challenge: build suffix array from labels  */

        ASSERT(sb != NULL, "No sequence Buffer");

        total_sequence_len = 0;
        L = 0;
        min_L = 1000000;
        for(i = 0 ;i < sb->num_seq;i++){
                total_sequence_len += sb->sequences[i]->seq_len;
                sb->sequences[i]->label[sb->sequences[i]->seq_len] = -1;
                for(j = 0; j < sb->sequences[i]->seq_len;j++){
                        if(sb->sequences[i]->label[j] > L){
                                L = sb->sequences[i]->label[j];
                        }
                        if(sb->sequences[i]->label[j] < min_L){
                                min_L = sb->sequences[i]->label[j];
                        }
                }

        }

        total_sequence_len += sb->num_seq;
        DPRINTF1("Total length is %d", total_sequence_len);
        L = L - min_L;
        MMALLOC(alphabet, sizeof(int) * L);

        for(i = 0; i < L;i++){
                alphabet[i] = i+min_L;
        }
        MMALLOC(suffix_array, sizeof(int*) * total_sequence_len);

        /* assign pointers  */
        suffix_len = 0;
        for(i = 0 ;i < sb->num_seq;i++){
                for(j = 0; j < sb->sequences[i]->seq_len;j++){
                        suffix_array[suffix_len] = &sb->sequences[i]->label[j];
                        suffix_len++;
                }
        }



        qsort(suffix_array, suffix_len, sizeof(int *), qsort_ints_cmp);


        for(i = 0; i < suffix_len;i++){
                tmp= suffix_array[i];
                fprintf(stdout,"i:%d\t",i);
                while(*tmp !=  -1){
                        fprintf(stdout," %2d", *tmp);
                        tmp++;
                }
                fprintf(stdout,"\n");

        }


        pst = malloc(sizeof(struct pst));
        //pst->current_suffix_size = param->num_query* 64;
        pst->suffix_array = suffix_array;
        pst->suffix_len = suffix_len;

        //pst->L = MAX_PST_LEN;
        pst->alpha = 0.0f;
        pst->p_min = 0.0001f;
        pst->lamba = 0.001f;
        pst->r = 1.05f;
        pst->total_len = 0;
        //pst->pst_root = alloc_node(pst->pst_root,"",0);
        //pst->ppt_root = alloc_node(pst->ppt_root,"",0);
        pst->rank_array = 0;
        int* p = NULL;
        MMALLOC(p, sizeof(int)* 5);
        for(i = 0; i < 5;i++){
                p[i] = 11;
        }
        //p[4] = 10;
        fprintf(stdout,"11 5 times occurs %d times.\n",count_int_string(p, pst->suffix_array, pst->suffix_len-1, 5));

        /*

        cStartClock = clock();

        //init root - removes if statement in recursion...
        sum = 0.0;
        for(i = 0;i < 5;i++){
                tmp[0] = alphabet[i];
                tmp[1] = 0;//alphabet[i];
                c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,1);
                pst->pst_root->nuc_probability[i] = c;
                pst->ppt_root->nuc_probability[i] = c;
                sum+= c;
        }
        for(i = 0;i < 5;i++){

                pst->pst_root->nuc_probability[i] =  pst->pst_root->nuc_probability[i]/ sum;
                pst->ppt_root->nuc_probability[i] =  pst->ppt_root->nuc_probability[i]/ sum;
                //fprintf(stderr,"%c\t%f\n",alphabet[i], n->nuc_probability[i]);
        }


        pst->pst_root = build_pst(pst,pst->pst_root );
        pst->ppt_root = build_ppt(pst,pst->ppt_root );

        fprintf(stderr,"built PST in \n");


        pst->pst_root  = alloc_bit_occ_pst(pst->pst_root , numseq);
        pst->ppt_root = alloc_bit_occ_pst(pst->ppt_root, numseq);
        ri =  scan_read_with_pst( ri, pst);
        fprintf(stderr,"scanned  in seconds\n");

        print_pst(pst,pst->pst_root,ri);

        */

        MFREE(suffix_array);
        MFREE(p);
        MFREE(alphabet);
        MFREE(pst);

        //free_pst(pst);
        return OK;
ERROR:
        return FAIL;
}

int count_int_string( int*p, int** suffix,int h,int len)
{
        int a,b;
        //for(i = 0; i < 1000000;i++){
        a = binsearch_down(p,suffix,h,len);
        b = binsearch_up(p,suffix,h,len);
        if(b == h && a != h){
                b++;
        }
        fprintf(stdout,"a: %d b:%d\n",a,b);
        return b-a;
}

int binsearch_down( int*p, int** suffix,int h,int len)
{
        int m = 0;
        int l = 0;
        /*if (t_long_strncmp(p,text+suffix[l],len)<= 0){
          l = l;
          }else */
        if(cmp_int_array(p,suffix[h],len) >  0){
                return h;
        }else{
                while(h-l > 1){
                        //m = (l+h)/2;
                        m = (l + h) >> 1;
                        if(cmp_int_array(p,suffix[m],len) <= 0){
                                h = m;
                        }else{
                                l = m;
                        }
                }
        }
        return l+1;
}

int binsearch_up( int*p, int** suffix,int h,int len)
{
        int m = 0;
        int l = 0;
        /*if (t_long_strncmp(p,text+suffix[l],len)<= 0){
          l = l;
          }else*/
        if(cmp_int_array(p,suffix[h],len) >  0){
                return h;
        }else{
                while(h-l > 1){
                        //m = (l+h)/2;
                        m = (l + h) >> 1;
                        if(cmp_int_array(p,suffix[m],len) < 0){
                                h = m;
                        }else{
                                l = m;
                        }
                }
        }
        return l+1;
}


int cmp_int_array(int* a, int*b, int len)
{
        int i;

        for(i =0; i < len;i++){
                if(a[i] != b[i]){
                        return a[i]-b[i];
                }

                if(a[i] == -1){
                        return 0;
                }
                if(b[i] == -1){
                        return 0;
                }
        }
        return 0;
}


int qsort_ints_cmp(const void *a, const void *b)
{
        const int* one = *(const int**) a;
        const int* two = *(const int**) b;
        int i;
        for (i = 0; one[i] == two[i];i++){
                if(one[i] == -1){
                        return 0;
                }
                if(two[i] == -1){
                        return 0;
                }
        }

        return one[i] - two[i];
}


void free_pst(struct pst_node* n)
{
        int i;
        for(i = 0;i < 5;i++){
                if(n->next[i]){
                        free_pst(n->next[i]);
                }
        }
        if(n->bit_occ){
                free(n->bit_occ);
        }
        free(n->label);
        free(n);

}

void print_pst(struct pst* pst,struct pst_node* n )
{
        int i;
        int internal;

        double p ,e;

        internal = 0;
        for(i = 0;i < 5;i++){
                if(n->next[i]){
                        internal++;
                }
        }
        if(!internal){


                p =  (double) n->occ / (double)  pst->numseq;

                e = p * log2(p);

                p = 1.0 -  (double) n->occ / (double)  pst->numseq;
                e+=  p * log2(p);
                e *= -1.0;
                fprintf(stderr,"%s	%d	%d	Entropy:%f	%f\n", n->label,n->in_T,n->occ,  e,pst->numseq);
        }
        for(i = 0;i < 5;i++){
                if(n->next[i]){
                        if(n->next[i]->in_T){
                                //fprintf(stderr,"Going:%d\n",i);
                                print_pst(pst,n->next[i]);
                        }
                }
        }
}



int main(const int argc,const char * argv[])
{
        struct seq_buffer* sb = NULL;


        int i;
        char *tmp_seq[119] = {
                "ACAGGCTAAAGGAGGGGGCAGTCCCCA",
                "AGGCTAAAGGAGGGGGCAGTCCCCACC",
                "AGGCTAAAGGAGGGGGCAGTCCCCACC",
                "AGTCCCCACCATATTTGAGTCTTTCTC",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGT",
                "CTAAAGGAGGGGGCAGTCCCCACCATA",
                "GAGGCTAAAGGAGGGGGCAGTCCCCAT",
                "GAGGGGGCAGTCCCCACCATATTTGAG",
                "GAGGGGGCAGTCCCCACCATATTTGAG",
                "GAGGGGGCAGTCCCCACCATATTTGAG",
                "GAGGGGGCAGTCCCCACCATATTTGAT",
                "GAGGGGGCAGTCCCCACCATATTTGAT",
                "GAGGGGGCAGTCCCCACCATATTTGAT",
                "GAGTCTTTCTCCAAGTTGCGCCGGACA",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTGGATATCACAGGCTAAAGGAGT",
                "GCTAAAGGAGGGGGCAGTCCCCACCAT",
                "GGAGGGGGCAGTCCCCACCATATTTGA",
                "GGAGGGGGCAGTCCCCACCATATTTGA",
                "GGATATCACAGGCTAAAGGAGGGGGCA",
                "GGCAGTCCCCACCATATTTGAGTCTTC",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGGCAGTCACCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGGCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGGCAGTCCCCACCATATTTGAGTCT",
                "GGGGCAGTCCCCACCATATTTGAGTCT",
                "GGGGCAGTCCCCACCATATTTGAGTCT",
                "GGGGGCAGTCCCCACCATATTTGAGTC",
                "GGGGGCAGTCCCCACCATATTTGAGTC",
                "GTCCCCACCATATTTGAGTCTTTCTCT",
                "TCCCCACCATATTTGAGTCTTTCTCCA",
                "TCCCCACCATATTTGAGTCTTTCTCCA",
                "TGAGTCTTTCTCCAAGTTGCGCCGGAT",
                "TGGATATCACAGGCTAAAGGAGGGGGC"};

        //119l
        RUNP(sb = create_ihmm_sequences_mem(tmp_seq ,119));
        RUN(random_label_ihmm_sequences(sb,10, 0.3f));
        RUN(analyze_label_sequences_with_pst(sb));
        free_ihmm_sequences(sb);



        return EXIT_SUCCESS;
ERROR:
        free_ihmm_sequences(sb);
        return EXIT_FAILURE;

}
//#endif




/*
void pst_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
{
        char alphabet[] = "ACGTN";
        char tmp[MAX_PST_LEN+5];
        struct read_info** ri = 0;
        int i,j,c;
        int numseq;
        clock_t cStartClock;
        FILE* file = 0;
        double sum;
        int num_patterns = 0;
        struct pst* pst = 0;

        pst = malloc(sizeof(struct pst));
        pst->current_suffix_size = param->num_query* 64;
        pst->suffix_array = malloc(sizeof(char*)* pst->current_suffix_size);

        pst->L = MAX_PST_LEN;
        pst->alpha = 0.0f;
        pst->p_min = 0.0001f;
        pst->lamba = 0.001f;
        pst->r = 1.05f;
        pst->total_len = 0;
        pst->pst_root = alloc_node(pst->pst_root,"",0);
        pst->ppt_root = alloc_node(pst->ppt_root,"",0);
        pst->rank_array = 0;

        ri = malloc(sizeof(struct read_info*) * param->num_query);

        for(i = 0; i < param->num_query;i++){
                ri[i] = malloc(sizeof(struct read_info));
                ri[i]->seq = 0;
                ri[i]->name = 0;
                ri[i]->qual = 0;
                ri[i]->len = 0;
                ri[i]->md = 0;
                //ri[i]->xp = 0;
                ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
                ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
                //ri[i]->read_start = -1;
                //ri[i]->read_end = -1;
        }
        file =  io_handler(file, file_num,param);

        struct pst_node** all_patterns = 0;


        while ((numseq = fp(ri, param,file)) != 0){
                pst->total_len = 0;


                for(i = 0; i < numseq;i++){
                        for(j = 0; j < ri[i]->len;j++){
                                ri[i]->seq[j] = alphabet[(int)ri[i]->seq[j]];
                        }
                        //ri[i]->seq[10] = 0;
                        //fprintf(stderr,"%d ",ri[i]->len);
                        pst->total_len += ri[i]->len;
                }


                cStartClock = clock();


                if(pst->current_suffix_size < pst->total_len){
                        pst->suffix_array = realloc(pst->suffix_array , sizeof(char*)* (pst->total_len+64));
                        pst->current_suffix_size =  (pst->total_len+64);
                }
                c = 0;
                pst->mean_length = 0.0;
                for(i = 0; i < numseq;i++){
                        for(j = 0; j < ri[i]->len;j++){//ri[i]->len;j++){
                                pst->suffix_array[c] = ri[i]->seq +j;
                                c++;
                        }
                        pst->mean_length +=  ri[i]->len;

                }

                pst->mean_length /= (float)numseq;


                pst->suffix_len = c;
                pst->numseq = numseq;
                pst->rank_array = malloc(sizeof(struct ranks*) * (int)(numseq));
                for(i = 0; i  < numseq;i++){
                        pst->rank_array[i] = malloc(sizeof(struct ranks));
                }


                fprintf(stderr,"%d\t%d\t%d	%f\n",c,pst->total_len, pst->suffix_len*(4 + 4*5),pst->mean_length);
                ///exit(0);
                qsort(pst->suffix_array, pst->suffix_len, sizeof(char *), qsort_string_cmp);
                fprintf(stderr,"built SA in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);

                cStartClock = clock();


                //init root - removes if statement in recursion...
                sum = 0.0;
                for(i = 0;i < 5;i++){
                        tmp[0] = alphabet[i];
                        tmp[1] = 0;//alphabet[i];
                        c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,1);
                        pst->pst_root->nuc_probability[i] = c;
                        pst->ppt_root->nuc_probability[i] = c;
                        sum+= c;
                }
                for(i = 0;i < 5;i++){

                        pst->pst_root->nuc_probability[i] =  pst->pst_root->nuc_probability[i]/ sum;
                        pst->ppt_root->nuc_probability[i] =  pst->ppt_root->nuc_probability[i]/ sum;
                        //fprintf(stderr,"%c\t%f\n",alphabet[i], n->nuc_probability[i]);
                }


                pst->pst_root = build_pst(pst,pst->pst_root );
                pst->ppt_root = build_ppt(pst,pst->ppt_root );

                fprintf(stderr,"built PST in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);


                cStartClock = clock();
                pst->pst_root  = alloc_bit_occ_pst(pst->pst_root , numseq);
                pst->ppt_root = alloc_bit_occ_pst(pst->ppt_root, numseq);
                ri =  scan_read_with_pst( ri, pst);
                fprintf(stderr,"scanned  in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);

                print_pst(pst,pst->pst_root,ri);


                num_patterns = 0;
                num_patterns = count_patterns(pst->pst_root, num_patterns);
                num_patterns = count_patterns(pst->ppt_root,num_patterns);
                //
                //exit(0);
                all_patterns = malloc(sizeof(struct pst_node*)  *num_patterns );

                num_patterns = 0;
                num_patterns = add_patterns(all_patterns,pst->pst_root, num_patterns);
                num_patterns = add_patterns(all_patterns,pst->ppt_root,num_patterns);

                qsort((void *)  all_patterns, num_patterns, sizeof(struct pst_node* ),(compfn) sort_pst_nodel_according_to_label);

                fprintf(stderr,"%d numpatterns\n",num_patterns );
                c = 0;
                for(i = 0; i < num_patterns-1;i++){
                        if(strcmp(all_patterns[i]->label, all_patterns[i+1]->label)){
                                all_patterns[c] = all_patterns[i];
                                c++;
                        }
                        //fprintf(stderr,"%s	%d\n",all_patterns[i]->label,i );
                }

                //for(i = 0; i < c;i++){
                //	fprintf(stderr,"%s	%d\n",all_patterns[i]->label,i );
                //}
                // got all strings combined and ready for dbscan.....

                cStartClock = clock();
                cluster_reads_based_on_pst_patterns(all_patterns,c,numseq,ri);
                fprintf(stderr,"clustered   in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);

                //print_pst(pst,pst->pst_root,ri);

                exit(0);
                run_gmm_on_sequences(ri,numseq);



                print_pst(pst,pst->pst_root,ri);
                fprintf(stderr,"happily got here...\n");
                print_pst(pst,pst->ppt_root,ri);
                exit(0);
                //cStartClock = clock();
                //ri =  scan_read_with_pst( ri, pst);
                //fprintf(stderr,"scanned  in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);

                for(i = 0; i  <pst->numseq;i++){
                        free(pst->rank_array[i]);// = malloc(sizeof(struct ranks));
                }
                free(pst->rank_array);// = malloc(sizeof(struct ranks*) * (int)(numseq));
                exit(0);


                cStartClock = clock();
                //pst->pst_root = count_patterns(ri,pst,pst->pst_root);
                fprintf(stderr,"counted  in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);
                print_pst(pst,pst->pst_root,ri);

                exit(0);



        }


        for(i = 0; i < param->num_query;i++){
                free(ri[i]->strand);
                free(ri[i]->hits);
                if(ri[i]->md){
                        free(ri[i]->md);
                }

                free(ri[i]);
        }
        free(ri);
        //fprintf(stderr,"%p\n",file);
        if(param->sam == 2 || param->sam == 1 || param->gzipped ){
                pclose(file);
        }else{
                //if(file_num != -1){
                fclose(file);
                //}
        }

}


struct pst_node* build_pst(struct pst* pst,struct pst_node* n )
{
        char alphabet[] = "ACGTN";

        char tmp[MAX_PST_LEN+5];
        int i;
        int j;
        int c;
        int add;
        int len = (int) strlen(n->label);
        double sum = 0.0f;



        double tmp_counts_s[5];



        //fprintf(stderr,"NODE: %s\n", n->label);
        //for(i = 0;i < 5;i++){
        //	fprintf(stderr,"%c+%s\t%f\n",alphabet[i], n->label, n->nuc_probability[i]);
        //}

        //step 2 test expansion

        //loop though letters at present node
        if(len + 1 < MAX_PST_LEN ){
                /// search for all strings and record probabilities S+ACGT...
                /// don't search rare strings...
                /// - super mega simple ...


                for(i = 0; i < 5;i++){
                        if(n->nuc_probability[i] >= pst->p_min){ /// string (suffix  + letter is not rare compared to all other letters  ) // e.g. acgC = 0.2 and not acgC = 0.0000001. hence it makes sense to "extend the suffix - test [A,C,G,T] - acgC





                                //init longer suffix
                                tmp[0] = alphabet[i];
                                for(j = 1; j < len+1;j++){
                                        tmp[j] = n->label[j-1];
                                }

                                sum = 0.0;
                                for(j = 0; j < 5;j++){
                                        tmp[len+1]  = alphabet[j];
                                        tmp[len+2] = 0;
                                        c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,len+2);
                                        tmp_counts_s[j] = c;
                                        sum+= c;
                                }

                                add = 0;
                                for(j = 0; j < 5;j++){
                                        if(tmp_counts_s[j] /  (pst->numseq * (pst->mean_length - (len+1))) > 0.001){
                                                add = 1;
                                                break;
                                        }
                                }
                                if(add){
                                        // here I know that probablility of 'X' is non-neglible AND that there existsa string 'X' - [A,C,G,T,N] which is frequent in the data - hence I add...
                                        n->next[i] = alloc_node(n->next[i] ,tmp,len+1);
                                        add = 0;
                                        for(j = 0; j < 5;j++){
                                                if((tmp_counts_s[j]/sum) / n->nuc_probability[j] >= pst->r){
                                                        add++;
                                                }

                                                if((tmp_counts_s[j]/sum) / n->nuc_probability[j] <= 1.0/ pst->r){
                                                        add++;
                                                }

                                                n->next[i]->nuc_probability[j] = tmp_counts_s[j]/sum;

                                        }
                                        if(add){
                                                n->next[i]->in_T = 1;
                                        }
                                        n->next[i] = build_pst(pst,n->next[i]  );

                                }
                        }
                }
        }
        c= 0;
        for(i = 0; i < 5;i++){
                if(n->next[i]){
                        c+= n->next[i]->in_T;
                }
        }
        if(c){
                n->in_T = 1;
        }
        return n;
}



struct pst_node* build_ppt(struct pst* pst,struct pst_node* n )
{
        char alphabet[] = "ACGTN";

        char tmp[MAX_PST_LEN+5];
        int i;
        int j;
        int c;
        int add;
        int len = (int) strlen(n->label);
        double sum = 0.0f;



        double tmp_counts_s[5];



        //fprintf(stderr,"NODE: %s\n", n->label);
        //for(i = 0;i < 5;i++){
        //	fprintf(stderr,"%c+%s\t%f\n",alphabet[i], n->label, n->nuc_probability[i]);
        //}

        //step 2 test expansion

        //loop though letters at present node
        if(len + 1 < MAX_PST_LEN ){
                /// search for all strings and record probabilities S+ACGT...
                /// don't search rare strings...
                /// - super mega simple ...


                for(i = 0; i < 5;i++){
                        if(n->nuc_probability[i] >= pst->p_min){ /// string (suffix  + letter is not rare compared to all other letters  ) // e.g. acgC = 0.2 and not acgC = 0.0000001. hence it makes sense to "extend the suffix - test [A,C,G,T] - acgC





                                //init longer prefix!!!!

                                for(j = 0; j < len;j++){
                                        tmp[j+1] = n->label[j];
                                }
                                tmp[len+1] = alphabet[i];

                                sum = 0.0;
                                for(j = 0; j < 5;j++){
                                        tmp[0]  = alphabet[j];
                                        tmp[len+2] = 0;
                                        c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,len+2);
                                        tmp_counts_s[j] = c;
                                        sum+= c;
                                }

                                add = 0;
                                for(j = 0; j < 5;j++){
                                        if(tmp_counts_s[j] /  (pst->numseq * (pst->mean_length - (len+1))) > 0.001){
                                                add = 1;
                                                break;
                                        }
                                }
                                if(add){
                                        // here I know that probablility of 'X' is non-neglible AND that there existsa string 'X' - [A,C,G,T,N] which is frequent in the data - hence I add...
                                        n->next[i] = alloc_node(n->next[i],tmp+1,len+1);
                                        add = 0;
                                        for(j = 0; j < 5;j++){
                                                if((tmp_counts_s[j]/sum) / n->nuc_probability[j] >= pst->r){
                                                        add++;
                                                }

                                                if((tmp_counts_s[j]/sum) / n->nuc_probability[j] <= 1.0/ pst->r){
                                                        add++;
                                                }

                                                n->next[i]->nuc_probability[j] = tmp_counts_s[j]/sum;

                                        }
                                        if(add){
                                                n->next[i]->in_T = 1;
                                        }
                                        n->next[i] = build_ppt(pst,n->next[i]);

                                }
                        }
                }
        }
        c= 0;
        for(i = 0; i < 5;i++){
                if(n->next[i]){
                        c+= n->next[i]->in_T;
                }
        }
        if(c){
                n->in_T = 1;
        }

        return n;
}


struct pst_node* alloc_node(struct pst_node* n,char* string,int len)
{
        int i;
        n = malloc(sizeof(struct pst_node));

        assert(n!=0);

        n->label = malloc(sizeof(char) *(len+1));
        assert(n->label != 0);

        for(i = 0; i < len;i++){
                n->label[i] =string[i];
        }
        n->label[len] = 0;
        n->in_T = 0;
        n->occ = 0;
        n->last_seen = -1;
        //n->bit_occ = 0;
        for(i =0; i < 5;i++){
                n->next[i] = 0;
                n->nuc_probability[i] = 0.2f;
        }

        return n;
}


struct read_info**  scan_read_with_pst(struct read_info** ri,struct pst* pst)
{
        int i,j;
        float P_T;
        //float P_PT;
        float P_R = 0;

        float* base_p = pst->pst_root->nuc_probability;
        char* seq;
        char* qual;

        float total_T = prob2scaledprob(1.0);
        float total_R = prob2scaledprob(1.0);

        float A,B;

        //scan to cpount occurances...
        for(i = 0; i < pst->numseq;i++){

                seq = ri[i]->seq;
                for(j = 0; j < ri[i]->len; j++ ){
                        pst->pst_root = count_pst_lables(pst->pst_root, seq,  j, i);
                        pst->ppt_root = count_ppt_lables(pst->ppt_root, seq, j , i);

                }
        }

        //fprintf(stderr,"%f\n",pst->numseq);

        for(i = 0; i < pst->numseq;i++){
                if(!ri[i]->qual){
                        ri[i]->qual = malloc(sizeof(char)* (ri[i]->len+1));
                }
                for(j = 0; j < ri[i]->len; j++ ){
                        ri[i]->qual[j] = 48;
                }
                qual = ri[i]->qual;
                seq = ri[i]->seq;

                //if(!i){

                //}

                P_T = prob2scaledprob(1.0);
                P_R = prob2scaledprob(1.0);
                //P_PT =  prob2scaledprob(1.0);
                //fprintf(stdout,"%s\n",seq );

                for(j = 0; j < ri[i]->len; j++ ){
                        P_R = P_R + prob2scaledprob(base_p[nuc_code[(int)seq[j]]]);


                        //get_occ
                        A = get_pst_prob(pst->pst_root, seq,  nuc_code[(int)seq[j]], j, i);
                        B = get_ppt_prob(pst->ppt_root, seq,  nuc_code[(int)seq[j]], j, i);


                        P_T = P_T + prob2scaledprob(max(A,B));

                        //fprintf(stdout,"%d	%f	%f\n", j , A,B);
                        ri[i]->qual[j] = 48 + (int) (10.0* (max(A,B)));
                }

                total_T = total_T + P_T;
                total_R = total_R + P_R;

                ri[i]->mapq = P_T-P_R;
        }
        return ri;
}



int add_patterns(struct pst_node** all_patterns, struct pst_node* n,int num)
{
        int i;

        int internal = 0;
        for(i = 0;i < 5;i++){
                if(n->next[i]){
                        internal++;
                }
        }
        if (!internal){
                if(strlen(n->label) > 2){
                        if(n->in_T){
                                all_patterns[num] = n;
                                num += 1;
                                //fprintf(stderr,"ADDED: %s	%d\n",n->label,num);
                        }
                }

        }

        for(i = 0;i < 5;i++){
                if(n->next[i]){
                        //fprintf(stderr,"Going:%d\n",i);
                        num = add_patterns(all_patterns,n->next[i],num);
                }
        }

        return num;
}

int count_patterns(struct pst_node* n,int num)
{
        int i;
        if(strlen(n->label) > 2){
                if(n->in_T){
                        num += 1;
                        //fprintf(stderr,"ADDED: %s	%d\n",n->label,num);
                }
        }

        for(i = 0;i < 5;i++){
                if(n->next[i]){
                        //fprintf(stderr,"Going:%d\n",i);
                        num = count_patterns(n->next[i],num);
                }
        }

        return num;
}

float get_pst_prob(struct pst_node* n, char* string,int target, int pos,int seq_id)
{
        int c;
        //if(!seq_id){
        //fprintf(stderr,"%d	%s	%f	(%d)	- OCC:%d\n", target,n->label,n->nuc_probability[target],n->in_T, n->occ);
        //}

        if(pos == 0){

                return n->nuc_probability[target];
        }
        pos = pos -1;
        c = nuc_code[(int)string[pos]];
        if(n->next[c]){
                return get_pst_prob(n->next[c], string, target,pos,seq_id);
        }else{

                return n->nuc_probability[target];
        }
}




float get_ppt_prob(struct pst_node* n, char* string,int target, int pos,int seq_id)
{
        int c;
        //if(!seq_id){
        //fprintf(stderr,"%d	%s	%f	(%d)	- OCC:%d\n", target,n->label,n->nuc_probability[target],n->in_T, n->occ);
        //}

        if(string[pos+1] == 0){

                return n->nuc_probability[target];
        }
        pos = pos +1;
        c = nuc_code[(int)string[pos]];
        if(n->next[c]){
                return get_ppt_prob(n->next[c], string, target,pos,seq_id);
        }else{

                return n->nuc_probability[target];
        }
}

int get_occ(struct pst_node* n, char* string,int target, int pos,int seq_id)
{
        //if(!seq_id){
        //fprintf(stderr,"%d	%s	%f	(%d)	- OCC:%d\n", target,n->label,n->nuc_probability[target],n->in_T, n->occ);
        //}

        if(pos == 0){
                return n->occ;
        }
        pos = pos -1;
        int c = nuc_code[(int)string[pos]];
        if(n->next[c]){
                return get_occ(n->next[c], string, target,pos,seq_id);
        }else{
                return n->occ;
        }
}


struct pst_node*  count_pst_lables(struct pst_node* n, char* string, int pos,int seq_id)
{
        //if(!seq_id){
        //fprintf(stderr,"%d	%s	%f	(%d)	- OCC:%d\n", target,n->label,n->nuc_probability[target],n->in_T, n->occ);
        //}

        int c = nuc_code[(int)string[pos]];
        if(n->next[c]){
                if(n->next[c]->last_seen != seq_id){
                        if(!bit_test(n->next[c]->bit_occ , seq_id)){
                                n->next[c]->occ++;
                                n->next[c]->bit_occ = bit_set(n->next[c]->bit_occ , seq_id);
                        }
                        n->next[c]->last_seen = seq_id;
                }
                pos = pos -1;
                if(pos != -1){
                        n->next[c] =  count_pst_lables(n->next[c], string,pos,seq_id);
                }
        }else{
                return n;
        }


        return n;
}


struct pst_node*  count_ppt_lables(struct pst_node* n, char* string, int pos,int seq_id)
{
        int c = nuc_code[(int)string[pos]];
        if(!string[pos]){
                return n;
        }
        if(n->next[c]){
                if(n->next[c]->last_seen != seq_id){
                        if(!bit_test(n->next[c]->bit_occ , seq_id)){
                                n->next[c]->occ++;
                                n->next[c]->bit_occ = bit_set(n->next[c]->bit_occ , seq_id);
                        }
                        n->next[c]->last_seen = seq_id;
                }
                pos = pos +1;
                //n = n->next;
                //if(pos){
                n->next[c] =  count_ppt_lables(n->next[c], string,pos,seq_id);
                //}
        }else{
                return n;
        }


        return n;
}






struct pst_node* alloc_bit_occ_pst(struct pst_node* n, int num)
{
        int i;
        n->bit_occ = _mm_malloc(sizeof(int)* (1+ num / BITSPERWORD) , 16);

        assert(n->bit_occ != 0);

        for(i = 0; i < (1+ num / BITSPERWORD);i++){
                n->bit_occ[i] = 0;
        }

        for(i = 0;i < 5;i++){
                if(n->next[i]){
                        n->next[i] = alloc_bit_occ_pst(n->next[i],num);
                }
        }
        return n;
}


int* bit_set(int*a, int i)
{
        a[i >> SHIFT] |= (1 <<(i & MASK));
        return a;
}


int* bit_clr(int*a, int i)
{
        a[i >> SHIFT] &= ~(1 <<(i & MASK));
        return a;
}

int bit_test(int*a, int i)
{
        return a[i >> SHIFT] & (1 <<(i & MASK));
}


int establish_rank(const void *a, const void *b)
{
        struct ranks* const *one = a;
        struct ranks* const *two = b;

        if((*one)->value >  (*two)->value){
                return -1;
        }else{
                return 1;
        }
}

int sort_pst_nodel_according_to_label(const void *a, const void *b)
{
        struct pst_node * const *one = a;
        struct pst_node* const *two = b;

        return strcmp((*one)->label,(*two)->label) ;
}






*/
