#include "ihmm_seq.h"

int random_label_ihmm_sequences(struct ihmm_sequences* iseq, int k)
{
        int i,j;
        //RUN(set_random_seed(void));

        ASSERT(iseq != NULL,"No sequences");

        for(i = 0;i< iseq->numseq;i++){
                for(j = 0;j < iseq->len[i];j++){
                        iseq->label[i][j] = random_int_zero_to_x(k-1) + 2;
                        DPRINTF1("%d",iseq->label[i][j]);
                }
        }
        return OK;
ERROR:
        return FAIL;
}





struct ihmm_sequences* create_ihmm_sequences(char** seq, int numseq)
{
        struct ihmm_sequences* iseq = NULL;
        int i,j;
        ASSERT(seq != NULL, "No sequences");
        MMALLOC(iseq, sizeof(struct ihmm_sequences));
        iseq->label = NULL;
        iseq->len = NULL;
        iseq->seq = NULL;
        iseq->u = NULL;
        iseq->numseq = numseq;
        iseq->max_len = 0;
        MMALLOC(iseq->len,sizeof(uint32_t) * numseq);

        iseq->max_len = 0;
        for(i = 0;i < numseq;i++){
                iseq->len[i] = (uint32_t) strlen(seq[i]);
                if(iseq->len[i] > iseq->max_len){
                        iseq->max_len = iseq->len[i];
                }
        }

        ASSERT(iseq->max_len != 0, "Weird length of sequences is zero");

        /* alloc arrays */
        RUNP(iseq->u = malloc_2d_float(iseq->u, numseq, iseq->max_len, 0.0f));
        RUNP(iseq->label = malloc_2d_int(iseq->label, numseq, iseq->max_len, 0.0f));
        RUNP(iseq->seq = malloc_2d_char(iseq->seq, numseq, iseq->max_len, 0.0f));


        /* fill seq */
        for(i = 0; i < numseq;i++){
                for(j = 0; j < iseq->len[i];j++){
                        switch(seq[i][j]){
                        case 'A':
                        case 'a':
                                iseq->seq[i][j] = 0;
                                break;
                        case 'C':
                        case 'c':
                                iseq->seq[i][j] = 1;
                                break;
                        case 'G':
                        case 'g':
                                iseq->seq[i][j] = 2;
                                break;
                        case 'T':
                        case 't':
                                iseq->seq[i][j] = 3;
                                break;
                        default:
                                ERROR_MSG("Non ACGT letter in sequence:%d %s.",i,seq[i]);
                                break;
                        }
                }
        }
        return iseq;
ERROR:
        return NULL;
}


void free_ihmm_sequences(struct ihmm_sequences* iseq)
{
        if(iseq){
                if(iseq->label){
                        free_2d((void**)iseq->label);
                }

                if(iseq->u){
                        free_2d((void**)iseq->u);
                }

                if(iseq->seq){
                        free_2d((void**)iseq->seq);
                }

                MFREE(iseq->len);
                MFREE(iseq);
        }
}






#ifdef ITESTSEQ

int main(const int argc,const char * argv[])
{
        struct ihmm_sequences* iseq = NULL;
        
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
        RUNP(iseq = create_ihmm_sequences(tmp_seq ,119));
                
        free_ihmm_sequences(iseq);
        return EXIT_SUCCESS;
ERROR:
        free_ihmm_sequences(iseq);
        return EXIT_FAILURE;
        
}
#endif
