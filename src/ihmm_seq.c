#include "ihmm_seq.h"


static struct ihmm_sequence* alloc_ihmm_seq(void);
static int realloc_ihmm_seq(struct ihmm_sequence* sequence);
static void free_ihmm_sequence(struct ihmm_sequence* sequence);

int random_label_ihmm_sequences(struct seq_buffer* sb, int k)
{
        int* label;
        int i,j;
        int len;
        //RUN(set_random_seed(void));

        ASSERT(sb != NULL,"No sequences");

        for(i = 0;i< sb->num_seq;i++){

                label = sb->sequences[i]->label;
                len = sb->sequences[i]->seq_len;
                DPRINTF3("Seq%d len %d ",i ,len);
                for(j = 0;j < len;j++){
                        label[j] =  random_int_zero_to_x(k-1) + 2;
                        DPRINTF3("%d",label[j]);
                }
        }
        return OK;
ERROR:
        return FAIL;
}

struct seq_buffer* create_ihmm_sequences_mem(char** seq, int numseq)
{
        struct seq_buffer* sb  = NULL;
        struct ihmm_sequence* sequence = NULL;
        
        int i,j,c,len;
        ASSERT(seq != NULL, "No sequences");
        ASSERT(numseq != 0, "No sequences");
        
        MMALLOC(sb,sizeof(struct seq_buffer));
        sb->malloc_num = 1024;
        sb->num_seq = -1; 
        sb->sequences = NULL;
        sb->max_len = 0;
        
        while(sb->malloc_num <= numseq){
                sb->malloc_num = sb->malloc_num << 1; 
        }
        
        MMALLOC(sb->sequences, sizeof(struct ihmm_sequence*) * sb->malloc_num );
        for(i = 0; i < sb->malloc_num;i++){
                sb->sequences[i] = NULL;
                RUNP(sb->sequences[i] = alloc_ihmm_seq());
        }
        
        /* fill seq */
        for(i = 0; i < numseq;i++){
                
                sequence = sb->sequences[i];
                len = (int) strlen(seq[i]);
                /* Fill name for fun  */
                snprintf(sequence->name, 256, "Sequence_%d", i+1);
                
                sequence->seq_len = 0;
                c = sequence->seq_len;
                for(j = 0; j <  len;j++){
                        switch(seq[i][j]){
                        case 'A':
                        case 'a':
                                sequence->seq[c] = 0;
                                break;
                        case 'C':
                        case 'c':
                                sequence->seq[c] = 1;
                                break;
                        case 'G':
                        case 'g':
                                sequence->seq[c] = 2;
                                break;
                        case 'T':
                        case 't':
                                sequence->seq[c] = 3;
                                break;
                        default:
                                ERROR_MSG("Non ACGT letter in sequence:%d %s.",i,seq[i]);
                                break;
                        }
                        sequence->u[c] = 1.0f;
                        sequence->label[c] = 2;
                        c++;
                        if(c == sequence->malloc_len){
                                RUN(realloc_ihmm_seq(sequence));
                        }
                        
                }
                sequence->seq_len = c;
                if(c > sb->max_len){
                        sb->max_len = c;
                }
        }
        sb->num_seq = numseq;
        return sb;
ERROR:
        free_ihmm_sequences(sb);
        return NULL;
}

struct seq_buffer* load_sequences(char* in_filename)
{
        struct seq_buffer* sb  = NULL;
        struct ihmm_sequence* sequence = NULL;
        FILE* f_ptr = NULL;
        char line[LINE_LEN];
        int i, seq_p;
        

        ASSERT(in_filename != NULL,"No input file specified - this should have been caught before!");


        seq_p = 0;
        MMALLOC(sb,sizeof(struct seq_buffer));
        sb->malloc_num = 1024;
        sb->num_seq = -1; 
        sb->sequences = NULL;
        sb->max_len = 0;
        
               
        MMALLOC(sb->sequences, sizeof(struct chromosome*) *sb->malloc_num );
        for(i = 0; i < sb->malloc_num;i++){
                sb->sequences[i] = NULL;
                RUNP(sb->sequences[i] = alloc_ihmm_seq());                
        }
        
        RUNP(f_ptr = fopen(in_filename, "r" ));
        while(fgets(line, LINE_LEN, f_ptr)){
                if(line[0] == '@' || line[0] == '>'){
                        line[strlen(line)-1] = 0;                       
                        for(i =0 ; i<strlen(line);i++){
                                if(isspace(line[i])){
                                        line[i] = 0;
                                                
                                }
                                        
                        }
                        sb->num_seq++; 

                        if(sb->num_seq == sb->malloc_num){
                                sb->malloc_num = sb->malloc_num << 1;
                                MREALLOC(sb->sequences,sizeof(struct ihmm_sequence*) * sb->malloc_num);
                                for(i = sb->num_seq; i < sb->malloc_num;i++){
                                        sb->sequences[i] = NULL;
                                        RUNP(sb->sequences[i] = alloc_ihmm_seq());     
                                }
                        }
                        sequence = sb->sequences[sb->num_seq];
                        snprintf(sequence->name,256,"%s",line+1);
                        seq_p = 1;
                }else if(line[0] == '+'){
                        seq_p = 0;
                }else{	
                        if(seq_p){
                                for(i = 0;i < LINE_LEN;i++){
                                        if(isalpha((int)line[i])){
                                                switch(line[i]){
                                                case 'A':
                                                case 'a':
                                                        sequence->seq[sequence->seq_len] = 0;
                                                        break;
                                                case 'C':
                                                case 'c':
                                                        sequence->seq[sequence->seq_len] = 1;
                                                        break;
                                                case 'G':
                                                case 'g':
                                                        sequence->seq[sequence->seq_len] = 2;
                                                        break;
                                                case 'T':
                                                case 't':
                                                        sequence->seq[sequence->seq_len] = 3;
                                                        break;
                                                default:
                                                        ERROR_MSG("Non ACGT letter in sequence:%d %s.",i,line[i]);
                                                        break;
                                                }
                                                sequence->u[sequence->seq_len] = 1.0f;
                                                sequence->label[sequence->seq_len] = 2;
                                                sequence->seq_len++;                                              
                                                if(sequence->seq_len == sequence->malloc_len){
                                                        RUN(realloc_ihmm_seq(sequence));
                                                }
                                        }
                                        if(iscntrl((int)line[i])){
                                                if(sequence->seq_len > sb->max_len ){
                                                        sb->max_len = sequence->seq_len; 
                                                }
                                                sequence->seq[sequence->seq_len] = 0;
                                                break;

                                        }
                                }
                        } /* here I would look for quality values.... */
                      
                }
        }
        fclose(f_ptr);
        sb->num_seq++; 
        return sb;
ERROR:
        free_ihmm_sequences(sb);
        if(f_ptr){
                fclose(f_ptr);
        }
        return NULL;
}


struct ihmm_sequence* alloc_ihmm_seq(void)
{
        struct ihmm_sequence* sequence = NULL;
        MMALLOC(sequence,sizeof(struct ihmm_sequence));
        sequence->seq = NULL;
        sequence->u = NULL;
        sequence->label = NULL;
        sequence->name = NULL;
        sequence->malloc_len = 128;
        sequence->seq_len = 0;
        MMALLOC(sequence->seq, sizeof(uint8_t) * sequence->malloc_len);
        MMALLOC(sequence->u, sizeof(float) * (sequence->malloc_len+1));
        MMALLOC(sequence->label , sizeof(int) * sequence->malloc_len);
        MMALLOC(sequence->name, sizeof(char) * 256);
        return sequence;
ERROR:
        free_ihmm_sequence(sequence);
        return NULL;              
}

int realloc_ihmm_seq(struct ihmm_sequence* sequence)
{
        ASSERT(sequence != NULL, "No Sequence.");
        
        sequence->malloc_len = sequence->malloc_len << 1;
        MREALLOC(sequence->seq, sizeof(uint8_t) *sequence->malloc_len);
        MREALLOC(sequence->u, sizeof(float) * (sequence->malloc_len+1));
        MREALLOC(sequence->label , sizeof(int) * sequence->malloc_len);
        
        return OK;
ERROR:
        return FAIL;
}

int print_labelled_ihmm_seq(struct ihmm_sequence* sequence)
{
        int i;
        
        ASSERT(sequence != NULL, "No Sequence.");
        for(i = 0; i < sequence->seq_len;i++){
                fprintf(stdout,"%3c","ACGT"[sequence->seq[i]]);
        }
        fprintf(stdout,"\n");
        for(i = 0; i < sequence->seq_len;i++){
                fprintf(stdout,"%3d",sequence->label[i]);
        }
        fprintf(stdout,"\n");

        return OK;
ERROR:
        return FAIL;
}

void free_ihmm_sequence(struct ihmm_sequence* sequence)
{
        if(sequence){
                if(sequence->seq){
                        MFREE(sequence->seq);
                }
                if(sequence->name){
                        MFREE(sequence->name);
                }
                if(sequence->u){
                        MFREE(sequence->u);
                }
                if(sequence->label){
                        MFREE(sequence->label);
                }
                MFREE(sequence);
        }
}


void free_ihmm_sequences(struct seq_buffer* sb)
{
        int i;
        struct ihmm_sequence* seq = NULL; 
        if(sb){
                for(i =0; i < sb->malloc_num;i++){
                        free_ihmm_sequence(sb->sequences[i]);
                }
                MFREE(sb->sequences);
                MFREE(sb);
        }
}





#ifdef ITESTSEQ

int main(const int argc,const char * argv[])
{
        struct ihmm_sequences* iseq = NULL;
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
        RUNP(iseq = create_ihmm_sequences_mem(tmp_seq ,119));
        free_ihmm_sequences(iseq);
        FILE* f_ptr = NULL;
        RUNP( f_ptr = fopen("ihmm_seq_itest_read_test.fa", "w"));
        for(i = 0; i< 119;i++){
                fprintf(f_ptr,">SEQ_%d\n%s\n",i, tmp_seq[i]);
        }
        fclose(f_ptr);

        iseq = NULL;

        RUNP(iseq = load_sequences("ihmm_seq_itest_read_test.fa"));
        
        free_ihmm_sequences(iseq);
        return EXIT_SUCCESS;
ERROR:
        free_ihmm_sequences(iseq);
        return EXIT_FAILURE;
        
}
#endif
