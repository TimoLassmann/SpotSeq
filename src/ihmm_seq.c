#include "ihmm_seq.h"


static struct ihmm_sequence* alloc_ihmm_seq(void);
static int realloc_ihmm_seq(struct ihmm_sequence* sequence);
static void free_ihmm_sequence(struct ihmm_sequence* sequence);


int detect_alphabet(struct seq_buffer* sb);

int translate_DNA_to_internal(struct seq_buffer* sb);
int translate_internal_to_DNA(struct seq_buffer* sb);

int translate_internal_to_PROTEIN(struct seq_buffer* sb);
int translate_PROTEIN_to_internal(struct seq_buffer* sb);

int compare_sequence_buffers(struct seq_buffer* a, struct seq_buffer* b);



/* The idea here is to assume that there are K states with emission
 * probabilities samples from a dirichlet distribution. Labels will be chiosen
 * by sampling from the distributions samples by the dirichlets. */

int dirichlet_emission_label_ihmm_sequences(struct seq_buffer* sb, int k, float alpha)
{
        float** emission = NULL;
        int i,j,c;
        uint8_t* seq;
        int* label;
        int len;
        rk_state rndstate;
        float sum = 0;

        float r;


        

        ASSERT(sb != NULL, "No sequence buffer");

        rk_randomseed(&rndstate);


        //allocfloat** malloc_2d_float(float**m,int newdim1, int newdim2,float fill_value)

        RUNP(emission = malloc_2d_float(emission, k+1,  sb->L , 0.0f));

        for(i = 0; i < k;i++){
                sum = 0.0;
                for(j = 0;j < sb->L;j++){
                        emission[i][j] = rk_gamma(&rndstate,alpha , 1.0);
                        sum += emission[i][j];
                }
                for(j = 0;j < sb->L;j++){
                        emission[i][j] /= sum;
                        emission[k][j] += emission[i][j]; /* Last row has sums of all emissions of Letter 0, 1, 2, .. L  *\/ */
                }
                
        }
       
        for(i = 0;i< sb->num_seq;i++){

                label = sb->sequences[i]->label;
                len = sb->sequences[i]->seq_len;
                seq = sb->sequences[i]->seq;
                
                for(j = 0;j < len;j++){
                        r = random_float_zero_to_x(emission[k][seq[j]]);
                        for(c = 0; c < k;c++){
                                r -= emission[c][seq[j]];
                                if(r <= 0.0f){
                                        label[j] = c+2;
                                        break;
                                }
                                                 
                        }
                }
        }

        free_2d((void**) emission);
        
        return OK;
ERROR:
        if(emission){
                free_2d((void**) emission);
        }
        return FAIL;
}


int label_ihmm_sequences_based_on_guess_hmm(struct seq_buffer* sb, int k, float alpha)
{
        float** transition = NULL;
        float** emission = NULL;
        float* tmp = NULL;
        int i,j,c;
        uint8_t* seq;
        int* label;
        int len;
        rk_state rndstate;
        float sum = 0;
        float sanity; 

        float r;
        int n;
        int cur_state;
        

        ASSERT(sb != NULL, "No sequence buffer");

        n = sb->L;
        if(sb->L == ALPHABET_DNA){
                n = 4;
        }
        rk_randomseed(&rndstate);


        //allocfloat** malloc_2d_float(float**m,int newdim1, int newdim2,float fill_value)

        RUNP(emission = malloc_2d_float(emission, k+1,  sb->L , 0.0f));
        RUNP(transition = malloc_2d_float(transition, k+1,  k , 0.0f));

        MMALLOC(tmp, sizeof(float) * k);
        //fprintf(stdout,"Emission\n");
        for(i = 0; i < k;i++){
                sum = 0.0;
                for(j = 0;j < n;j++){
                        emission[i][j] = rk_gamma(&rndstate,alpha , 1.0);
                        sum += emission[i][j];
                }
                sanity = 0.0f;
                for(j = 0;j < n;j++){
                        emission[i][j] /= sum;
                        emission[k][j] += emission[i][j]; /* Last row has sums of all emissions of Letter 0, 1, 2, .. L  *\/ */
                        sanity += emission[i][j];
                        //fprintf(stdout,"%f ",emission[i][j]);
                }
                //fprintf(stdout,"sum: %f\n",sanity); 
        }
        //fprintf(stdout,"Transition\n");
        for(i = 0; i < k;i++){
                sum = 0.0;
                for(j = 0; j < k;j++){
                       
                        transition[i][j] = rk_gamma(&rndstate,alpha , 1.0);
                        
                        sum += transition[i][j];
                }
                sanity = 0.0f;
                for(j = 0; j < k;j++){
                        transition[i][j] /= sum;
                        transition[k][j] += transition[i][j];
                        sanity += transition[i][j];
                        // fprintf(stdout,"%f ",transition[i][j]);
                }
                // fprintf(stdout,"sum: %f\n", sanity); 
        }
        //exit(0);
        
        for(i = 0;i< sb->num_seq;i++){

                label = sb->sequences[i]->label;
                len = sb->sequences[i]->seq_len;
                seq = sb->sequences[i]->seq;
                r = random_float_zero_to_x(emission[k][seq[0]]);
                for(c = 0; c < k;c++){
                        r -= emission[c][seq[0]];
                        if(r <= 0.0f){
                                label[0] = c+2;
                                cur_state = c;
                                break;
                        }
                                                 
                }
                
                for(j = 1;j < len;j++){
                        sum = 0.0f;
                        for(c = 0; c < k;c++){
                                tmp[c] = transition[cur_state][c] * emission[c][seq[j]];
                                sum += tmp[c];
                        }
                        r = random_float_zero_to_x(sum);
                        for(c = 0; c < k;c++){
                                r -= tmp[c];
                                if(r <= 0.0f){
                                        label[j] = c+2;
                                        cur_state = c;
                                        break;
                                }
                                                 
                        }
                }
        }
        for(i = 0;i< sb->num_seq;i++){
                label = sb->sequences[i]->label;
                len = sb->sequences[i]->seq_len;
                seq = sb->sequences[i]->seq;
                for(j = 0; j < 5;j++){
                        label[j] = 3;
                }
                 for(j = len; j != len-5;j--){
                        label[j] = 3;
                }
        }
        free_2d((void**) emission);
        free_2d((void**) transition );
        MFREE(tmp);
        return OK;
ERROR:
        if(emission){
                free_2d((void**) emission);
        }

        if(transition){
                free_2d((void**) transition);
        }
        if(tmp){
                MFREE(tmp);
        }
        return FAIL;
}


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
        sb->L = -1;
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
                        /*switch(seq[i][j]){
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
                                }*/
                        sequence->seq[c] = seq[i][j];
                        sequence->u[c] = 1.0f;
                        sequence->label[c] = 2;
                        c++;
                        if(c == sequence->malloc_len){
                                RUN(realloc_ihmm_seq(sequence));
                        }
                        
                }
                sequence->seq[len] = 0;
                
                
                sequence->seq_len = c;
                
                if(c > sb->max_len){
                        sb->max_len = c;
                }
        }
        sb->num_seq = numseq;
        RUN(detect_alphabet(sb));
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
        sb->L = -1;
               
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
                                                /*switch(line[i]){
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
                                                }*/
                                                sequence->seq[sequence->seq_len] = line[i];
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
        RUN(detect_alphabet(sb));
        return sb;
ERROR:
        free_ihmm_sequences(sb);
        if(f_ptr){
                fclose(f_ptr);
        }
        return NULL;
}

int write_ihmm_sequences(struct seq_buffer* sb, char* filename, char* comment)
{
        FILE* f_ptr = NULL;
        int i,j,c;
        int has_names;
        int max_label;
        int digits;
        int block;

        char** dwb = NULL; /* digit write buffer */
        
        
        ASSERT(sb!= NULL, "No sequences.");
        ASSERT(filename != NULL, "No filename given");


        if(sb->L == ALPHABET_DNA){
                RUN(translate_internal_to_DNA(sb));
                DPRINTF1("DNA");
        }
        if(sb->L == ALPHABET_PROTEIN){
                RUN(translate_internal_to_PROTEIN(sb));
                DPRINTF1("PROT");
        }
        DPRINTF1("L: %d",sb->L);
       
        /* Check if sequence names are present.. */
        has_names = 0;
        for (i = 0;i < sb->num_seq;i++){
                if(sb->sequences[i]->name != NULL){
                        has_names++;
                }
        }
        if(has_names != sb->num_seq){
                if(has_names == 0){
                        /* create new sequence names */
                        for (i = 0;i < sb->num_seq;i++){
                                snprintf(sb->sequences[i]->name,256,"Seq%d",i+1);
                        }
                        
                }
                if(has_names != 0 && has_names != sb->num_seq){
                        ERROR_MSG("Some sequences have names; other don't.");
                }
        }

        
        /* search for largest label  */
        max_label = 0;
        for (i = 0;i < sb->num_seq;i++){
                for(j = 0; j < sb->sequences[i]->seq_len; j++){
                        if(sb->sequences[i]->label[j] > max_label){
                                max_label = sb->sequences[i]->label[j];
                        }
                }
        }
        ASSERT(max_label != 0, "No labels found!");
        i = 10;
        digits = 1;
        while(max_label / i != 0){
                i = i * 10;
                digits++;
                if(digits == 20){
                        ERROR_MSG("more than 20 digits in sequence labels seem a bit too much.");
                }
        }
        //DPRINTF1("max_len:%d",sb->max_len );
        RUNP(dwb = malloc_2d_char(dwb, sb->max_len , 20, 0));
        f_ptr = NULL;
        /* open file and write */
        
        RUNP(f_ptr = fopen(filename, "w"));

        /* print header  */

        fprintf(f_ptr,"# Labelled sequence format:\n");
        fprintf(f_ptr,"# \n");
        fprintf(f_ptr,"# Sequences entries are divided in to block of length 70. The first\n");
        fprintf(f_ptr,"# block begins with >NAME and subsequent blocks with ^NAME (to\n");
        fprintf(f_ptr,"# avoid confusion with fastq). The first line of each block is the\n");
        fprintf(f_ptr,"# protein / nucleotide sequence. The next lines encode the\n");
        fprintf(f_ptr,"# numerical label:\n");

        fprintf(f_ptr,"# \n");
        fprintf(f_ptr,"# AC\n");
        fprintf(f_ptr,"# 01\n");
        fprintf(f_ptr,"# 12\n");
        fprintf(f_ptr,"# 03\n");
        fprintf(f_ptr,"# \n");
        fprintf(f_ptr,"# means : 'A' is labelled 010 ->  10\n");
        fprintf(f_ptr,"# means : 'C' is labelled 123 -> 123\n");
        fprintf(f_ptr,"# Comment lines begin with '#'\n");

        if(comment){
                fprintf(f_ptr,"# %s\n", comment);
        }
        fprintf(f_ptr,"# \n");
        
        
        for (i = 0;i < sb->num_seq;i++){
                for(j = 0; j < sb->sequences[i]->seq_len; j++){
                        snprintf(dwb[j],20,"%019d",sb->sequences[i]->label[j]);
                }

                for(block = 0; block <= sb->sequences[i]->seq_len / BLOCK_LEN;block++){
                        if(!block){
                                fprintf(f_ptr,">%s\n", sb->sequences[i]->name);
                        }else{
                                fprintf(f_ptr,"^%s\n", sb->sequences[i]->name);
                        }
                        for(j = block * BLOCK_LEN; j < MACRO_MIN((block +1 ) * BLOCK_LEN,sb->sequences[i]->seq_len) ; j++){
                                fprintf(f_ptr,"%c", sb->sequences[i]->seq[j]);
                        }
                        fprintf(f_ptr,"\n");
                        
                        //fprintf(f_ptr,"%s\n", sb->sequences[i]->seq);
                        for(c = digits;c > 0;c--){
                                //for(j = 0; j < sb->sequences[i]->seq_len; j++){
                                for(j = block * BLOCK_LEN; j < MACRO_MIN((block +1 ) * BLOCK_LEN,sb->sequences[i]->seq_len) ; j++){
                                        fprintf(f_ptr,"%c", dwb[j][19-c]);
                                }
                                fprintf(f_ptr,"\n");
                        }
                }
                
        }


        

        fclose(f_ptr);
        if(sb->L == ALPHABET_DNA){
                RUN(translate_DNA_to_internal(sb));                
        }
        if(sb->L == ALPHABET_PROTEIN){
                RUN(translate_PROTEIN_to_internal(sb));
        }


        free_2d((void**) dwb);
        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        if(dwb){
                free_2d((void**) dwb);
        }
        return FAIL;
}


struct seq_buffer* load_ihmm_sequences(char* in_filename)
{
        struct seq_buffer* sb  = NULL;
        struct ihmm_sequence* sequence = NULL;
        FILE* f_ptr = NULL;
        char** digit_buffer = NULL;
        char line[LINE_LEN];
        int i, seq_p;
        int label_pos;
        int old_label_pos;
        int digit;

        
        

        ASSERT(in_filename != NULL,"No input file specified - this should have been caught before!");

        RUNP(digit_buffer = malloc_2d_char(digit_buffer,BLOCK_LEN,20,0 ) );
        

        seq_p = 0;
        MMALLOC(sb,sizeof(struct seq_buffer));
        sb->malloc_num = 1024;
        sb->num_seq = -1; 
        sb->sequences = NULL;
        sb->max_len = 0;
        sb->L = -1;


        label_pos = 0;
        old_label_pos = 0;
        digit = 0;
               
        MMALLOC(sb->sequences, sizeof(struct chromosome*) *sb->malloc_num );
        for(i = 0; i < sb->malloc_num;i++){
                sb->sequences[i] = NULL;
                RUNP(sb->sequences[i] = alloc_ihmm_seq());                
        }
        
        RUNP(f_ptr = fopen(in_filename, "r" ));
        while(fgets(line, LINE_LEN, f_ptr)){
                DPRINTF1("%d (labpos: %d -> %d) %s",seq_p,  old_label_pos,label_pos, line);
                if(line[0] == '>'){

                        if(sb->num_seq != -1){ /* i.e. I read in the first sequence -> 0 */
                                //DPRINTF1("NUMSEQ: %d %d -> %d",sb->num_seq, old_label_pos,label_pos );
                                for(i = 0;i < label_pos - old_label_pos;i++){
                                        digit_buffer[i][digit] = 0;
                                        sb->sequences[sb->num_seq]->label[old_label_pos+i] = atoi(digit_buffer[i]);
                                        //DPRINTF1("translate: %s -> %d",digit_buffer[i],atoi(digit_buffer[i]));
                                }
                        }
                        
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
                        old_label_pos = 0;
                        label_pos = 0;
                        digit = 0;
                }else if(line[0] == '^'){
                        if(sb->num_seq != -1){ /* i.e. I read in the first sequence -> 0 */
                                for(i = 0;i < label_pos - old_label_pos;i++){
                                        digit_buffer[i][digit] = 0;
                                        sb->sequences[sb->num_seq]->label[old_label_pos+i] = atoi(digit_buffer[i]);
                                }
                        }
                        digit = 0;
                        seq_p = 1;
                        old_label_pos = label_pos;
                }else{	
                        if(seq_p == 1){
                                for(i = 0;i < LINE_LEN;i++){
                                        if(isalpha((int)line[i])){
                                                sequence->seq[sequence->seq_len] = line[i];
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
                                seq_p++;
                        } else if (seq_p  > 1){
                                label_pos = old_label_pos;
                                for(i = 0;i < LINE_LEN;i++){
                                        if(isdigit((int) line[i])){
                                                digit_buffer[i][digit] = line[i];
                                                label_pos++;
                                        }
                                        if(iscntrl((int)line[i])){
                                                digit++;
                                                break;  
                                        }
                                }
                                seq_p++;

                        }/* here I would look for quality values.... */
                      

                }
        }

        if(sb->num_seq){ /* i.e. I read in the first sequence -> 0 */
                //DPRINTF1("NUMSEQ: %d %d -> %d",sb->num_seq, old_label_pos,label_pos );
                for(i = 0;i < label_pos - old_label_pos;i++){
                        digit_buffer[i][digit] = 0;
                        sb->sequences[sb->num_seq]->label[old_label_pos+i] = atoi(digit_buffer[i]);
                        //DPRINTF1("translate: %s -> %d",digit_buffer[i],atoi(digit_buffer[i]));
                }
        }
        fclose(f_ptr);
        sb->num_seq++;
        RUN(detect_alphabet(sb));
        free_2d((void**) digit_buffer);
        return sb;
ERROR:
        
        free_ihmm_sequences(sb);
        if(digit_buffer){
                free_2d((void**)digit_buffer);
        }
        if(f_ptr){
                fclose(f_ptr);
        }
        return NULL;
}

int print_states_per_sequence(struct seq_buffer* sb)
{

        float* counts = NULL;
        int num_states;
        float sum;
        int i,j;
        ASSERT(sb != NULL, "No sequence buffer");
        num_states = 0;
        for(i =0; i < sb->num_seq;i++){
                for(j = 0; j < sb->sequences[i]->seq_len;j++){
                        if(sb->sequences[i]->label[j] > num_states){
                                num_states = sb->sequences[i]->label[j];
                        }
                }
        }
        ASSERT(num_states != 0,"No states found");

        num_states++;

        MMALLOC(counts, sizeof(float) * num_states);
        for(i =0; i < sb->num_seq;i++){
                sum = 0;
                for(j = 0; j < num_states;j++){
                        counts[j] = 0;
                }
                for(j = 0; j < sb->sequences[i]->seq_len;j++){
                        counts[sb->sequences[i]->label[j]]++;
                }
                fprintf(stdout,"Seq %d\t",i);
                for(j = 0; j < num_states;j++){
                        sum += counts[j];
                   
                }
                for(j = 0; j < num_states;j++){
                        fprintf(stdout,"%4.1f ",counts[j] / sum * 100.0f);
                        
                }
                fprintf(stdout,"\n");
                
        }
        MFREE(counts);
        return OK;
ERROR:
        if(counts){
             MFREE(counts);   
        }
        
        return FAIL;
}


/* The purpose is to automatically detect whether the sequences are DNA /
 * protein based on match to IUPAC codes. */
int detect_alphabet(struct seq_buffer* sb)
{
        struct ihmm_sequence* sequence = NULL;
        int i,j;
        int min,c;
        uint8_t DNA[256];
        uint8_t protein[256];
        uint8_t query[256];
        int diff[3];
        char DNA_letters[]= "acgtACGTnN";
        char protein_letters[] = "ACDEFGHIKLMNPQRSTVWY";

        
        ASSERT(sb != NULL, "No sequence buffer.");

        for(i = 0; i <256;i++){
                DNA[i] = 0;
                protein[i] = 0;
                query[i] = 0;
        }

        for(i = 0 ; i < strlen(DNA_letters);i++){
                DNA[(int) DNA_letters[i]] = 1;
        }
        
        for(i = 0 ; i < strlen(protein_letters);i++){
                protein[(int) protein_letters[i]] = 1;
        }
        
        for(i = 0; i < sb->num_seq;i++){
                sequence = sb->sequences[i];
                for(j = 0; j < sequence->seq_len;j++){
                        query[(int)sequence->seq[j]] = 1;
                }
        }
       
        diff[0] = 0;
        diff[1] = 0;
        for(i = 0; i < 256;i++){
                if(query[i] != DNA[i]){
                        diff[0]++;
                }
                if(query[i] != protein[i]){
                        diff[1]++;
                }
        }

        c = -1;
        min = 2147483647;
        for(i = 0; i < 2;i++){
                if(diff[i] < min){
                        min = diff[i];
                        c = i;
                }                
        }
        if(c == 0){
                LOG_MSG("Detected DNA sequences.");
                sb->L = ALPHABET_DNA;
                RUN(translate_DNA_to_internal(sb));
        }else if(c == 1){
                LOG_MSG("Detected protein sequences.");
                sb->L = ALPHABET_PROTEIN;
                RUN(translate_PROTEIN_to_internal(sb));
        }else{
                ERROR_MSG("Alphabet not recognized.");
        }
        return OK;
ERROR:
        return FAIL;
}

int translate_DNA_to_internal(struct seq_buffer* sb )
{
        struct ihmm_sequence* sequence = NULL;
        int i,j;

        ASSERT(sb != NULL,"No sequence buffer.");
        
        for(i = 0; i < sb->num_seq;i++){
                sequence = sb->sequences[i];
                for(j = 0; j < sequence->seq_len;j++){
                        switch(sequence->seq[j]){
                        case 'A':
                        case 'a':
                                sequence->seq[j] = 0;
                                break;
                        case 'C':
                        case 'c':
                                sequence->seq[j] = 1;
                                break;
                        case 'G':
                        case 'g':
                                sequence->seq[j] = 2;
                                break;
                        case 'T':
                        case 't':
                                sequence->seq[j] = 3;
                                break;
                        case 'N':
                        case 'n':
                                sequence->seq[j] = random_int_zero_to_x(3);
                                break;
                                
                        default:
                                ERROR_MSG("Non ACGTN letter in sequence:%d %s.",i,sequence->seq);
                                break;
                        }

                }
                sequence->seq[sequence->seq_len] = 0;
        }
        return OK;
ERROR:
        return FAIL;
}


int translate_internal_to_DNA(struct seq_buffer* sb )
{
        struct ihmm_sequence* sequence = NULL;
        int i,j;
        
        ASSERT(sb != NULL,"No sequence buffer.");
        
        for(i = 0; i < sb->num_seq;i++){
                sequence = sb->sequences[i];
                for(j = 0; j < sequence->seq_len;j++){
                        sequence->seq[j] = "ACGT"[sequence->seq[j]];
                }
                sequence->seq[sequence->seq_len] = 0;
        }
        return OK;
ERROR:
        return FAIL;
}



int translate_PROTEIN_to_internal(struct seq_buffer* sb )
{
        struct ihmm_sequence* sequence = NULL;
        int i,j;

        ASSERT(sb != NULL,"No sequence buffer.");
        
        for(i = 0; i < sb->num_seq;i++){
                sequence = sb->sequences[i];
                for(j = 0; j < sequence->seq_len;j++){
                        switch(sequence->seq[j]){

                                
                        case 'A':
                                sequence->seq[j] = 0;
                                break;
                        case 'C':
                                sequence->seq[j] = 1;
                                break;
                        case 'D':
                                sequence->seq[j] = 2;
                                break;
                        case 'E':
                                sequence->seq[j] = 3;
                                break;
                        case 'F':
                                sequence->seq[j] = 4;
                                break;
                        case 'G':
                                sequence->seq[j] = 5;
                                break;
                        case 'H':
                                sequence->seq[j] = 6;
                                break;
                        case 'I':
                                sequence->seq[j] = 7;
                                break;
                        case 'K':
                                sequence->seq[j] = 8;
                                break;
                        case 'L':
                                sequence->seq[j] = 9;
                                break;
                        case 'M':
                                sequence->seq[j] = 10;
                                break;
                        case 'N':
                                sequence->seq[j] = 11;
                                break;
                                // ACDEFGHIKLMNPQRSTVWY"
                        case 'P':
                                sequence->seq[j] = 12;
                                break;
                        case 'Q':
                                sequence->seq[j] = 13;
                                break;
                        case 'R':
                                sequence->seq[j] = 14;
                                break;
                        case 'S':
                                sequence->seq[j] = 15;
                                break;
                        case 'T':
                                sequence->seq[j] = 16;
                                break;
                        case 'V':
                                sequence->seq[j] = 17;
                                break;
                        case 'W':
                                sequence->seq[j] = 18;
                                break;
                        case 'Y':
                                sequence->seq[j] = 19;
                                break;
                        default:
                                ERROR_MSG("Non ACGTN letter in sequence:%d %s.",i,sequence->seq);
                                break;
                        }

                }
        }
        return OK;
ERROR:
        return FAIL;
}


int translate_internal_to_PROTEIN(struct seq_buffer* sb )
{
        struct ihmm_sequence* sequence = NULL;
        int i,j;

        ASSERT(sb != NULL,"No sequence buffer.");
        for(i = 0; i < sb->num_seq;i++){
                sequence = sb->sequences[i];
                for(j = 0; j < sequence->seq_len;j++){
                        sequence->seq[j] = "ACDEFGHIKLMNPQRSTVWY"[sequence->seq[j]];
                }
        }
        return OK;
ERROR:
        return FAIL;
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

struct ihmm_sequence* add_spacer_ihmm_seq(struct ihmm_sequence* sequence, int space_len, int L)
{
        int i;
        
        struct ihmm_sequence* tmp = NULL;
        ASSERT(sequence != NULL, "No Sequence.");
        RUNP(tmp = alloc_ihmm_seq());
        if(sequence->seq_len+10 == sequence->malloc_len){
                RUN(realloc_ihmm_seq(tmp));
        }
        /* First X spacer letters  */
        for(i = 0; i < space_len;i++){
                tmp->seq[sequence->seq_len] = L;
                tmp->u[sequence->seq_len] = 1.0f;
                tmp->label[sequence->seq_len] = 2;
                tmp->seq_len++;    

        }
        /* copy sequence  */
        for(i = 0 ; i < sequence->seq_len;i++){
                tmp->seq[sequence->seq_len] = sequence->seq[i]; 
                tmp->u[sequence->seq_len] = 1.0f;
                tmp->label[sequence->seq_len] = 2;
                tmp->seq_len++;    
        }
        /* Last X spacer letters  */
        for(i = 0; i < space_len;i++){
 
                tmp->seq[sequence->seq_len] = L;
                tmp->u[sequence->seq_len] = 1.0f;
                tmp->label[sequence->seq_len] = 2;
                tmp->seq_len++;    

        }
        free_ihmm_sequence(sequence);
        
        return tmp;
ERROR:
        if(tmp){
                free_ihmm_sequence(tmp);
        }
        return NULL;
}

int print_labelled_ihmm_buffer(struct seq_buffer* sb)
{
        struct ihmm_sequence* sequence = NULL;
        int i,j;
        ASSERT(sb != NULL, "No sequence buffer");
        if(sb->L == ALPHABET_DNA){
                RUN(translate_internal_to_DNA(sb));                
        }
        if(sb->L == ALPHABET_PROTEIN){
                RUN(translate_internal_to_PROTEIN(sb));
        }

        for(i = 0; i < sb->num_seq;i++){
                sequence = sb->sequences[i];
                for(j = 0; j < sequence->seq_len;j++){
                        fprintf(stdout,"%3c",sequence->seq[j]);
                }
                fprintf(stdout,"\n");
                for(j = 0; j < sequence->seq_len;j++){
                        fprintf(stdout,"%3d",sequence->label[j]);
                }
                fprintf(stdout,"\n");
        }
        if(sb->L == ALPHABET_DNA){
                RUN(translate_DNA_to_internal(sb));                
        }
        if(sb->L == ALPHABET_PROTEIN){
                RUN(translate_PROTEIN_to_internal(sb));
        }
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
        if(sb){
                for(i =0; i < sb->malloc_num;i++){
                        free_ihmm_sequence(sb->sequences[i]);
                }
                MFREE(sb->sequences);
                MFREE(sb);
        }
}


int compare_sequence_buffers(struct seq_buffer* a, struct seq_buffer* b)
{
        int i,j;


        ASSERT(a != NULL, "No a");
        ASSERT(b != NULL, "No b");

        ASSERT(a->max_len == b->max_len , "max len differs");
        ASSERT(a->num_seq == b->num_seq, "number of sequences differ");

        for(i = 0; i < a->num_seq;i++){
                ASSERT(strcmp(a->sequences[i]->name,b->sequences[i]->name) == 0 , "Names differ");
                ASSERT(a->sequences[i]->seq_len == b->sequences[i]->seq_len, "Sequence lengths differ" );
                for(j = 0; j < a->sequences[i]->seq_len;j++){
                       ASSERT(a->sequences[i]->seq[j]  == b->sequences[i]->seq[j], "Sequences  differ" ); 
                }

                for(j = 0; j < a->sequences[i]->seq_len;j++){
                       ASSERT(a->sequences[i]->label[j]  == b->sequences[i]->label[j], "Labels differ" ); 
                }
        }
        
        return OK;
ERROR:
        return FAIL;
}


#ifdef ITESTSEQ

int main(const int argc,const char * argv[])
{
        struct seq_buffer* iseq = NULL;
        struct seq_buffer* iseq_b = NULL;
        
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
        RUN(write_ihmm_sequences(iseq,"test.lfa","generated by iseq_ITEST"));

      

        RUNP(iseq_b  =load_ihmm_sequences("test.lfa"));

        RUN(compare_sequence_buffers(iseq,iseq_b));
        free_ihmm_sequences(iseq);
        free_ihmm_sequences(iseq_b);




        //Protein test...
        char *tmp_seq2[18] = {
"RRRAHTQAEQKRRDAIKRGYDDLQTIVPTCQQQDFSIGSQKLSKAIVLQKTIDYIQFLH",
"RREAHTQAEQKRRDAIKKGYDSLQELVPRCQPNDSSGYKLSKALILQKSIEYIGYL",
"RRITHISAEQKRRFNIKLGFDTLHGLVSTLSAQPSLKVSKATTLQKTAEYILMLQ",
"RRAGHIHAEQKRRYNIKNGFDTLHALIPQLQQNPNAKLSKAAMLQKGADHIKQLR",
"KRILHLHAEQNRRSALKDGFDQLMDIIPDLYSGGVKPTNAVVLAKSADHIRRLQ",
"KKATHLRCERQRREAINSGYSDLKDLIPQTTTSLGCKTTNAAILFRACDFMSQLK",
"LRTSHKLAERKRRKEIKELFDDLKDALPLDKSTKSSKWGLLTRAIQYIEQLK",
"YRRTHTANERRRRGEMRDLFEKLKITLGLLHSSKVSKSLILTRAFSEIQGLT",
"TRKSVSERKRRDEINELLENLKTIVQNPSDSNEKISHETILFRVFERVSGVD",
"GHRSETEKQRRDDTNDLLNEFKKIVQKSESEKLSKEEVLFRIVKLLSGIQ",
"KRAHHNALERKRRDHIKDSFHSLRDSVPSLQGEKASRAQILDKATEYIQYMR",
"RRAHHNELERRRRDHIKDHFTILKDAIPLLDGEKSSRALILKRAVEFIHVMQ",
"KRAHHNALERRRRDHIKESFTNLREAVPTLKGEKASRAQILKKTTECIQTMR",
"GRHVHNELEKRRRAQLKRCLEQLRQQMPLGVDHTRYTTLSLLRGARMHIQKLE",
"NRSSHNELEKHRRAKLRLYLEQLKQLVPLGPDSTRHTTLSLLKRAKVHIKKLE",
"SRSTHNEMEKNRRAHLRLCLEKLKGLVPLGPESSRHTTLSLLTKAKLHIKKLE",
"NRTSHNELEKNRRAHLRNCLDGLKAIVPLNQDATRHTTLGLLTQARALIENLK",
"NRSTHNELEKNRRAHLRLCLERLKVLIPLGPDCTRHTTLGLLNKAKAHIKKLE"};
        RUNP(iseq = create_ihmm_sequences_mem(tmp_seq2 ,18));

        RUN(random_label_ihmm_sequences(iseq, 123));
        RUN(print_labelled_ihmm_buffer(iseq));

        RUN(write_ihmm_sequences(iseq,"test.lfa","generated by iseq_ITEST"));

      

        RUNP(iseq_b  =load_ihmm_sequences("test.lfa"));

        RUN(compare_sequence_buffers(iseq,iseq_b));
        
        RUN(print_labelled_ihmm_buffer(iseq));
        free_ihmm_sequences(iseq);
        free_ihmm_sequences(iseq_b);

        RUNP(iseq = create_ihmm_sequences_mem(tmp_seq ,18));
        LOG_MSG("alpha:100");
        RUN(dirichlet_emission_label_ihmm_sequences(iseq, 2, 100));
        RUN(print_labelled_ihmm_buffer(iseq));

        

        LOG_MSG("alpha: 0.3");
        RUN(dirichlet_emission_label_ihmm_sequences(iseq, 2, 0.3));
        RUN(print_labelled_ihmm_buffer(iseq));

        LOG_MSG("alpha: 0.05");
        RUN(dirichlet_emission_label_ihmm_sequences(iseq, 2, 0.05));
        RUN(print_labelled_ihmm_buffer(iseq));
        


        free_ihmm_sequences(iseq);
                
        return EXIT_SUCCESS;
ERROR:
        free_ihmm_sequences(iseq);
        free_ihmm_sequences(iseq_b);
        return EXIT_FAILURE;
        
}
#endif
