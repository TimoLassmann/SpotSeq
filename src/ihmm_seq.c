#include "ihmm_seq.h"

static void free_ihmm_sequence(struct ihmm_sequence* sequence);


int detect_alphabet(struct seq_buffer* sb);

int translate_DNA_to_internal(struct seq_buffer* sb);
int translate_internal_to_DNA(struct seq_buffer* sb);

int translate_internal_to_PROTEIN(struct seq_buffer* sb);
int translate_PROTEIN_to_internal(struct seq_buffer* sb);

int compare_sequence_buffers(struct seq_buffer* a, struct seq_buffer* b);

int add_background_to_sequence_buffer(struct seq_buffer* sb);

int reverse_complement(struct ihmm_sequence*  org,struct ihmm_sequence* dest);

struct seq_buffer* concatenate_sequences(struct seq_buffer* sb)
{
        struct seq_buffer* sb_concat = NULL;
        struct ihmm_sequence* s = NULL;
        struct ihmm_sequence* a = NULL;
        int i,j;
        int new_len = 0;

        ASSERT(sb!= NULL, "No sequence buffer");
        for(i = 0; i < sb->num_seq;i++){
                new_len += sb->sequences[i]->seq_len;

        }

        MMALLOC(sb_concat,sizeof(struct seq_buffer));
        sb_concat->malloc_num = 1 ;
        sb_concat->num_seq = 1;
        sb_concat->sequences = NULL;
        sb_concat->max_len = new_len;
        sb_concat->L = sb->L;

        sb_concat->background = NULL;

        MMALLOC(sb_concat->background, sizeof(float)*  sb->L);
        for(i = 0;i < sb->L;i++){
                sb_concat->background[i] = sb->background[i];
        }


        MMALLOC(sb_concat->sequences, sizeof(struct ihmm_sequence*) *sb_concat->malloc_num );
        for(i = 0; i < sb_concat->num_seq;i++){
                sb_concat->sequences[i] = NULL;
                RUNP(sb_concat->sequences[i] = alloc_ihmm_seq());
                while(sb_concat->sequences[i]->malloc_len <= new_len){
                        RUN(realloc_ihmm_seq(sb_concat->sequences[i]));
                }
        }
        a = sb_concat->sequences[0];
        for(i = 0; i < sb->num_seq;i++){
                s = sb->sequences[i];
                for(j = 0; j < s->seq_len;j++){
                        a->seq[a->seq_len] = s->seq[j];
                        a->seq_len++;
                }

        }
        free_ihmm_sequences(sb);
        return sb_concat;
ERROR:
        return NULL;
}


struct seq_buffer* get_sequences_from_hdf5_model(char* filename)
{
        struct hdf5_data* hdf5_data = NULL;
        struct seq_buffer* sb = NULL;
        char** name = NULL;
        char** seq = NULL;
        int** label = NULL;
        float* background;

        int num_seq;
        int max_len;
        int max_name_len;
        int local_L;
        int i,j;

        ASSERT(filename != NULL, "No filename");
        ASSERT(my_file_exists(filename) != 0,"File %s does not exist.",filename);

        hdf5_data = hdf5_create();

        hdf5_open_file(filename,hdf5_data);
        hdf5_read_attributes(hdf5_data,hdf5_data->file);
        //print_attributes(hdf5_data);
        get_group_names(hdf5_data);

        hdf5_open_group("SequenceInformation",hdf5_data);
        hdf5_read_attributes(hdf5_data, hdf5_data->group);
        ASSERT(hdf5_data->num_attr != 0 , "Could not find attributes");
        //print_attributes(hdf5_data);

        num_seq = -1;
        max_len = -1;
        max_name_len = -1;
        local_L = -1;
        for(i = 0; i < hdf5_data->num_attr;i++){
                if(!strncmp("NumSeq", hdf5_data->attr[i]->attr_name, 6)){
                        num_seq = hdf5_data->attr[i]->int_val;
                }

                if(!strncmp("MaxLen", hdf5_data->attr[i]->attr_name, 6)){
                        max_len = hdf5_data->attr[i]->int_val;
                }

                if(!strncmp("MaxNameLen", hdf5_data->attr[i]->attr_name, 10)){
                        max_name_len = hdf5_data->attr[i]->int_val;
                }

                if(!strncmp("Alphabet", hdf5_data->attr[i]->attr_name, 8)){
                        local_L = hdf5_data->attr[i]->int_val;
                }
        }

        ASSERT(num_seq != -1, "No numseq");
        ASSERT(max_len != -1, "No maxlen");
        ASSERT(max_name_len != -1, "No maxnamelen");
        ASSERT(local_L != -1,"No Alphabet");



        hdf5_read_dataset("Names",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read beta");
        name = (char**)hdf5_data->data;


        hdf5_read_dataset("Sequences",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read transition_counts");
        seq = (char**)hdf5_data->data;


        hdf5_read_dataset("Background",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 1, "Could not read background");
        background = (float*)hdf5_data->data;

        hdf5_read_dataset("Labels",hdf5_data);
        ASSERT(hdf5_data->data != NULL && hdf5_data->rank == 2, "Could not read emission_counts");
        label = (int**)hdf5_data->data;
        hdf5_close_group(hdf5_data);
        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);


        MMALLOC(sb,sizeof(struct seq_buffer));
        sb->malloc_num = num_seq ;
        sb->num_seq = num_seq;
        sb->sequences = NULL;
        sb->max_len = max_len;
        sb->L = local_L;

        sb->background = background;


        MMALLOC(sb->sequences, sizeof(struct ihmm_sequence*) *sb->malloc_num );
        for(i = 0; i < sb->num_seq;i++){
                sb->sequences[i] = NULL;
                RUNP(sb->sequences[i] = alloc_ihmm_seq());
                while(sb->sequences[i]->malloc_len <= sb->max_len){
                        RUN(realloc_ihmm_seq(sb->sequences[i]));
                }
        }

        /* copy stuff over */
        for(i = 0; i < sb->num_seq;i++){
                for (j = 0; j < max_name_len;j++){
                        if(!name[i][j]){
                                sb->sequences[i]->name[j] = 0;
                                break;
                        }
                        sb->sequences[i]->name[j] = name[i][j];
                }
        }

        /* make sequence matrix */
        for(i = 0; i < sb->num_seq;i++){
                for (j = 0; j < sb->max_len;j++){
                        if(seq[i][j] == -1){
                                break;
                        }
                        sb->sequences[i]->seq[j] = seq[i][j];
                        sb->sequences[i]->label[j] = label[i][j];
                }
                sb->sequences[i]->seq_len = j;
        }

        free_2d((void**)   label);
        free_2d((void**)   name);
        free_2d((void**)   seq);

        return sb;
ERROR:
        if(hdf5_data){
                hdf5_close_file(hdf5_data);
                hdf5_free(hdf5_data);
        }
        if(label){
                free_2d((void**)   label);
        }
        if(name){
                free_2d((void**)   name);
        }
        if(seq){
                free_2d((void**)   seq);
        }


        return NULL;
}

/* I want to save the labelled sequences in the hdf5model file to be able to resume jobs.  */
int add_sequences_to_hdf5_model(char* filename,struct seq_buffer* sb)
{

        struct hdf5_data* hdf5_data = NULL;
        int i,j,len;

        char** name = NULL;
        char** seq = NULL;
        int** label = NULL;

        int max_name_len;
        int has_seq_info;


        ASSERT(sb!=NULL, "No sequence buffer");

        /* make sequence name matrix */
        max_name_len = -1;
        for(i = 0; i < sb->num_seq;i++){
                len = strlen(sb->sequences[i]->name);
                if(len > max_name_len){
                        max_name_len = len;
                }
        }
        max_name_len+=1;

        RUNP(name = malloc_2d_char(name, sb->num_seq, max_name_len, 0));
        for(i = 0; i < sb->num_seq;i++){
                len = strlen(sb->sequences[i]->name);
                for (j = 0; j < len;j++){
                        name[i][j] = sb->sequences[i]->name[j];
                }
        }

        /* make sequence matrix */

        RUNP(seq = malloc_2d_char(seq, sb->num_seq, sb->max_len, -1));
        for(i = 0; i < sb->num_seq;i++){
                len = sb->sequences[i]->seq_len;
                for (j = 0; j < len;j++){
                        seq[i][j] = sb->sequences[i]->seq[j];
                }
        }

        /* make  label matrix */

        label = malloc_2d_int(label, sb->num_seq, sb->max_len, -1);
        for(i = 0; i < sb->num_seq;i++){
                len = sb->sequences[i]->seq_len;
                for (j = 0; j < len;j++){
                        label[i][j] = sb->sequences[i]->label[j];
                }
        }

        // RUNP(name = alloc_2)

        RUNP(hdf5_data = hdf5_create());

        hdf5_open_file(filename,hdf5_data);

        get_group_names(hdf5_data);

        has_seq_info = 0;
        fprintf(stdout,"Groups:\n");
        for(i = 0; i < hdf5_data->grp_names->num_names;i++){
                if(strncmp("SequenceInformation", hdf5_data->grp_names->names[i],19) == 0){
                        has_seq_info = 1;
                }
                fprintf(stdout,"%d %s\n",i,hdf5_data->grp_names->names[i]);
        }
        if(!has_seq_info ){
                hdf5_create_group("SequenceInformation",hdf5_data);
        }else{
                hdf5_open_group("SequenceInformation",hdf5_data);
        }

        hdf5_data->num_attr = 0;

        hdf5_add_attribute(hdf5_data, "NumSeq",     "", sb->num_seq , 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "MaxLen",     "", sb->max_len , 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "MaxNameLen", "", max_name_len, 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "Alphabet",   "", sb->L       , 0.0f, HDF5GLUE_INT);

        hdf5_write_attributes(hdf5_data, hdf5_data->group);

        hdf5_data->rank = 2;
        hdf5_data->dim[0] = sb->num_seq;
        hdf5_data->dim[1] = max_name_len;
        hdf5_data->chunk_dim[0] = sb->num_seq;
        hdf5_data->chunk_dim[1] =max_name_len;
        hdf5_data->native_type = H5T_NATIVE_CHAR;
        hdf5_write("Names",&name[0][0], hdf5_data);

        hdf5_data->rank = 2;
        hdf5_data->dim[0] = sb->num_seq;
        hdf5_data->dim[1] = sb->max_len;
        hdf5_data->chunk_dim[0] = sb->num_seq;
        hdf5_data->chunk_dim[1] = sb->max_len;
        hdf5_data->native_type = H5T_NATIVE_CHAR;
        hdf5_write("Sequences",&seq[0][0], hdf5_data);


        hdf5_data->rank = 2;
        hdf5_data->dim[0] = sb->num_seq;
        hdf5_data->dim[1] = sb->max_len;
        hdf5_data->chunk_dim[0] = sb->num_seq;
        hdf5_data->chunk_dim[1] = sb->max_len;
        hdf5_data->native_type = H5T_NATIVE_INT;
        hdf5_write("Labels",&label[0][0], hdf5_data);




        hdf5_data->rank = 1;
        hdf5_data->dim[0] = sb->L;
        hdf5_data->dim[1] = 0;
        hdf5_data->chunk_dim[0] = sb->L;
        hdf5_data->chunk_dim[1] = 0;
        hdf5_data->native_type = H5T_NATIVE_FLOAT;
        hdf5_write("Background",&sb->background[0], hdf5_data);

        hdf5_close_group(hdf5_data);
        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);

        free_2d((void**)   label);
        free_2d((void**)   name);
        free_2d((void**)   seq);

        return OK;
ERROR:
        if(hdf5_data){
                hdf5_close_file(hdf5_data);
                hdf5_free(hdf5_data);
        }
        if(label){
                free_2d((void**)   label);
        }
        if(name){
                free_2d((void**)   name);
        }
        if(seq){
                free_2d((void**)   seq);
        }
        return FAIL;
}


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
                        emission[i][j] = rk_gamma(&rndstate,alpha + 1.0 + ((float)i) / (float)(k)*10.0f, 1.0);
                        sum += emission[i][j];
                }
                sanity = 0.0f;
                for(j = 0;j < n;j++){
                        emission[i][j] /= sum;
                        emission[k][j] += emission[i][j]; /* Last row has sums of all emissions of Letter 0, 1, 2, .. L  *\/ */
                        sanity += emission[i][j];
                        //          fprintf(stdout,"%f ",emission[i][j]);
                }
                //       fprintf(stdout,"sum: %f\n",sanity);
        }

        //fprintf(stdout,"Transition\n");
        for(i = 0; i < k;i++){
                sum = 0.0;
                for(j = 0; j < k;j++){

                        transition[i][j] = rk_gamma(&rndstate,alpha , 1.0);
                        //if(i == j){
                        //        transition[i][j] = rk_gamma(&rndstate,alpha+((float)j) / (float) k , 1.0);
                        //}
                        sum += transition[i][j];
                }
                sanity = 0.0f;
                for(j = 0; j < k;j++){
                        transition[i][j] /= sum;
                        transition[k][j] += transition[i][j];
                        sanity += transition[i][j];
                        //              fprintf(stdout,"%f ",transition[i][j]);
                }
                //      fprintf(stdout,"sum: %f\n", sanity);
        }
        /* for(i = 0; i < 1000;i++){ */
        /*         fprintf(stdout,"%d %d\n", i, rk_interval(10,&rndstate)); */
        /* } */
        //exit(0);
        cur_state = -1;
        for(i = 0;i< sb->num_seq;i++){

                label = sb->sequences[i]->label;
                len = sb->sequences[i]->seq_len;
                seq = sb->sequences[i]->seq;
                r = rk_double(&rndstate) * emission[k][seq[0]];//  random_float_zero_to_x(emission[k][seq[0]]);
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
                        r = rk_double(&rndstate) *  sum;//random_float_zero_to_x(sum);
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
        /*for(i = 0;i< sb->num_seq;i++){
                label = sb->sequences[i]->label;
                len = sb->sequences[i]->seq_len;
                seq = sb->sequences[i]->seq;
                for(j = 0; j < 5;j++){
                        label[j] = 3;
                }
                 for(j = len; j != len-5;j--){
                        label[j] = 3;
                }
                }*/
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


int random_label_ihmm_sequences(struct seq_buffer* sb, int k,float alpha)
{
        int* label;
        int i,j,c;
        int len;
        float* state_prob = NULL;
        float sum,r;
        rk_state rndstate;

        ASSERT(alpha > 0.0,"alpha should not be zero");
        //RUN(set_random_seed(void));

        ASSERT(sb != NULL,"No sequences");
        rk_randomseed(&rndstate);
        MMALLOC(state_prob, sizeof(float) * k);

        sum = 0.0;
        for(i = 0; i< k; i++){
                state_prob[i] = rk_gamma(&rndstate, alpha, 1.0);
                sum += state_prob[i];
        }
        for(i = 0; i< k; i++){
                state_prob[i] = state_prob[i] / sum;
        }


        for(i = 0;i< sb->num_seq;i++){

                label = sb->sequences[i]->label;
                len = sb->sequences[i]->seq_len;
                DPRINTF3("Seq%d len %d ",i ,len);
                for(j = 0;j < len;j++){
                        r = rk_double(&rndstate);
                        sum = 0;
                        for(c = 0; c < k;c++){
                                sum += state_prob[c];
                                if(r < sum){
                                        label[j] = c+2;
                                        break;
                                }
                        }
                        label[j] =  random_int_zero_to_x(k-1) + 2;
                        DPRINTF3("%d",label[j]);
                }
        }
        MFREE(state_prob);
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
        sb->background = NULL;
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
        RUN(add_background_to_sequence_buffer(sb));
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
        sb->background = NULL;

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

        RUN(add_background_to_sequence_buffer(sb));
        return sb;
ERROR:
        free_ihmm_sequences(sb);
        if(f_ptr){
                fclose(f_ptr);
        }
        return NULL;
}

int add_background_to_sequence_buffer(struct seq_buffer* sb)
{
        int i,j;
        float sum = 0.0f;
        ASSERT(sb!= NULL, "No sequence buffer");
        ASSERT(sb->L != 0, "No alphabet detected");

        MMALLOC(sb->background, sizeof(float)* sb->L);
        for(i = 0; i < sb->L;i++){
               sb->background[i] = 0.0f;
        }

        for(i = 0; i < sb->num_seq;i++){
                for(j = 0;j < sb->sequences[i]->seq_len;j++){
                        sb->background[sb->sequences[i]->seq[j]]++;
                        sum++;
                }
        }
        ASSERT(sum != 0.0f,"No sequence counts found");
        for(i = 0; i < sb->L;i++){
                //LOG_MSG("%d %f %f",i,sb->background[i], sum);
                 sb->background[i] /= sum;
        }
        // exit(0);
        return OK;
ERROR:
        return FAIL;
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
        DPRINTF1("numseq: %d",sb->num_seq );

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
                        fprintf(stdout,"%d %d %d\n",i,j,sb->sequences[i]->label[j]);
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
int write_ihmm_sequences_fasta(struct seq_buffer* sb, char* filename)
{
        FILE* f_ptr = NULL;
        int i,j,c;
        int has_names;
        int max_label;
        int digits;
        int block;


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
        DPRINTF1("numseq: %d",sb->num_seq );

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



        //DPRINTF1("max_len:%d",sb->max_len );

        f_ptr = NULL;
        /* open file and write */

        RUNP(f_ptr = fopen(filename, "w"));


        for (i = 0;i < sb->num_seq;i++){


                for(block = 0; block <= sb->sequences[i]->seq_len / BLOCK_LEN;block++){
                        if(!block){
                                fprintf(f_ptr,">%s\n", sb->sequences[i]->name);
                        }//else{
                         //       fprintf(f_ptr,"^%s\n", sb->sequences[i]->name);
                        //}
                        for(j = block * BLOCK_LEN; j < MACRO_MIN((block +1 ) * BLOCK_LEN,sb->sequences[i]->seq_len) ; j++){
                                fprintf(f_ptr,"%c", sb->sequences[i]->seq[j]);
                        }
                        fprintf(f_ptr,"\n");

                        //fprintf(f_ptr,"%s\n", sb->sequences[i]->seq);
                        //for(c = digits;c > 0;c--){
                                //for(j = 0; j < sb->sequences[i]->seq_len; j++){
                        //        for(j = block * BLOCK_LEN; j < MACRO_MIN((block +1 ) * BLOCK_LEN,sb->sequences[i]->seq_len) ; j++){
                        //               fprintf(f_ptr,"%c", dwb[j][19-c]);
                        //       }
                        //       fprintf(f_ptr,"\n");
                        //}
                }

        }




        fclose(f_ptr);
        if(sb->L == ALPHABET_DNA){
                RUN(translate_DNA_to_internal(sb));
        }
        if(sb->L == ALPHABET_PROTEIN){
                RUN(translate_PROTEIN_to_internal(sb));
        }
        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
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
        sb->background = NULL;

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
        RUN(add_background_to_sequence_buffer(sb));
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

int add_reverse_complement_sequences_to_buffer(struct seq_buffer* sb)
{
        struct ihmm_sequence* sequence = NULL;
        char buffer[BUFFER_LEN];
        int i,j;
        int old_numseq;
        ASSERT(sb!= NULL, "No sequence buffer");
        ASSERT(sb->L == ALPHABET_DNA, "No DNA sequences in buffer");
        /* remember old sequence count. */
        old_numseq = sb->num_seq;
        /* main loop */
        for(i = 0; i < old_numseq;i++){
                sequence = sb->sequences[sb->num_seq];
                /* alloc space for revcomp sequence */
                while(sequence->malloc_len < sb->sequences[i]->seq_len){
                        RUN(realloc_ihmm_seq(sequence));
                }
                /* copy sequence name */
                snprintf(sequence->name, 256 , "%s_rev", sb->sequences[i]->name);
                /* copy & reverse_complement sequences */
                reverse_complement(sb->sequences[i], sequence);
                sequence->seq_len = sb->sequences[i]->seq_len;
                sb->num_seq++;
                if(sb->num_seq == sb->malloc_num){
                        sb->malloc_num = sb->malloc_num << 1;
                        MREALLOC(sb->sequences,sizeof(struct ihmm_sequence*) * sb->malloc_num);
                        for(i = sb->num_seq; i < sb->malloc_num;i++){
                                sb->sequences[i] = NULL;
                                RUNP(sb->sequences[i] = alloc_ihmm_seq());
                        }
                }
        }

        return OK;
ERROR:
        return FAIL;
}

int reverse_complement(struct ihmm_sequence*  org,struct ihmm_sequence* dest)
{
        int rev[5] = {3,2,1,0,4};

        int i,c;
        int len;
        uint8_t* org_s;
        uint8_t* dest_s;
        len = org->seq_len;
        org_s = org->seq;
        dest_s = dest->seq;
        c = 0;
        for(i = len-1; i >= 0;i--){
                dest_s[c] = rev[org_s[i]];
                dest->label[c] = org->label[c];
                dest->u[c] =   org->u[c];
                c++;
        }
        return OK;
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
        int r;
        ASSERT(sb != NULL,"No sequence buffer.");

        /* A 	Adenine */
        /* C 	Cytosine */
        /* G 	Guanine */
        /* T (or U) 	Thymine (or Uracil) */
        /* R 	A or G */
        /* Y 	C or T */
        /* S 	G or C */
        /* W 	A or T */
        /* K 	G or T */
        /* M 	A or C */
        /* B 	C or G or T */
        /* D 	A or G or T */
        /* H 	A or C or T */
        /* V 	A or C or G */
        /* N 	any base */
        /* . or - 	gap */

        for(i = 0; i < sb->num_seq;i++){
                sequence = sb->sequences[i];
                for(j = 0; j < sequence->seq_len;j++){
                        switch(sequence->seq[j]){
                        case 'A':
                        case 'a':
                                /* A 	Adenine */
                                sequence->seq[j] = 0;
                                break;
                        case 'C':
                        case 'c':
                                /* C 	Cytosine */
                                sequence->seq[j] = 1;
                                break;
                        case 'G':
                        case 'g':
                                /* G 	Guanine */
                                sequence->seq[j] = 2;
                                break;
                        case 'T':
                        case 't':
                        case 'U':
                        case 'u':
                                /* T (or U) 	Thymine (or Uracil) */
                                sequence->seq[j] = 3;
                                break;
                        case 'R':
                        case 'r':
                                /* R 	A or G */
                                r = random_int_zero_to_x(1);
                                if(r == 0){
                                        sequence->seq[j] = 0;
                                }else{
                                        sequence->seq[j] = 2;
                                }
                                break;
                        case 'Y':
                        case 'y':
                                /* Y 	C or T */
                                r = random_int_zero_to_x(1);
                                if(r == 0){
                                        sequence->seq[j] = 1;
                                }else{
                                        sequence->seq[j] = 3;
                                }
                                break;
                        case 'S':
                        case 's':
                                /* S 	G or C */
                                r = random_int_zero_to_x(1);
                                if(r == 0){
                                        sequence->seq[j] = 1;
                                }else{
                                        sequence->seq[j] = 2;
                                }
                                break;
                        case 'W':
                        case 'w':
                                /* W 	A or T */
                                r = random_int_zero_to_x(1);
                                if(r == 0){
                                        sequence->seq[j] = 0;
                                }else{
                                        sequence->seq[j] = 3;
                                }
                                break;
                        case 'K':
                        case 'k':
                                /* K 	G or T */
                                r = random_int_zero_to_x(1);
                                if(r == 0){
                                        sequence->seq[j] = 2;
                                }else{
                                        sequence->seq[j] = 3;
                                }
                                break;
                        case 'M':
                        case 'm':
                                /* M 	A or C */
                                r = random_int_zero_to_x(1);
                                if(r == 0){
                                        sequence->seq[j] = 0;
                                }else{
                                        sequence->seq[j] = 2;
                                }
                                break;
                        case 'B':
                        case 'b':
                                /* B 	C or G or T */
                                r = random_int_zero_to_x(2);
                                if(r == 0){
                                        sequence->seq[j] = 1;
                                }else if(r == 1){
                                        sequence->seq[j] = 2;
                                }else{
                                        sequence->seq[j] = 3;
                                }
                                break;
                        case 'D':
                        case 'd':
                                /* D 	A or G or T */
                                r = random_int_zero_to_x(2);
                                if(r == 0){
                                        sequence->seq[j] = 0;
                                }else if(r == 1){
                                        sequence->seq[j] = 2;
                                }else{
                                        sequence->seq[j] = 3;
                                }
                                break;
                        case 'H':
                        case 'h':
                                /* H 	A or C or T */
                                r = random_int_zero_to_x(2);
                                if(r == 0){
                                        sequence->seq[j] = 0;
                                }else if(r == 1){
                                        sequence->seq[j] = 1;
                                }else{
                                        sequence->seq[j] = 3;
                                }
                                break;
                        case 'V':
                        case 'v':
                                /* V 	A or C or G */
                                r = random_int_zero_to_x(2);
                                if(r == 0){
                                        sequence->seq[j] = 0;
                                }else if(r == 1){
                                        sequence->seq[j] = 1;
                                }else{
                                        sequence->seq[j] = 2;
                                }
                                break;
                        case 'N':
                        case 'n':
                                sequence->seq[j] = random_int_zero_to_x(3);
                                break;
                        default:
                                ERROR_MSG("Non ACGTN letter in sequence:%s (%d) %s (%d out of %d).",sequence->name, i,sequence->seq,j, sequence->seq_len);
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
        int r;

        ASSERT(sb != NULL,"No sequence buffer.");

        for(i = 0; i < sb->num_seq;i++){
                sequence = sb->sequences[i];
                for(j = 0; j < sequence->seq_len;j++){
                        /* Amino Acid Code: Three letter Code: Amino Acid: */
                        /* ---------------- ------------------ ----------- */
                        /* A.................Ala.................Alanine */
                        /* B.................Asx.................Aspartic acid or Asparagine */
                        /* C.................Cys.................Cysteine */
                        /* D.................Asp.................Aspartic Acid */
                        /* E.................Glu.................Glutamic Acid */
                        /* F.................Phe.................Phenylalanine */
                        /* G.................Gly.................Glycine */
                        /* H.................His.................Histidine */
                        /* I.................Ile.................Isoleucine */
                        /* K.................Lys.................Lysine */
                        /* L.................Leu.................Leucine */
                        /* M.................Met.................Methionine */
                        /* N.................Asn.................Asparagine */
                        /* P.................Pro.................Proline */
                        /* Q.................Gln.................Glutamine */
                        /* R.................Arg.................Arginine */
                        /* S.................Ser.................Serine */
                        /* T.................Thr.................Threonine */
                        /* V.................Val.................Valine */
                        /* W.................Trp.................Tryptophan */
                        /* X.................Xaa.................Any amino acid */
                        /* Y.................Tyr.................Tyrosine */
                        /* Z.................Glx.................Glutamine or Glutamic acid */
                        switch(toupper(sequence->seq[j])){
                        case 'A':
                                sequence->seq[j] = 0;
                                break;
                        case 'B'://D or N - 2 or 11
                                r = random_int_zero_to_x(1);
                                if(r == 0){
                                        sequence->seq[j] = 2;
                                }else{
                                        sequence->seq[j] = 11;
                                }
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
                        case 'X':
                                sequence->seq[j] = random_int_zero_to_x(19);
                                break;
                        case 'Z':
                                r = random_int_zero_to_x(1);
                                if(r == 0){
                                        sequence->seq[j] = 3;
                                }else{
                                        sequence->seq[j] = 13;
                                }
                                break;
                        default:
                                WARNING_MSG("Non ACGTN letter in sequence:%d %c.",i,sequence->seq[j]);
                                sequence->seq[j] = random_int_zero_to_x(19);
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
        sequence->score = 1.0f;
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
        ASSERT(sb != NULL, "No sequence buffer");
        if(sb->L == ALPHABET_DNA){
                RUN(translate_internal_to_DNA(sb));
        }
        if(sb->L == ALPHABET_PROTEIN){
                RUN(translate_internal_to_PROTEIN(sb));
        }
        /*
        for(i = 0; i < sb->num_seq;i++){
                sequence = sb->sequences[i];
                for(j = 0; j < sequence->seq_len;j++){
                        fprintf(stdout,"%4c",sequence->seq[j]);
                }
                fprintf(stdout,"\n");
                for(j = 0; j < sequence->seq_len;j++){
                        fprintf(stdout,"%4d",sequence->label[j]);
                }
                fprintf(stdout,"\n");
                for(j = 0; j < sequence->seq_len;j++){
                        fprintf(stdout,"%4.0f",scaledprob2prob(sequence->u[j])*100.0f);
                }
                fprintf(stdout,"\n");

        }*/

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
                if(sb->background){
                        MFREE(sb->background);
                }
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
                        ASSERT(a->sequences[i]->label[j]  == b->sequences[i]->label[j], "Labels differ");
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
        RUN(add_reverse_complement_sequences_to_buffer(iseq));
        RUN(write_ihmm_sequences(iseq,"test_dna.lfa","generated by iseq_ITEST"));



        RUNP(iseq_b  = load_ihmm_sequences("test_dna.lfa"));

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

        RUN(random_label_ihmm_sequences(iseq, 123,0.3));
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
