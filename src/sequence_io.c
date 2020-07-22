#include <string.h>

#include "tlseqbuffer.h"
#include "tldevel.h"
#include "tlseqio.h"
#include "tlhdf5wrap.h"
#include "tlmisc.h"

#include "sequence_struct.h"
#include "sequence_alloc.h"

#define MAX_SEQ_READ 100000

#define SEQUENCE_IO_IMPORT
#include "sequence_io.h"

int read_sequences_file(struct seq_buffer** seq_buf,char* filename )
{
        struct seq_buffer* sb = NULL; /* structure used in seqer */
        struct file_handler* f_hand = NULL;
        struct tl_seq_buffer* tlsb = NULL; /* structure used for genetic fasta/ fastq reader  */
        //int i,j;
        //int printed;
        RUN(open_fasta_fastq_file(&f_hand, filename, TLSEQIO_READ));

        RUN(read_fasta_fastq_file(f_hand, &tlsb,MAX_SEQ_READ));
        if(tlsb->num_seq == MAX_SEQ_READ){
                WARNING_MSG("Currently seqer only works with the first %d sequences in an input file.",MAX_SEQ_READ);
        }
        if(tlsb->num_seq == 0){
                ERROR_MSG("No sequences found in file: %s.",filename);
        }

        RUN(detect_format(tlsb));

        RUN(convert_tl_seq_buf_into_ihmm_seq_buf(tlsb, &sb));

        free_tl_seq_buffer(tlsb);
        RUN(close_seq_file(&f_hand));
        *seq_buf = sb;

        return OK;
ERROR:
        return FAIL;
}

int convert_tl_seq_buf_into_ihmm_seq_buf(struct tl_seq_buffer* tlsb, struct seq_buffer** ret)
{
        struct seq_buffer* sb = NULL;
        int i,j;
        int printed;
        RUN(alloc_sequence_buffer(&sb, tlsb->num_seq));


        sb->max_len = 0;
        for(i = 0; i < tlsb->num_seq ;i++){
                if(sb->sequences[i]->malloc_len < tlsb->sequences[i]->len){
                        RUN(realloc_ihmm_seq(sb->sequences[i], tlsb->sequences[i]->len));
                }
                printed = snprintf(sb->sequences[i]->name, MAX_SEQUENCE_NAME_LEN, "%s", tlsb->sequences[i]->name);
                ASSERT(printed < MAX_SEQUENCE_NAME_LEN,"characters printed entirely fills buffer");

                sb->sequences[i]->seq_len = tlsb->sequences[i]->len;
                for(j = 0; j < tlsb->sequences[i]->len;j++){
                        sb->sequences[i]->seq[j] = tlsb->sequences[i]->seq[j];
                }

                if(sb->sequences[i]->seq_len > sb->max_len){
                        sb->max_len = sb->sequences[i]->seq_len;
                }
        }
        //LOG_MSG("MAXLEN: %d ", sb->max_len);
        sb->num_seq = tlsb->num_seq;
        sb->L = tlsb->L;

        *ret = sb;
        return OK;
ERROR:
        return FAIL;
}

int convert_ihmm_seq_buf_into_tl_seq_buf(struct seq_buffer* in,struct tl_seq_buffer** ret)
{

        struct tl_seq_buffer* sb = NULL;
        int i,j;
        int printed;

        RUN(alloc_tl_seq_buffer(&sb, in->num_seq));
        //RUN(alloc_sequence_buffer(&sb, tlsb->num_seq));


        sb->max_len;
        for(i = 0; i < in->num_seq ;i++){
                while(sb->sequences[i]->malloc_len <  in->sequences[i]->seq_len){
                        RUN(resize_tl_seq(sb->sequences[i]));
                        //RUN(realloc_ihmm_seq(sb->sequences[i], tlsb->sequences[i]->len));
                }
                printed = snprintf(sb->sequences[i]->name , TL_SEQ_MAX_NAME_LEN, "%s", in->sequences[i]->name );
                ASSERT(printed < TL_SEQ_MAX_NAME_LEN,"characters printed entirely fills buffer");

                sb->sequences[i]->len = in->sequences[i]->seq_len;
                for(j = 0; j < in->sequences[i]->seq_len ;j++){
                        sb->sequences[i]->seq[j] = in->sequences[i]->seq[j];
                }

                if(sb->sequences[i]->len > sb->max_len){
                        sb->max_len = sb->sequences[i]->len;
                }
        }
        //LOG_MSG("MAXLEN: %d ", sb->max_len);
        sb->num_seq = in->num_seq;
        sb->L = in->L;
        *ret = sb;
        return OK;
ERROR:
        return FAIL;
}





int add_sequences_to_hdf5_model(char* filename,struct seq_buffer* sb, int num_models)
{
        struct hdf5_data* hdf5_data = NULL;
        int i,j,c,len;
        int pos;

        char** name = NULL;
        char** seq = NULL;
        int** label = NULL;
        double** scores = NULL;
        int max_name_len;
        //int has_seq_info;

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

        RUN(galloc(&name, sb->num_seq, max_name_len));
        for(i = 0; i < sb->num_seq;i++){
                len = strlen(sb->sequences[i]->name);
                for (j = 0; j < len;j++){
                        name[i][j] = sb->sequences[i]->name[j];
                }
                for(j = len;j < max_name_len;j++){
                        name[i][j] = 0;
                }
        }

        /* make sequence matrix */

        //RUNP(seq = galloc(seq, sb->num_seq, sb->max_len, -1));
        RUN(galloc(&seq, sb->num_seq, sb->max_len));
        for(i = 0; i < sb->num_seq;i++){
                len = sb->sequences[i]->seq_len;
                for (j = 0; j < len;j++){
                        seq[i][j] = sb->sequences[i]->seq[j];
                }
                /* NEW  */
                for(j = len;j < sb->max_len;j++){
                        seq[i][j] = -1;
                }

        }

        /* make  label matrix */

        //RUNP(label = galloc(label, sb->num_seq, (sb->max_len+1)* num_models, -1));
        RUN(galloc(&label, sb->num_seq, (sb->max_len+1)* num_models));

        for(i = 0; i < sb->num_seq;i++){
                len = sb->sequences[i]->seq_len;
                pos = 0;
                for(c = 0; c < num_models;c++){
                        for (j = 0; j < sb->max_len+1;j++){
                                label[i][pos] = sb->sequences[i]->label_arr[c][j];
                                pos++;
                        }
                }
                for(c = pos; c < (sb->max_len+1)* num_models;c++){
                        label[i][c] = -1;
                }
        }

        /* Score matrix */

        //RUNP(scores = galloc(scores, sb->num_seq, num_models,0.0));
        RUN(galloc(&scores, sb->num_seq, num_models));

        for(i = 0; i < sb->num_seq;i++){
                len = sb->sequences[i]->seq_len;
                pos = 0;
                for(j = 0; j < num_models;j++){
                        scores[i][j] = sb->sequences[i]->score_arr[j];

                }
        }

        open_hdf5_file(&hdf5_data, filename);

        //hdf5_open_file(struct hdf5_data *hdf5_data)
        //RUNP(hdf5_data = hdf5_create());

        //hdf5_open_file(filename,hdf5_data);
        /*
        get_group_names(hdf5_data);

        has_seq_info = 0;

        for(i = 0; i < hdf5_data->grp_names->num_names;i++){
                if(strncmp("SequenceInformation", hdf5_data->grp_names->names[i],19) == 0){
                        has_seq_info = 1;
                }
        }
        if(!has_seq_info ){
                hdf5_create_group("SequenceInformation",hdf5_data);
        }else{
                hdf5_open_group("SequenceInformation",hdf5_data);
        }
        hdf5_data->num_attr = 0;
        */
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data, "/SequenceInformation", "Numseq", sb->num_seq));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data, "/SequenceInformation", "MaxLen", sb->max_len));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data, "/SequenceInformation", "MaxNameLen",max_name_len));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data, "/SequenceInformation", "Alphabet", sb->L));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data, "/SequenceInformation", "NumModels", num_models));

        RUN(HDFWRAP_WRITE_DATA(hdf5_data, "/SequenceInformation", "Names", name));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data, "/SequenceInformation", "Sequences", seq));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data, "/SequenceInformation", "Labels", label));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data, "/SequenceInformation", "CompetitiveScores", scores));
        //RUN(HDFWRAP_WRITE_DATA(hdf5_data, "/SequenceInformation", "Background", sb->background));
        RUN(close_hdf5_file(&hdf5_data));
        gfree(label);
        gfree(name);
        gfree(seq);
        return OK;
ERROR:
        if(hdf5_data){
                RUN(close_hdf5_file(&hdf5_data));

        }
        if(label){
                gfree(label);
        }
        if(name){
                gfree(name);
        }
        if(seq){
                gfree(seq);
        }
        return FAIL;
}

struct seq_buffer* get_sequences_from_hdf5_model(char* filename, int mode)
{
        struct hdf5_data* hdf5_data = NULL;
        struct seq_buffer* sb = NULL;
        char** name = NULL;
        char** seq = NULL;
        int** label = NULL;
        double** scores = NULL;
        //double* background;

        int num_seq;
        int max_len;
        int max_name_len;
        int local_L;
        int i,j,c;
        int num_models;
        int pos;
        ASSERT(filename != NULL, "No filename");
        ASSERT(my_file_exists(filename) != 0,"File %s does not exist.",filename);


        open_hdf5_file(&hdf5_data, filename);
        //hdf5_data = hdf5_create();
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, "/SequenceInformation", "Numseq", &num_seq));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, "/SequenceInformation", "MaxLen", &max_len));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, "/SequenceInformation", "MaxNameLen",&max_name_len));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, "/SequenceInformation", "Alphabet",&local_L));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, "/SequenceInformation", "NumModels", &num_models));

        ASSERT(num_seq != -1, "No numseq");
        ASSERT(max_len != -1, "No maxlen");
        ASSERT(max_name_len != -1, "No maxnamelen");
        ASSERT(local_L != -1,"No Alphabet");
        ASSERT(num_models > 0, "No models");


        RUN(HDFWRAP_READ_DATA(hdf5_data, "/SequenceInformation", "Names", &name));
        RUN(HDFWRAP_READ_DATA(hdf5_data, "/SequenceInformation", "Sequences", &seq));

        RUN(HDFWRAP_READ_DATA(hdf5_data, "/SequenceInformation", "CompetitiveScores", &scores));
        //RUN(HDFWRAP_READ_DATA(hdf5_data, "/SequenceInformation", "Background", &background));


        if(mode == IHMM_SEQ_READ_ALL){
                RUN(HDFWRAP_READ_DATA(hdf5_data, "/SequenceInformation", "Labels", &label));
        }else{
                label= NULL;
        }
        RUN(close_hdf5_file(&hdf5_data));
        //hdf5_close_file(hdf5_data);
        //hdf5_free(hdf5_data);

        MMALLOC(sb,sizeof(struct seq_buffer));

        sb->malloc_num = num_seq;
        sb->num_seq = num_seq;
        sb->org_num_seq = -1;
        sb->sequences = NULL;
        sb->max_len = max_len;
        sb->L = local_L;
        sb->alphabet = NULL;
        sb->num_state_arr = NULL;

        //sb->background = background;

        MMALLOC(sb->sequences, sizeof(struct ihmm_sequence*) *sb->malloc_num );
        for(i = 0; i < sb->num_seq;i++){
                sb->sequences[i] = NULL;
                RUN(alloc_ihmm_seq(&sb->sequences[i]));
                while(sb->sequences[i]->malloc_len <= sb->max_len){
                        RUN(realloc_ihmm_seq(sb->sequences[i],sb->max_len));
                }
        }
        if(mode == IHMM_SEQ_READ_ALL){
                RUN(add_multi_model_label_and_u(sb, num_models));
                /* copy stuff over */
                for(i = 0; i < sb->num_seq;i++){
                        pos = 0;
                        for(c = 0; c < num_models;c++){
                                for (j = 0; j < sb->max_len+1;j++){
                                        sb->sequences[i]->label_arr[c][j] = label[i][pos];
                                        pos++;
                                }
                                sb->sequences[i]->score_arr[c] = scores[i][c];
                        }
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
                }
                sb->sequences[i]->seq_len = j;
        }


        if(mode == IHMM_SEQ_READ_ALL){
                gfree(label);
        }
        gfree(scores);
        gfree(name);
        gfree(seq);
        return sb;
ERROR:
        if(hdf5_data){
                RUN(close_hdf5_file(&hdf5_data));
        }
        if(label){
                gfree(label);
        }
        if(name){
                gfree(name);
        }
        if(seq){
                gfree(seq);
        }

        return NULL;
}
