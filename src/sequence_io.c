
#include "tldevel.h"
#include "tlseqio.h"

#include "sequence_struct.h"

#include "sequence_alloc.h"

#define MAX_SEQ_READ 1000000

#define SEQUENCE_IO_IMPORT
#include "sequence_io.h"


int read_sequences_file(struct seq_buffer** seq_buf,char* filename )
{
        struct seq_buffer* sb = NULL; /* structure used in seqer */
        struct file_handler* f_hand = NULL;
        struct tl_seq_buffer* tlsb = NULL; /* structure used for genetic fasta/ fastq reader  */

        int i,j;
        RUN(open_fasta_fastq_file(&f_hand, filename, TLSEQIO_READ));

        RUN(read_fasta_fastq_file(f_hand, &tlsb,MAX_SEQ_READ));
        if(tlsb->num_seq == MAX_SEQ_READ){
                WARNING_MSG("Currently seqer only works with the first %d sequences in an input file.",MAX_SEQ_READ);
        }
        if(tlsb->num_seq == 0){
                ERROR_MSG("No sequences found in file: %s.",filename);
        }

        RUN(detect_format(tlsb));

        RUN(alloc_sequence_buffer(&sb, tlsb->num_seq));

        for(i = 0; i < tlsb->num_seq ;i++){
                if(sb->sequences[i]->malloc_len < tlsb->sequences[i]->len){
                        RUN(realloc_ihmm_seq(sb->sequences[i], tlsb->sequences[i]->len));
                }
                snprintf(sb->sequences[i]->name, MAX_SEQUENCE_NAME_LEN, "%s", tlsb->sequences[i]->name);


                snprintf(sb->sequences[i]->seq , sb->sequences[i]->malloc_len, "%s", tlsb->sequences[i]->seq);
                sb->sequences[i]->seq_len = tlsb->sequences[i]->len;
        }

        sb->max_len = tlsb->max_len;
        sb->num_seq = tlsb->num_seq;

        free_tl_seq_buffer(tlsb);
        RUN(close_seq_file(&f_hand));

        *seq_buf = sb;

        return OK;
ERROR:
        return FAIL;
}
