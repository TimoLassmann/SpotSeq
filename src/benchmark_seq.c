

#include <inttypes.h>
#include "benchmark_seq.h"
#include "tldevel.h"

struct seq_buffer* alloc_seq_buffer(int num_seq)
{
        struct seq_buffer* sb  = NULL;
        struct sequence* sequence = NULL;
        int i;

        ASSERT(num_seq > 0, "No parameters.");

        MMALLOC(sb,sizeof(struct seq_buffer));
        sb->malloc_num = num_seq;
        sb->num_seq = 0;
        sb->seqs = NULL;
        MMALLOC(sb->seqs, sizeof(struct ihmm_sequence*) *sb->malloc_num );
        for(i = 0; i < sb->malloc_num;i++){
                sb->seqs[i] = NULL;
                sequence = NULL;
                MMALLOC(sequence,sizeof(struct sequence));
                sequence->seq = NULL;
                sequence->name = NULL;
                sequence->malloc_len = 128;
                sequence->seq_len = 0;
                MMALLOC(sequence->seq, sizeof(uint8_t) * sequence->malloc_len);
                MMALLOC(sequence->name, sizeof(char) * 256);
                sb->seqs[i] = sequence;
        }
        return sb;

        /*sequence->seq[sequence->seq_len] = line[i];
        sequence->seq_len++;
        if(sequence->seq_len == sequence->malloc_len){
                sequence->malloc_len = sequence->malloc_len << 1;
                MREALLOC(sequence->seq, sizeof(uint8_t) *sequence->malloc_len);
                }*/

ERROR:
        free_sb(sb);
        return NULL;
}

int write_sequences_to_file(struct seq_buffer* sb,char* filename)
{
        FILE* f_ptr = NULL;
        int i;
        ASSERT(sb != NULL," no sequences.");
        ASSERT(filename != NULL ," No filename.");


        RUNP(f_ptr = fopen(filename, "w"));

        for(i = 0; i < sb->num_seq;i++){
                fprintf(f_ptr,">Seq_%d\n%s\n",i+1,sb->seqs[i]->seq);
        }
        fclose(f_ptr);
        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        return FAIL;

}

int reset_sb(struct seq_buffer* sb)
{
        int i;
        sb->num_seq = 0;
        for(i = 0; i < sb->malloc_num;i++){
                sb->seqs[i]->seq_len = 0;
        }
        return OK;
}

void free_sb(struct seq_buffer* sb)
{
        int i;
        struct sequence* seq = NULL;
        if(sb){
                for(i =0; i < sb->malloc_num;i++){
                        seq = sb->seqs[i];
                        MFREE(seq->seq);
                        MFREE(seq->name);
                        MFREE(seq);
                }
                MFREE(sb->seqs);
                MFREE(sb);
        }
}
