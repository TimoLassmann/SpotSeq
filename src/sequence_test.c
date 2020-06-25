#include "tldevel.h"

#include "sequence_alloc.h"
#include "sequence_io.h"
#include "sequence_struct.h"
int main(int argc, char *argv[])
{
        struct seq_buffer* sb = NULL;
        int i;


        RUN(read_sequences_file(&sb, argv[1]));

        for(i = 0;i < sb->num_seq;i++){
                fprintf(stdout,"%s\n%s\n", sb->sequences[i]->name, sb->sequences[i]->seq);
        }

        free_ihmm_sequences(sb);


        return EXIT_SUCCESS;
ERROR:

        return EXIT_FAILURE;
}
