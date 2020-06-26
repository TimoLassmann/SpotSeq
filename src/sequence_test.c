#include "tldevel.h"
#include "tlrng.h"

#include "sequence_alloc.h"
#include "sequence_io.h"
#include "sequence_struct.h"
#include "sequence_prep.h"


int main(int argc, char *argv[])
{
        struct rng_state* rng = NULL;
        struct seq_buffer* sb = NULL;
        int i;
        int l;
        const int p_limit = 100;

        RUNP(rng = init_rng(0));

        RUN(read_sequences_file(&sb, argv[1]));


        l = MACRO_MIN(p_limit, sb->num_seq);
        for(i = 0;i < l;i++){
                fprintf(stdout,"%s\n%s\n", sb->sequences[i]->name, sb->sequences[i]->seq);
        }





        RUN(prep_sequences(sb, rng, 1,0,0.0));

        free_rng(rng);
        free_ihmm_sequences(sb);


        return EXIT_SUCCESS;
ERROR:

        return EXIT_FAILURE;
}
