#include "emit_random.h"

#include <getopt.h>

static int print_help(char **argv);


int main(const int argc, char * argv[])
{
        struct seq_buffer* sb = NULL;
        struct seq_buffer* sb_in = NULL;
        char* in = NULL;
        char* out = NULL;
        int c;
        int num_seq = 10000;

        print_program_header(argv, "Emit random sequences.");


        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"out",required_argument,0,'o'},

                        {"num",required_argument,0,'n'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"m:o:n:h",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 'm':
                        in = optarg;
                        break;
                case 'n':
                        num_seq = atoi(optarg);
                        break;
                case 'o':
                        out = optarg;
                        break;

                case 'h':
                        RUN(print_help(argv));
                        exit(EXIT_SUCCESS);
                        break;
                default:
                        ERROR_MSG("not recognized");
                        break;
                }
        }
        if(!in){
                RUN(print_help(argv));
                ERROR_MSG("No input file! use -m <blah.h5>");
        }
        if(!my_file_exists(in)){
                ERROR_MSG("The file <%s> does not exist.",in);
        }
        if(!out){
                RUN(print_help(argv));
                ERROR_MSG("No output file! use -o <blah.fa>");
        }
        if(num_seq <= 0){
                RUN(print_help(argv));
                ERROR_MSG("No sequences requested...! use -n XXXX ");
        }

        RUNP(sb_in = get_sequences_from_hdf5_model(in, IHMM_SEQ_READ_ONLY_SEQ));
        RUNP(sb = emit_sequences_from_random_model(sb_in, num_seq));

        RUN( write_ihmm_sequences_fasta(sb, out));
        free_ihmm_sequences(sb);
        free_ihmm_sequences(sb_in);
        return EXIT_SUCCESS;
ERROR:
        free_ihmm_sequences(sb);
        free_ihmm_sequences(sb_in);
        return EXIT_FAILURE;
}

int print_help(char **argv)
{
        const char usage[] = " -m <h5 model> -out <fasta output file> -n <number of sequences>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);

        return OK;
}
